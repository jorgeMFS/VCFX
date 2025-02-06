#include "VCFX_hwe_tester.h"
#include <getopt.h>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cstdlib>

// ---------------------------------------------------
// run
// ---------------------------------------------------
int VCFXHWETester::run(int argc, char* argv[]) {
    // No real arguments for now, just handle --help
    int opt;
    bool showHelp = false;
    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };
    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
            default:
                showHelp = true;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Perform HWE on stdin
    performHWE(std::cin);
    return 0;
}

// ---------------------------------------------------
// displayHelp
// ---------------------------------------------------
void VCFXHWETester::displayHelp() {
    std::cout 
        << "VCFX_hwe_tester: Perform Hardy-Weinberg Equilibrium (HWE) tests on a biallelic VCF.\n\n"
        << "Usage:\n"
        << "  VCFX_hwe_tester [options] < input.vcf > output\n\n"
        << "Description:\n"
        << "  Reads each variant line, ignoring multi-allelic calls. For biallelic lines,\n"
        << "  collects genotypes as 0/0, 0/1, 1/1, then uses an exact test to produce\n"
        << "  a p-value for HWE.\n\n"
        << "Example:\n"
        << "  VCFX_hwe_tester < input.vcf > results.txt\n";
}

// ---------------------------------------------------
// isBiallelic: returns true if alt has no commas
// ---------------------------------------------------
bool VCFXHWETester::isBiallelic(const std::string &alt) {
    // If ALT has a comma, it implies multiple alt alleles. We'll skip
    if (alt.find(',') != std::string::npos) {
        return false;
    }
    return true;
}

// ---------------------------------------------------
// parseGenotypes
//   we interpret 0/1,1/0 => het, 0/0 => homRef, 1/1 => homAlt
//   skip or ignore any other combos
// ---------------------------------------------------
bool VCFXHWETester::parseGenotypes(const std::vector<std::string>& genotypes,
                                   int& homRef, int& het, int& homAlt) {
    homRef=0;
    het=0;
    homAlt=0;
    for (auto &gt : genotypes) {
        if (gt.empty() || gt=="." || gt=="./." || gt==".|.") {
            // missing => skip
            continue;
        }
        // unify '|' with '/'
        std::string g = gt;
        for (char &c : g) {
            if (c=='|') c='/';
        }
        // split
        std::size_t delim = g.find('/');
        if (delim==std::string::npos) {
            // not diploid => skip
            continue;
        }
        std::string a1 = g.substr(0, delim);
        std::string a2 = g.substr(delim+1);
        if (a1=="0" && a2=="0") {
            homRef++;
        } else if ((a1=="0" && a2=="1") || (a1=="1" && a2=="0")) {
            het++;
        } else if (a1=="1" && a2=="1") {
            homAlt++;
        } else {
            // e.g. a2>1 => multi-allelic, or unknown => skip
            // we might do an error or skip
            return false; 
        }
    }
    return true;
}

// ---------------------------------------------------
// genotypeProbability
//   core exact test approach
// ---------------------------------------------------
double VCFXHWETester::genotypeProbability(int homRef, int het, int homAlt) {
    int N = homRef + het + homAlt;
    int x = 2*homAlt + het; // # alt alleles
    int y = 2*homRef + het; // # ref
    // We do an enumerative approach
    // log factorial caching
    static const int MAXN = 20000; // arbitrary
    static std::vector<double> logFac; 
    static int cachedN=0;
    if (N*2>cachedN) {
        // rebuild up to 2*N
        logFac.resize(2*N+1, 0.0);
        double cum=0.0;
        for (int i=1; i<=2*N; i++) {
            cum += std::log((double)i);
            logFac[i]=cum;
        }
        cachedN=2*N;
    }
    auto logFactorial = [&](int n)->double {
        if (n<=0) return 0.0;
        if (n>cachedN) return 0.0; 
        return logFac[n];
    };

    // The function to compute log prob of a genotype distribution "a"
    // where a is # of het
    auto logProb = [&](int a)->double {
        // homRef= (y-a)/2
        // homAlt= (x-a)/2
        if ((x-a)<0 || (y-a)<0) return -INFINITY;
        if (((x-a)%2)!=0 || ((y-a)%2)!=0) return -INFINITY;
        int hr = (y-a)/2;
        int ha = (x-a)/2;
        if (hr<0 || ha<0) return -INFINITY;
        // let's define M = hr+ha+a = N
        // log multinomial = log(N!) - log(hr!) - log(a!) - log(ha!)
        double lcoef = logFactorial(N) - (logFactorial(hr)+logFactorial(ha)+logFactorial(a));
        // p= y/(x+y), q= x/(x+y)
        double p= (double)y/(double)(x+y);
        double q= (double)x/(double)(x+y);
        if (p<=0.0 || q<=0.0) return -INFINITY;
        double logTerm= hr*2*std::log(p)+ a*std::log(2.0*p*q)+ ha*2*std::log(q);
        return lcoef + logTerm;
    };

    int observedHet = het;
    double logObs = logProb(observedHet);
    if (logObs==-INFINITY) {
        return 1.0; 
    }

    // sum over all possible "a"
    // find min log for stability
    double minLog = logObs;
    int maxA = std::min(x, y);
    std::vector<double> logVals(maxA+1, -INFINITY);
    logVals[observedHet] = logObs;
    if (logObs<minLog) minLog= logObs;

    // fill others
    for (int a=0; a<=maxA; a++) {
        if (a==observedHet) continue;
        double lp = logProb(a);
        logVals[a] = lp;
        if (lp<minLog && lp!=-INFINITY) {
            minLog= lp;
        }
    }

    // exponentiate relative
    double sum=0.0, obsExp=0.0;
    for (int a=0; a<=maxA; a++) {
        if (logVals[a]==-INFINITY) continue;
        double rel= logVals[a]- minLog;
        double e= std::exp(rel);
        sum+= e;
        if (a==observedHet) obsExp= e;
    }

    // pvalue => sum of e <= obsExp
    // but the typical approach is to do a "two-sided" test or "mid-p" approach
    // we'll do a basic approach: all states with prob <= probObs => pval
    double probObs= obsExp/sum;
    double pVal=0.0;
    for (int a=0; a<=maxA; a++) {
        if (logVals[a]==-INFINITY) continue;
        double rel= logVals[a]- minLog;
        double e= std::exp(rel);
        double thisProb= e/sum;
        if (thisProb <= probObs+1e-12) {
            pVal+= thisProb;
        }
    }
    if (pVal>1.0) pVal=1.0;
    if (pVal<0.0) pVal=0.0;
    return pVal;
}

// ---------------------------------------------------
// calculateHWE
//   if multi-allelic => skip
//   we do the exact test
// ---------------------------------------------------
double VCFXHWETester::calculateHWE(int homRef, int het, int homAlt) {
    // total # samples
    int N= homRef+ het+ homAlt;
    if (N<1) {
        return 1.0;
    }
    // for biallelic
    int x= 2*homAlt + het;
    int y= 2*homRef + het;
    if (x+y != 2*N) {
        // not biallelic or doesn't match
        return 1.0;
    }
    // do exact test
    return genotypeProbability(homRef, het, homAlt);
}

// ---------------------------------------------------
// performHWE
//   we skip multi-allelic lines and parse GT=0/0,0/1,1/1
//   print pvalue
// ---------------------------------------------------
void VCFXHWETester::performHWE(std::istream& in) {
    std::string line;
    // print header
    std::cout << "CHROM\tPOS\tID\tREF\tALT\tHWE_pvalue\n";

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0]=='#') {
            // skip original headers
            continue;
        }

        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string f;
            while(std::getline(ss,f,'\t')) {
                fields.push_back(f);
            }
        }
        if (fields.size()<10) {
            // not enough columns
            continue;
        }
        std::string chrom= fields[0];
        std::string pos  = fields[1];
        std::string id   = fields[2];
        std::string ref  = fields[3];
        std::string alt  = fields[4];
        // check biallelic
        if (!isBiallelic(alt)) {
            // skip
            continue;
        }

        // find GT index
        std::string format= fields[8];
        std::vector<std::string> fmts;
        {
            std::stringstream fs(format);
            std::string x;
            while(std::getline(fs,x,':')) {
                fmts.push_back(x);
            }
        }
        int gt_index= -1;
        for (size_t i=0; i<fmts.size(); i++) {
            if (fmts[i]=="GT") {
                gt_index= i;
                break;
            }
        }
        if (gt_index<0) {
            continue;
        }
        // gather genotypes
        std::vector<std::string> gts;
        for (size_t s=9; s<fields.size(); s++) {
            std::stringstream sampSS(fields[s]);
            std::vector<std::string> sampToks;
            {
                std::string xx;
                while(std::getline(sampSS,xx,':')) {
                    sampToks.push_back(xx);
                }
            }
            if ((size_t)gt_index< sampToks.size()) {
                gts.push_back(sampToks[gt_index]);
            } else {
                gts.push_back(".");
            }
        }

        int hr=0, h=0, ha=0;
        bool ok= parseGenotypes(gts, hr, h, ha);
        if (!ok) {
            // skip or do partial
            continue;
        }
        double pVal= calculateHWE(hr,h,ha);
        std::cout << chrom << "\t" << pos << "\t" << id << "\t" 
                  << ref   << "\t" << alt << "\t"
                  << std::fixed << std::setprecision(6) << pVal << "\n";
    }
}

// ---------------------------------------------------
// isBiallelic
// ---------------------------------------------------
bool VCFXHWETester::isBiallelic(const std::string &alt) {
    // if alt has a comma => multi-allelic => false
    if (alt.find(',')!=std::string::npos) {
        return false;
    }
    return true;
}

// ---------------------------------------------------
// main
// ---------------------------------------------------
int main(int argc, char* argv[]) {
    VCFXHWETester tester;
    return tester.run(argc, argv);
}
