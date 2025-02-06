#include "VCFX_inbreeding_calculator.h"
#include <getopt.h>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <iomanip>

// -------------------------------------------------------------------------
// displayHelp
// -------------------------------------------------------------------------
void VCFXInbreedingCalculator::displayHelp() {
    std::cout 
        << "VCFX_inbreeding_calculator: Compute individual inbreeding coefficients (F)\n"
        << "based on biallelic sites in a VCF.\n\n"
        << "Usage:\n"
        << "  VCFX_inbreeding_calculator [options] < input.vcf > output.txt\n\n"
        << "Description:\n"
        << "  Reads a VCF. For each biallelic site (no commas in ALT), collects each\n"
        << "  sample's genotype as 0/0 =>0, 0/1 =>1, 1/1 =>2. Then for each sample, we\n"
        << "  exclude that sample's allele(s) when computing the alt allele frequency.\n"
        << "  The expected heterozygosity for that site is 2 p(1-p). Summing across all\n"
        << "  included sites => sumExp. Observed number of heterozygous sites => obsHet.\n"
        << "  We define F = 1 - obsHet / sumExp. If no sites, F=NA.\n\n"
        << "Example:\n"
        << "  VCFX_inbreeding_calculator < input.vcf > inbreeding.txt\n";
}

// -------------------------------------------------------------------------
// parseGenotype: parse "0/1" or "1|0". We only accept 0 or 1 as alleles.
//   Return -1 if missing or multi-allelic. Return 0 => homRef, 1 => het, 2 => homAlt
// -------------------------------------------------------------------------
int VCFXInbreedingCalculator::parseGenotype(const std::string& s) {
    if (s.empty() || s==".") {
        return -1;
    }
    // unify '|' => '/'
    std::string g(s);
    for (char &c : g) {
        if (c=='|') c='/';
    }
    // split on '/'
    size_t slash= g.find('/');
    if (slash==std::string::npos) {
        return -1; // not diploid
    }
    std::string a1= g.substr(0, slash);
    std::string a2= g.substr(slash+1);
    if (a1=="." || a2=="." || a1.empty() || a2.empty()) {
        return -1; // missing
    }
    // must parse them as int. If either >1 => skip
    int i1=0, i2=0;
    try {
        i1= std::stoi(a1);
        i2= std::stoi(a2);
    } catch(...) {
        return -1;
    }
    if (i1<0 || i2<0) {
        return -1; // negative not valid
    }
    // if either allele >=2 => multi-allelic => skip
    if (i1>1 || i2>1) {
        return -1;
    }
    // so we have only 0 or 1
    if (i1==i2) {
        // if both 0 => code=0, if both 1 => code=2
        return (i1==0) ? 0 : 2;
    }
    // otherwise => 1
    return 1;
}

// -------------------------------------------------------------------------
// isBiallelic: if ALT has a comma => multi-ALT => skip
// -------------------------------------------------------------------------
bool VCFXInbreedingCalculator::isBiallelic(const std::string &alt) {
    if (alt.find(',')!=std::string::npos) {
        return false;
    }
    return true;
}

// -------------------------------------------------------------------------
// calculateInbreeding
//   1) read header => sample names
//   2) read variants => if biallelic, parse genotype => store code
//   3) after reading all, for each sample, do a loop over variants => compute F
// -------------------------------------------------------------------------
void VCFXInbreedingCalculator::calculateInbreeding(std::istream& in,
                                                   std::ostream& out)
{
    std::string line;
    bool foundChrom= false;
    std::vector<std::string> sampleNames;
    // We'll store all variants in a vector
    std::vector<InbreedingVariant> variants;
    int numSamples=0;

    // read header lines
    while (true) {
        auto pos= in.tellg();
        if(!std::getline(in, line)) {
            break; 
        }
        if (line.empty()) continue;
        if (line[0]=='#') {
            // check if #CHROM => parse sample columns
            if(!foundChrom && line.rfind("#CHROM",0)==0) {
                foundChrom= true;
                // parse
                std::stringstream ss(line);
                std::string t;
                std::vector<std::string> tokens;
                while(std::getline(ss,t,'\t')) {
                    tokens.push_back(t);
                }
                // from col=9 onward => samples
                for (size_t c=9; c<tokens.size(); c++){
                    sampleNames.push_back(tokens[c]);
                }
                numSamples= sampleNames.size();
            }
            // skip printing or storing header lines
            continue;
        } else {
            // not a header => revert and break
            in.seekg(pos);
            break;
        }
    }
    if(!foundChrom){
        std::cerr<<"Error: No #CHROM line found.\n";
        return;
    }

    // read data lines
    while(std::getline(in,line)) {
        if(line.empty()|| line[0]=='#') {
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
        if(fields.size()<10) {
            // skip
            continue;
        }
        std::string chrom= fields[0];
        std::string posStr= fields[1];
        std::string id= fields[2];
        std::string ref= fields[3];
        std::string alt= fields[4];
        if(!isBiallelic(alt)) {
            // skip multi-ALT
            continue;
        }
        int posVal=0;
        try {
            posVal= std::stoi(posStr);
        } catch(...) {
            continue;
        }

        // The format col= fields[8], we won't parse subfields, we assume GT is first or something
        // We'll just parse the sample columns directly for genotype
        // i.e. fields[9..8+numSamples]
        InbreedingVariant v;
        v.chrom= chrom;
        v.pos= posVal;
        v.genotypeCodes.resize(numSamples, -1);
        int sampleCol=9;
        for(int s=0; s<numSamples; s++){
            if((size_t)(sampleCol+s)>= fields.size()) {
                // missing
                break;
            }
            int code= parseGenotype(fields[sampleCol+s]);
            v.genotypeCodes[s]= code;
        }
        variants.push_back(v);
    }

    if(variants.empty()){
        std::cerr<<"No biallelic variants found.\n";
        // we can output an empty result
        return;
    }

    // Now compute inbreeding for each sample
    // We'll do a 2 pass approach: in pass 1, for each variant we compute total # alt alleles among
    // all samples that have code>=0. We'll store that in an int altCount
    // We'll also store how many samples have code>=0 => count *2 => totalAlleles

    // store altCount, nGenotypes for each variant
    std::vector<int> altCounts(variants.size(),0);
    std::vector<int> validSampleCounts(variants.size(),0);

    // pass 1
    for(size_t vIdx=0; vIdx<variants.size(); vIdx++){
        int altC=0;
        int nG=0;
        const auto & codes= variants[vIdx].genotypeCodes;
        for(int s=0; s<numSamples; s++){
            int c= codes[s];
            if(c<0) continue; // skip
            // c in [0,1,2]
            // alt alleles => c
            // 0 => no alt, 1 => 1 alt, 2 => 2 alt
            altC+= c;
            nG++;
        }
        altCounts[vIdx]= altC;
        validSampleCounts[vIdx]= nG;
    }

    // pass 2: for each sample, we do sum of obsHet and sumExp
    std::vector<double> obsHet (numSamples, 0.0);
    std::vector<double> sumExp(numSamples, 0.0);
    std::vector<int> usedVariantCount(numSamples,0);

    for(size_t vIdx=0; vIdx<variants.size(); vIdx++){
        const auto & codes= variants[vIdx].genotypeCodes;
        int altC= altCounts[vIdx];
        int nG= validSampleCounts[vIdx];
        if(nG<2) {
            // not enough samples
            continue;
        }
        // for each sample with code>=0 => compute freq excluding sample
        // altEx= altC - code, nEx= nG-1 => p= altEx/(2*nEx)
        for(int s=0; s<numSamples; s++){
            int c= codes[s];
            if(c<0) {
                // skip
                continue;
            }
            int altEx= altC - c;
            int nEx= nG -1;
            if(nEx<1) {
                // can't do freq
                continue;
            }
            double p= (double)altEx/(2.0*nEx);
            // expected het => 2 p(1-p)
            double eHet= 2.0 * p * (1.0 - p);
            sumExp[s]+= eHet;
            if(c==1) {
                obsHet[s]+=1.0;
            }
            usedVariantCount[s]++;
        }
    }

    // now compute F => 1 - obsHet/sumExp
    // if sumExp=0 => "NA"
    // output
    std::cout<<"Sample\tInbreedingCoefficient\n";
    for(int s=0; s<numSamples; s++){
        if(usedVariantCount[s]==0) {
            std::cout<< sampleNames[s] << "\tNA\n";
            continue;
        }
        double e= sumExp[s];
        if(e<=0.0) {
            std::cout<< sampleNames[s] << "\tNA\n";
            continue;
        }
        double f= 1.0 - (obsHet[s]/ e);
        std::cout<< sampleNames[s] << "\t"<< std::fixed << std::setprecision(6) << f <<"\n";
    }
}

// -------------------------------------------------------------------------
// run
// -------------------------------------------------------------------------
int VCFXInbreedingCalculator::run(int argc, char* argv[]) {
    // parse args for help
    int opt;
    bool showHelp= false;
    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };
    while(true){
        int c= getopt_long(argc, argv, "h", long_opts, NULL);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp= true;
                break;
            default:
                showHelp= true;
        }
    }
    if(showHelp) {
        displayHelp();
        return 0;
    }

    // do main
    calculateInbreeding(std::cin, std::cout);
    return 0;
}

int main(int argc, char* argv[]){
    VCFXInbreedingCalculator calc;
    return calc.run(argc, argv);
}
