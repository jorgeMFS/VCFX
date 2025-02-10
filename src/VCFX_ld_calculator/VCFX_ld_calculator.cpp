#include "VCFX_ld_calculator.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <iomanip>


// ----------------------------------------------------------------------
// Display Help
// ----------------------------------------------------------------------
void VCFXLDCalculator::displayHelp() {
    std::cout
        << "VCFX_ld_calculator: Calculate pairwise LD (r^2) for variants in a VCF region.\n\n"
        << "Usage:\n"
        << "  VCFX_ld_calculator [options] < input.vcf\n\n"
        << "Options:\n"
        << "  --region <chr:start-end>  Only compute LD for variants in [start, end] on 'chr'.\n"
        << "  --help, -h                Show this help message.\n\n"
        << "Description:\n"
        << "  Reads a VCF from stdin (uncompressed) and collects diploid genotypes as code:\n"
        << "     0 => homRef (0/0), 1 => het (0/1), 2 => homAlt(1/1), -1 => missing/other.\n"
        << "  Then for each pair of variants in the region, compute pairwise r^2 ignoring samples with missing genotype.\n"
        << "  Output an MxM matrix of r^2 along with variant identifiers.\n\n"
        << "Example:\n"
        << "  VCFX_ld_calculator --region chr1:10000-20000 < input.vcf > ld_matrix.txt\n";
}

// ----------------------------------------------------------------------
// parseRegion( "chr1:10000-20000" ) => regionChrom="chr1", regionStart=10000, regionEnd=20000
// returns false if can't parse
// ----------------------------------------------------------------------
bool VCFXLDCalculator::parseRegion(const std::string &regionStr,
                                   std::string &regionChrom,
                                   int &regionStart,
                                   int &regionEnd)
{
    // find ':'
    auto colonPos= regionStr.find(':');
    if(colonPos==std::string::npos) {
        return false;
    }
    regionChrom = regionStr.substr(0, colonPos);
    auto dashPos= regionStr.find('-', colonPos+1);
    if(dashPos==std::string::npos) {
        return false;
    }
    std::string startStr= regionStr.substr(colonPos+1, dashPos-(colonPos+1));
    std::string endStr= regionStr.substr(dashPos+1);
    try {
        regionStart= std::stoi(startStr);
        regionEnd=   std::stoi(endStr);
    } catch(...) {
        return false;
    }
    if(regionStart> regionEnd) {
        return false;
    }
    return true;
}

// ----------------------------------------------------------------------
// parse a genotype string like "0/1" or "1|0"
//   => 0 => 0/0
//   => 1 => 0/1
//   => 2 => 1/1
//   => -1 => missing or multi-allelic
// ----------------------------------------------------------------------
int VCFXLDCalculator::parseGenotype(const std::string &s) {
    if(s.empty() || s=="." || s=="./." || s==".|.") {
        return -1;
    }
    // unify '|' => '/'
    std::string g(s);
    for (char &c : g) {
        if(c=='|') c='/';
    }
    // split
    auto slashPos= g.find('/');
    if(slashPos==std::string::npos) {
        return -1; // not diploid
    }
    auto a1= g.substr(0,slashPos);
    auto a2= g.substr(slashPos+1);
    if(a1.empty()||a2.empty()) return -1;
    if(a1=="."|| a2==".") return -1;
    // parse int
    int i1=0, i2=0;
    try {
        i1= std::stoi(a1);
        i2= std::stoi(a2);
    } catch(...) {
        return -1;
    }
    if(i1<0 || i2<0) return -1;
    if(i1>1 || i2>1) {
        // multi-allelic
        return -1;
    }
    if(i1==i2) {
        // 0 => 0/0 => code=0, 1 => 1/1 => code=2
        if(i1==0) return 0; 
        else return 2;
    } else {
        return 1;
    }
}

// ----------------------------------------------------------------------
// compute r^2 between genotype arrays g1,g2 ignoring -1
// standard formula
// Let n = count of samples where both g1,g2 != -1
// X[i]=g1[i], Y[i]=g2[i]
// meanX= average of X, meanY= average of Y
// cov= average(XY)-meanX*meanY
// varX= average(X^2)- meanX^2, varY likewise
// r= cov / sqrt(varX varY), r^2= r*r
// ----------------------------------------------------------------------
double VCFXLDCalculator::computeRsq(const std::vector<int> &g1,
                                    const std::vector<int> &g2)
{
    if(g1.size()!= g2.size()) return 0.0; 
    int n=0;
    long sumX=0, sumY=0, sumXY=0, sumX2=0, sumY2=0;
    for(size_t i=0; i<g1.size(); i++){
        int x= g1[i];
        int y= g2[i];
        if(x<0||y<0) continue; // skip missing
        n++;
        sumX+= x;
        sumY+= y;
        sumXY+= x*y;
        sumX2+= (x*x);
        sumY2+= (y*y);
    }
    if(n<2) return 0.0;
    double meanX= (double)sumX/(double)n;
    double meanY= (double)sumY/(double)n;
    double cov= ((double)sumXY/(double)n) - (meanX* meanY);
    double varX= ((double)sumX2/(double)n) - (meanX* meanX);
    double varY= ((double)sumY2/(double)n) - (meanY* meanY);
    if(varX<=0.0 || varY<=0.0) return 0.0;
    double r= cov/ (std::sqrt(varX)* std::sqrt(varY));
    return r*r;
}

// ----------------------------------------------------------------------
// computeLD: read lines, skip #, parse sample columns
//   store variants within region
//   then produce MxM matrix of r^2
// ----------------------------------------------------------------------
void VCFXLDCalculator::computeLD(std::istream &in,
                                 std::ostream &out,
                                 const std::string &regionChrom,
                                 int regionStart,
                                 int regionEnd)
{
    bool foundChromHeader= false;
    std::vector<std::string> sampleNames;
    std::vector<LDVariant> variants;
    int numSamples=0;

    // parse header => get sample col
    // pass header lines as is
    std::string line;
    while(true) {
        auto pos= in.tellg();
        if(!std::getline(in,line)) break;
        if(line.empty()) {
            out << line << "\n";
            continue;
        }
        if(line[0]=='#') {
            out << line << "\n";
            if(!foundChromHeader && line.rfind("#CHROM",0)==0) {
                // parse sample
                foundChromHeader= true;
                std::stringstream ss(line);
                std::string tok;
                std::vector<std::string> tokens;
                while(std::getline(ss,tok,'\t')) {
                    tokens.push_back(tok);
                }
                // from col=9 onward => sample
                for(size_t c=9; c<tokens.size(); c++){
                    sampleNames.push_back(tokens[c]);
                }
                numSamples= sampleNames.size();
            }
            continue;
        } 
        // data line
        if(!foundChromHeader) {
            std::cerr<<"Error: encountered data line before #CHROM.\n";
            // we can break or skip
            break;
        }
        // parse
        // CHROM POS ID REF ALT QUAL FILTER INFO ...
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string f;
            while(std::getline(ss,f,'\t')){
                fields.push_back(f);
            }
        }
        if(fields.size()<10) {
            // not enough columns => pass line => skip
            out << line << "\n";
            continue;
        }
        std::string &chrom= fields[0];
        std::string &posStr= fields[1];
        int posVal=0;
        try {
            posVal= std::stoi(posStr);
        } catch(...) {
            // pass
            out << line <<"\n";
            continue;
        }
        // if regionChrom not empty => check if chrom matches
        // if regionChrom empty => we do everything
        if(!regionChrom.empty()) {
            if(chrom!= regionChrom) {
                // not in region => pass line
                out << line <<"\n";
                continue;
            }
            if(posVal< regionStart || posVal> regionEnd) {
                // out of range
                out << line <<"\n";
                continue;
            }
        }
        // we keep this variant => parse genotype codes
        LDVariant v;
        v.chrom= chrom;
        v.pos=   posVal;
        v.genotype.resize(numSamples, -1);
        // sample columns => fields[9..]
        int sampleCol=9;
        for(int s=0; s<numSamples; s++){
            if((size_t)(sampleCol+s)>= fields.size()) break;
            v.genotype[s] = parseGenotype(fields[sampleCol+s]);
        }
        variants.push_back(std::move(v));
        // pass line unchanged
        out << line << "\n";
    }

    // after reading all => we compute pairwise r^2
    // if we have M variants => produce an MxM matrix
    // we do: #LD_MATRIX_START, then we print M lines
    // first line => #VariantIndices ...
    // or we print a table with variant index as row, col + r^2
    // We'll do a simple text approach
    size_t M= variants.size();
    if(M<2) {
        out << "#LD_MATRIX_START\n";
        out << "No or only one variant in the region => no pairwise LD.\n";
        out << "#LD_MATRIX_END\n";
        return;
    }

    out << "#LD_MATRIX_START\n";
    // Print header row => e.g. row0col0 is "", row0col i => i. We'll also print "Chrom:Pos" as col
    // first line => tab, var1, var2, ...
    out << "Index/Var";
    for (size_t j=0; j<M; j++){
        out << "\t" << variants[j].chrom << ":" << variants[j].pos;
    }
    out << "\n";

    // each row => i => row header => i-chrom:pos => r2 vs all j
    for(size_t i=0; i<M; i++){
        out << variants[i].chrom << ":" << variants[i].pos;
        for(size_t j=0; j<M; j++){
            if(j< i) {
                // symmetrical => we can do the same as [j][i]
                // but let's just compute anyway or store
                double r2= computeRsq(variants[i].genotype, variants[j].genotype);
                out << "\t" << std::fixed << std::setprecision(4) << r2;
            } else if(i==j) {
                // r2 with self => 1.0
                out << "\t1.0000";
            } else {
                // i< j => compute
                double r2= computeRsq(variants[i].genotype, variants[j].genotype);
                out << "\t" << std::fixed << std::setprecision(4) << r2;
            }
        }
        out << "\n";
    }

    out << "#LD_MATRIX_END\n";
}

// ----------------------------------------------------------------------
// main run
// ----------------------------------------------------------------------
int VCFXLDCalculator::run(int argc, char* argv[]) {
    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"region", required_argument, 0, 'r'},
        {0,0,0,0}
    };
    bool showHelp=false;
    std::string regionStr;
    while(true){
        int c= getopt_long(argc, argv, "hr:", long_opts, NULL);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 'r':
                regionStr= optarg;
                break;
            default:
                showHelp=true;
        }
    }
    if(showHelp) {
        displayHelp();
        return 0;
    }

    // parse region if any
    std::string regionChrom;
    int regionStart=0, regionEnd=0;
    if(!regionStr.empty()) {
        if(!parseRegion(regionStr, regionChrom, regionStart, regionEnd)){
            std::cerr<<"Error parsing region '"<<regionStr<<"'. Use e.g. chr1:10000-20000\n";
            return 1;
        }
    }

    // do main
    computeLD(std::cin, std::cout, regionChrom, regionStart, regionEnd);
    return 0;
}

int main(int argc, char* argv[]) {
    VCFXLDCalculator calc;
    return calc.run(argc, argv);
}
