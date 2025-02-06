#include "VCFX_haplotype_phaser.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

int VCFXHaplotypePhaser::run(int argc, char* argv[]) {
    // parse arguments
    int opt;
    bool showHelp = false;
    double ldThreshold = 0.8;

    static struct option long_options[] = {
        {"help",        no_argument,       0, 'h'},
        {"ld-threshold", required_argument, 0, 'l'},
        {0,             0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hl:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'l':
                try {
                    ldThreshold = std::stod(optarg);
                } catch(...) {
                    std::cerr << "Error: invalid LD threshold.\n";
                    displayHelp();
                    return 1;
                }
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || ldThreshold < 0.0 || ldThreshold > 1.0) {
        displayHelp();
        return 1;
    }

    phaseHaplotypes(std::cin, std::cout, ldThreshold);
    return 0;
}

void VCFXHaplotypePhaser::displayHelp() {
    std::cout << "VCFX_haplotype_phaser: Group variants into blocks by naive LD threshold.\n\n"
              << "Usage:\n"
              << "  VCFX_haplotype_phaser [options] < input.vcf > blocks.txt\n\n"
              << "Options:\n"
              << "  -h, --help               Show this help message\n"
              << "  -l, --ld-threshold <val> r^2 threshold [0..1], default 0.8\n\n"
              << "Example:\n"
              << "  VCFX_haplotype_phaser --ld-threshold 0.9 < input.vcf > blocks.txt\n";
}

void VCFXHaplotypePhaser::phaseHaplotypes(std::istream& in, std::ostream& out, double ldThreshold) {
    std::string line;
    bool foundHeader = false;

    // We'll store the final list of variants
    std::vector<VariantData> variantList;
    std::vector<std::string> sampleNames;

    // read header lines until first non-# line
    while (true) {
        auto pos = in.tellg();
        if (!std::getline(in, line)) {
            // no data
            break;
        }
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            // if it's #CHROM, parse sample columns
            if (!foundHeader && line.rfind("#CHROM",0)==0) {
                // parse
                std::stringstream ss(line);
                std::string f;
                std::vector<std::string> tokens;
                while (std::getline(ss, f, '\t')) {
                    tokens.push_back(f);
                }
                // from col=9 onward => samples
                for (size_t c=9; c<tokens.size(); c++) {
                    sampleNames.push_back(tokens[c]);
                }
                foundHeader = true;
            }
            // print the header line out as well
            out << line << "\n";
        } else {
            // not a header => revert get pointer and break
            in.seekg(pos);
            break;
        }
    }

    if (!foundHeader) {
        std::cerr << "Error: no #CHROM line found.\n";
        return;
    }

    // Now read the variants
    std::string chrom, posStr, id, ref, alt, qual, filter, info, format;
    while (std::getline(in, line)) {
        if (line.empty() || line[0]=='#') {
            continue;
        }
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string t;
            while (std::getline(ss, t, '\t')) {
                fields.push_back(t);
            }
        }
        if (fields.size()<10) {
            std::cerr << "Warning: skipping line with <10 fields: " << line << "\n";
            continue;
        }

        chrom = fields[0];
        posStr= fields[1];
        id   = fields[2];
        ref  = fields[3];
        alt  = fields[4];
        // skip reading qual,filter,info,format, etc. for this simplistic approach
        // We'll parse sample columns from fields[9..]
        int posVal=0;
        try {
            posVal = std::stoi(posStr);
        } catch(...) {
            std::cerr << "Warning: invalid pos => skip " << line << "\n";
            continue;
        }

        // Build the genotype vector for this variant
        std::vector<int> genotypeVec;
        genotypeVec.reserve(fields.size()-9);
        for (size_t s=9; s<fields.size(); s++) {
            // e.g. "0/1"
            std::string gt = fields[s];
            // find slash or pipe
            size_t delim = gt.find_first_of("/|");
            if (delim==std::string::npos) {
                genotypeVec.push_back(-1);
                continue;
            }
            std::string a1 = gt.substr(0, delim);
            std::string a2 = gt.substr(delim+1);
            if (a1.empty() || a2.empty() || a1=="." || a2==".") {
                genotypeVec.push_back(-1);
                continue;
            }
            int i1=0, i2=0;
            try {
                i1= std::stoi(a1);
                i2= std::stoi(a2);
            } catch(...) {
                genotypeVec.push_back(-1);
                continue;
            }
            genotypeVec.push_back(i1 + i2); // naive approach
        }
        VariantData v;
        v.chrom= chrom;
        v.pos=   posVal;
        v.genotype = genotypeVec;
        variantList.push_back(v);
    }

    if (variantList.empty()) {
        std::cerr << "Error: no variant data found.\n";
        return;
    }

    // group variants
    auto blocks = groupVariants(variantList, ldThreshold);

    // we also output # the line "#HAPLOTYPE_BLOCKS" or something
    out << "#HAPLOTYPE_BLOCKS_START\n";
    for (size_t b=0; b<blocks.size(); b++) {
        out << "Block " << (b+1) << ": ";
        // We'll list the variants with index, chrom, pos
        // e.g. "0:(chr1:1000), 1:(chr1:1050)"
        for (size_t i=0; i<blocks[b].size(); i++) {
            int idx = blocks[b][i];
            out << idx << ":(" << variantList[idx].chrom << ":" << variantList[idx].pos << ")";
            if (i+1<blocks[b].size()) out << ", ";
        }
        out << "\n";
    }
    out << "#HAPLOTYPE_BLOCKS_END\n";
}

std::vector<std::vector<int>> VCFXHaplotypePhaser::groupVariants(const std::vector<VariantData>& variants,
                                                                 double ldThreshold)
{
    std::vector<std::vector<int>> blocks;
    std::vector<int> currentBlock;
    for (size_t i=0; i<variants.size(); i++) {
        if (currentBlock.empty()) {
            currentBlock.push_back(i);
        } else {
            // check LD with last variant in current block
            int lastIdx = currentBlock.back();
            double r2 = calculateLD(variants[lastIdx], variants[i]);
            if (r2 >= ldThreshold) {
                currentBlock.push_back(i);
            } else {
                // new block
                blocks.push_back(currentBlock);
                currentBlock.clear();
                currentBlock.push_back(i);
            }
        }
    }
    if (!currentBlock.empty()) {
        blocks.push_back(currentBlock);
    }
    return blocks;
}

double VCFXHaplotypePhaser::calculateLD(const VariantData& v1, const VariantData& v2) {
    // We compute r^2
    // ignoring missing (-1)
    int n=0;
    long sumX=0, sumY=0, sumXY=0, sumX2=0, sumY2=0;
    const auto &g1= v1.genotype;
    const auto &g2= v2.genotype;
    if (g1.size()!=g2.size()) {
        return 0.0;
    }
    for (size_t s=0; s<g1.size(); s++) {
        int x = g1[s];
        int y = g2[s];
        if (x<0 || y<0) {
            continue; // missing
        }
        n++;
        sumX+= x;
        sumY+= y;
        sumXY += x*y;
        sumX2 += x*x;
        sumY2 += y*y;
    }
    if (n==0) return 0.0;
    double meanX= (double)sumX/n;
    double meanY= (double)sumY/n;
    double cov  = ((double)sumXY/n) - (meanX*meanY);
    double varX = ((double)sumX2/n) - (meanX*meanX);
    double varY = ((double)sumY2/n) - (meanY*meanY);
    if (varX<=0.0 || varY<=0.0) {
        return 0.0;
    }
    double r = cov/(std::sqrt(varX)*std::sqrt(varY));
    return r*r;
}

int main(int argc, char* argv[]) {
    VCFXHaplotypePhaser hp;
    return hp.run(argc, argv);
}
