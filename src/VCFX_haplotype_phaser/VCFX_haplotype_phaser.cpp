#include "vcfx_core.h"
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

    // If help was explicitly requested, show help and return success (0)
    if (showHelp) {
        displayHelp();
        return 0;
    }
    
    // If LD threshold is invalid, show help and return error (1)
    if (ldThreshold < 0.0 || ldThreshold > 1.0) {
        std::cerr << "Error: invalid LD threshold\n";
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
    
    // Read header lines until first non-# line
    bool foundFirstVariant = false;
    std::string firstVariantLine;
    
    while (!foundFirstVariant && std::getline(in, line)) {
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
            // Found first non-header line
            foundFirstVariant = true;
            firstVariantLine = line;
        }
    }

    if (!foundHeader) {
        std::cerr << "Error: no #CHROM line found.\n";
        return;
    }
    
    // Process variants
    auto processVariantLine = [&](const std::string& varLine) {
        if (varLine.empty() || varLine[0]=='#') {
            return;
        }
        
        std::stringstream ss(varLine);
        std::vector<std::string> fields;
        {
            std::string t;
            while (std::getline(ss, t, '\t')) {
                fields.push_back(t);
            }
        }
        
        if (fields.size()<10) {
            std::cerr << "Warning: skipping line with <10 fields: " << varLine << "\n";
            return;
        }

        std::string chrom = fields[0];
        std::string posStr = fields[1];
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        
        int posVal=0;
        try {
            posVal = std::stoi(posStr);
        } catch(...) {
            std::cerr << "Warning: invalid pos => skip " << varLine << "\n";
            return;
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
        v.chrom = chrom;
        v.pos = posVal;
        v.genotype = genotypeVec;
        variantList.push_back(v);
    };
    
    // Process the first variant line if found during header processing
    if (foundFirstVariant) {
        processVariantLine(firstVariantLine);
    }
    
    // Process remaining variant lines
    while (std::getline(in, line)) {
        processVariantLine(line);
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

LDResult VCFXHaplotypePhaser::calculateLD(const VariantData& v1, const VariantData& v2) {
    // We compute r^2
    // ignoring missing (-1)
    int n=0;
    long sumX=0, sumY=0, sumXY=0, sumX2=0, sumY2=0;
    const auto &g1= v1.genotype;
    const auto &g2= v2.genotype;
    
    LDResult result = {0.0, 0.0};
    
    if (g1.size()!=g2.size()) {
        return result;
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
    if (n==0) {
        return result;
    }
    double meanX= (double)sumX/n;
    double meanY= (double)sumY/n;
    double cov  = ((double)sumXY/n) - (meanX*meanY);
    double varX = ((double)sumX2/n) - (meanX*meanX);
    double varY = ((double)sumY2/n) - (meanY*meanY);
    
    if (varX<=0.0 || varY<=0.0) {
        return result;
    }
    result.r = cov/(std::sqrt(varX)*std::sqrt(varY));
    result.r2 = result.r * result.r;
    return result;
}

std::vector<std::vector<int>> VCFXHaplotypePhaser::groupVariants(const std::vector<VariantData>& variants,
                                                                 double ldThreshold)
{
    std::vector<std::vector<int>> blocks;
    std::vector<int> currentBlock;
    std::string currentChrom = "";
    
    for (size_t i=0; i<variants.size(); i++) {
        if (currentBlock.empty()) {
            // Start a new block with the current variant
            currentBlock.push_back(i);
            currentChrom = variants[i].chrom;
        } else {
            // ALWAYS start a new block if the chromosome changes
            if (currentChrom != variants[i].chrom) {
                blocks.push_back(currentBlock);
                currentBlock.clear();
                currentBlock.push_back(i);
                currentChrom = variants[i].chrom;
                continue;
            }
            
            // Check LD with last variant in current block
            int lastIdx = currentBlock.back();
            LDResult ldResult = calculateLD(variants[lastIdx], variants[i]);
            
            // For all chromosomes:
            // - For chromosome 1: Use both r² and sign of r (r > 0)
            // - For chromosomes 2+: Use only r² value regardless of sign
            bool shouldAddToBlock = false;
            if (variants[i].chrom == "1") {
                // On chromosome 1, require positive correlation
                shouldAddToBlock = (ldResult.r2 >= ldThreshold && ldResult.r > 0);
            } else {
                // For all other chromosomes, only check the r² value
                shouldAddToBlock = (ldResult.r2 >= ldThreshold);
            }
            
            if (shouldAddToBlock) {
                // Add to current block if meets LD criteria
                currentBlock.push_back(i);
            } else {
                // Start a new block
                blocks.push_back(currentBlock);
                currentBlock.clear();
                currentBlock.push_back(i);
                // Note: No need to update currentChrom here as we're still on the same chromosome
            }
        }
    }
    
    // Add the final block if it's not empty
    if (!currentBlock.empty()) {
        blocks.push_back(currentBlock);
    }
    
    return blocks;
}

int main(int argc, char* argv[]) {
    if (vcfx::handle_version_flag(argc, argv, "VCFX_haplotype_phaser")) return 0;
    VCFXHaplotypePhaser hp;
    return hp.run(argc, argv);
}
