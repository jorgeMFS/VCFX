#ifndef VCFX_HAPLOTYPE_EXTRACTOR_H
#define VCFX_HAPLOTYPE_EXTRACTOR_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Structure to represent a haplotype block
struct HaplotypeBlock {
    std::string chrom;
    int start;
    int end;
    std::vector<std::string> haplotypes;
};

// Function to extract haplotype blocks from VCF
bool extractHaplotypeBlocks(std::istream& in, std::ostream& out);

#endif // VCFX_HAPLOTYPE_EXTRACTOR_H
