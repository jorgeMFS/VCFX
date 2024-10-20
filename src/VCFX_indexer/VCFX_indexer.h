#ifndef VCFX_INDEXER_H
#define VCFX_INDEXER_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// Structure to hold variant information for indexing
struct VariantIndex {
    std::string chrom;
    int pos;
    long file_offset; // Byte offset in the file
};

// Function to display help message
void printHelp();

// Function to create an index for VCF file
void createVCFIndex(std::istream& in, std::ostream& out);

#endif // VCFX_INDEXER_H
