#ifndef VCFX_MULTIALLELIC_SPLITTER_H
#define VCFX_MULTIALLELIC_SPLITTER_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Structure to hold VCF variant information
struct VCFVariant {
    std::string chrom;
    int pos;
    std::string id;
    std::string ref;
    std::vector<std::string> alt;
    std::string qual;
    std::string filter;
    std::string info;
    std::vector<std::string> other_fields; // FORMAT and sample data
};

// Function to parse a VCF line into VCFVariant structure
bool parseVCFLine(const std::string& line, VCFVariant& variant);

// Function to split multi-allelic variants into bi-allelic
bool splitMultiAllelicVariants(std::istream& in, std::ostream& out);

#endif // VCFX_MULTIALLELIC_SPLITTER_H
