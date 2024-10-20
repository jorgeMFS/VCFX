#ifndef VCFX_DISTANCE_CALCULATOR_H
#define VCFX_DISTANCE_CALCULATOR_H

#include <iostream>
#include <string>

// Function to display help message
void printHelp();

// Structure to hold variant information
struct VCFVariant {
    std::string chrom;
    int pos;
};

// Function to parse a VCF line into VCFVariant structure
bool parseVCFLine(const std::string& line, VCFVariant& variant);

// Function to calculate distances between consecutive variants
bool calculateDistances(std::istream& in, std::ostream& out);

#endif // VCFX_DISTANCE_CALCULATOR_H
