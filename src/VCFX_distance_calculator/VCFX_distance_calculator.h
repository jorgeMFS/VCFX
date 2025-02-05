// VCFX_distance_calculator.h
#ifndef VCFX_DISTANCE_CALCULATOR_H
#define VCFX_DISTANCE_CALCULATOR_H

#include <iostream>
#include <string>

// Function to display help message
void printHelp();

// Structure to hold variant information (only CHROM and POS are used)
struct VCFVariant {
    std::string chrom;
    int pos;
};

// Function to parse a VCF line into a VCFVariant structure
bool parseVCFLine(const std::string& line, VCFVariant& variant);

// Function to calculate distances between consecutive variants
// Writes one output line per variant (with distance from the previous variant on that chromosome)
// and returns true on success.
bool calculateDistances(std::istream& in, std::ostream& out);

#endif // VCFX_DISTANCE_CALCULATOR_H
