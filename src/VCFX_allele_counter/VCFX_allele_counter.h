#ifndef VCFX_ALLELE_COUNTER_H
#define VCFX_ALLELE_COUNTER_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Structure to hold command-line arguments
struct AlleleCounterArguments {
    std::vector<std::string> samples;
};

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], AlleleCounterArguments& args);

// Function to split a string by a delimiter
std::vector<std::string> splitString(const std::string& str, char delimiter);

// Function to process VCF and count alleles
bool countAlleles(std::istream& in, std::ostream& out, const AlleleCounterArguments& args);

// Function to extract genotype from genotype string
std::string extractGenotype(const std::string& genotype_str);

#endif // VCFX_ALLELE_COUNTER_H
