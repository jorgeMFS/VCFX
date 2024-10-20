#ifndef VCFX_CONCORDANCE_CHECKER_H
#define VCFX_CONCORDANCE_CHECKER_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Structure to hold command-line arguments
struct ConcordanceArguments {
    std::string sample1;
    std::string sample2;
};

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], ConcordanceArguments& args);

// Function to split a string by a delimiter
std::vector<std::string> splitString(const std::string& str, char delimiter);

// Function to process VCF and calculate concordance
bool calculateConcordance(std::istream& in, std::ostream& out, const ConcordanceArguments& args);

// Function to extract genotype from genotype string
std::string extractGenotype(const std::string& genotype_str);

#endif // VCFX_CONCORDANCE_CHECKER_H
