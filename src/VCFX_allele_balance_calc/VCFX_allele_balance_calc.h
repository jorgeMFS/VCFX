#ifndef VCFX_ALLELE_BALANCE_CALC_H
#define VCFX_ALLELE_BALANCE_CALC_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Structure to hold command-line arguments
struct AlleleBalanceArguments {
    std::vector<std::string> samples;
};

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], AlleleBalanceArguments& args);

// Function to split a string by a delimiter
std::vector<std::string> splitString(const std::string& str, char delimiter);

// Function to process VCF and calculate allele balance
bool calculateAlleleBalance(std::istream& in, std::ostream& out, const AlleleBalanceArguments& args);

// Function to extract genotype from genotype string
std::string extractGenotype(const std::string& genotype_str);

// Function to calculate allele balance ratio
double computeAlleleBalance(const std::string& genotype);

#endif // VCFX_ALLELE_BALANCE_CALC_H
