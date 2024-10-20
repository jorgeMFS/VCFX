#ifndef VCFX_GENOTYPE_QUERY_H
#define VCFX_GENOTYPE_QUERY_H

#include <iostream>
#include <string>
#include <vector>

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], std::string& genotype_query);

// Function to display help message
void printHelp();

// Function to perform genotype query on VCF records
void genotypeQuery(std::istream& in, std::ostream& out, const std::string& genotype_query);

#endif // VCFX_GENOTYPE_QUERY_H
