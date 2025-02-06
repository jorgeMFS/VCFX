#ifndef VCFX_GENOTYPE_QUERY_H
#define VCFX_GENOTYPE_QUERY_H

#include <iostream>
#include <string>

// Parses command-line arguments to get the desired genotype query and a 'strict' flag
bool parseArguments(int argc, char* argv[], std::string& genotype_query, bool &strictCompare);

// Prints usage/help
void printHelp();

// Filters a VCF from 'in', writing only records that contain at least one sample with the requested genotype
void genotypeQuery(std::istream& in, std::ostream& out, const std::string& genotype_query, bool strictCompare);

#endif // VCFX_GENOTYPE_QUERY_H
