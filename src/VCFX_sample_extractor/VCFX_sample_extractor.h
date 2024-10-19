#ifndef VCFX_SAMPLE_EXTRACTOR_H
#define VCFX_SAMPLE_EXTRACTOR_H

#include <iostream>
#include <string>
#include <vector>

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], std::string& sample_name);

// Function to extract sample data from VCF
void extractSampleData(std::istream& in, std::ostream& out, const std::string& sample_name);

#endif // VCFX_SAMPLE_EXTRACTOR_H
