#ifndef VCFX_FORMAT_CONVERTER_H
#define VCFX_FORMAT_CONVERTER_H

#include <iostream>
#include <string>
#include <vector>

// Enumeration for output formats
enum class OutputFormat {
    BED,
    CSV
};

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], OutputFormat& format);

// Function to convert VCF to BED
void convertVCFtoBED(std::istream& in, std::ostream& out);

// Function to convert VCF to CSV
void convertVCFtoCSV(std::istream& in, std::ostream& out);

#endif // VCFX_FORMAT_CONVERTER_H
