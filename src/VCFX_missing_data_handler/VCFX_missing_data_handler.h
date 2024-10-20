#ifndef VCFX_MISSING_DATA_HANDLER_H
#define VCFX_MISSING_DATA_HANDLER_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Structure to hold command-line arguments
struct Arguments {
    bool fill_missing = false;
    std::string default_genotype = "./.";
    std::vector<std::string> input_files;
};

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], Arguments& args);

// Function to process VCF and handle missing data
bool handleMissingData(std::istream& in, std::ostream& out, const Arguments& args);

// Utility function to split a string by a delimiter
std::vector<std::string> splitString(const std::string& str, char delimiter);

#endif // VCFX_MISSING_DATA_HANDLER_H
