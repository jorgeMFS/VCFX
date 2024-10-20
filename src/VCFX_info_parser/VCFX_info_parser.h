#ifndef VCFX_INFO_PARSER_H
#define VCFX_INFO_PARSER_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], std::vector<std::string>& info_fields);

// Function to parse the INFO field and display selected fields
bool parseInfoFields(std::istream& in, std::ostream& out, const std::vector<std::string>& info_fields);

#endif // VCFX_INFO_PARSER_H
