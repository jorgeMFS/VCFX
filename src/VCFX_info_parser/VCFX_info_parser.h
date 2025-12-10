#ifndef VCFX_INFO_PARSER_H
#define VCFX_INFO_PARSER_H

#include <iostream>
#include <string>
#include <vector>

// Prints the help message
void printHelp();

// Parses command-line arguments into 'info_fields'.
// Returns true if at least one info field was parsed successfully.
bool parseArguments(int argc, char *argv[], std::vector<std::string> &info_fields);

// Splits a string by a delimiter
std::vector<std::string> split(const std::string &s, char delimiter);

// Parses the VCF, extracting the selected INFO fields and printing to 'out'
bool parseInfoFields(std::istream &in, std::ostream &out, const std::vector<std::string> &info_fields);

// Memory-mapped file processing (fast path)
bool parseInfoFieldsMmap(const char* filepath, std::ostream& out,
                         const std::vector<std::string>& info_fields, bool quiet = false);

#endif // VCFX_INFO_PARSER_H
