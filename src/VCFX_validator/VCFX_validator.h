#ifndef VCFX_VALIDATOR_H
#define VCFX_VALIDATOR_H

#include <iostream>
#include <string>

// Function to validate VCF header
bool validateVCFHeader(const std::string& line);

// Function to validate a single VCF record
bool validateVCFRecord(const std::string& line, int line_number);

// Function to validate the entire VCF file
bool validateVCF(std::istream& in, std::ostream& out);

#endif // VCFX_VALIDATOR_H
