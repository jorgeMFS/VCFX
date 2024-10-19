#ifndef VCFX_FIELD_EXTRACTOR_H
#define VCFX_FIELD_EXTRACTOR_H

#include <iostream>
#include <string>
#include <vector>

// Function to parse and extract specified fields from VCF records
std::vector<std::string> parseFields(const std::string& record, const std::vector<std::string>& fields);

// Function to process and extract fields from VCF
void extractFields(std::istream& in, std::ostream& out, const std::vector<std::string>& fields);

#endif // VCFX_FIELD_EXTRACTOR_H
