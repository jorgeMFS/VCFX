#ifndef VCFX_FIELD_EXTRACTOR_H
#define VCFX_FIELD_EXTRACTOR_H

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>

// Displays usage/help
void printHelp();

// Main extraction function that reads from 'in', writes to 'out', extracting
// user-specified fields (including standard fields, INFO subfields, and sample subfields).
void extractFields(std::istream& in, std::ostream& out, const std::vector<std::string>& fields);

#endif // VCFX_FIELD_EXTRACTOR_H
