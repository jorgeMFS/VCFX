#ifndef VCFX_FIELD_EXTRACTOR_H
#define VCFX_FIELD_EXTRACTOR_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// Displays usage/help
void printHelp();

// Main extraction function that reads from 'in', writes to 'out', extracting
// user-specified fields (including standard fields, INFO subfields, and sample subfields).
void extractFields(std::istream &in, std::ostream &out, const std::vector<std::string> &fields);

#endif // VCFX_FIELD_EXTRACTOR_H
