#ifndef VCFX_POSITION_SUBSETTER_H
#define VCFX_POSITION_SUBSETTER_H

#include <iostream>
#include <string>

// Function to display help message
void printHelp();

// Function to subset VCF records based on genomic range
bool subsetVCFByPosition(std::istream& in, std::ostream& out, const std::string& region);

#endif // VCFX_POSITION_SUBSETTER_H
