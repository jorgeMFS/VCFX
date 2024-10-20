#ifndef VCFX_ALLELE_FREQ_CALC_H
#define VCFX_ALLELE_FREQ_CALC_H

#include <iostream>
#include <string>

// Function to display help message
void printHelp();

// Function to perform allele frequency calculation on VCF records
void calculateAlleleFrequency(std::istream& in, std::ostream& out);

#endif // VCFX_ALLELE_FREQ_CALC_H
