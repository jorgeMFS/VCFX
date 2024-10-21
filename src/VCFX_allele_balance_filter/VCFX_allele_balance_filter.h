#ifndef VCFX_ALLELE_BALANCE_FILTER_H
#define VCFX_ALLELE_BALANCE_FILTER_H

#include <iostream>
#include <string>
#include <vector>
#include <sstream>

// VCFXAlleleBalanceFilter: Header file for Allele Balance Filter Tool
class VCFXAlleleBalanceFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input based on allele balance threshold
    void filterByAlleleBalance(std::istream& in, std::ostream& out, double threshold);

    // Parses the genotype to calculate allele balance
    double calculateAlleleBalance(const std::string& genotype);
};

#endif // VCFX_ALLELE_BALANCE_FILTER_H
