#ifndef VCFX_INBREEDING_CALCULATOR_H
#define VCFX_INBREEDING_CALCULATOR_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// VCFXInbreedingCalculator: Header file for Inbreeding Coefficient Calculator Tool
class VCFXInbreedingCalculator {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Calculates inbreeding coefficients (F-statistics) from VCF input
    void calculateInbreedingCoefficients(std::istream& in, std::ostream& out);

    // Parses genotype and returns allele counts
    bool parseGenotype(const std::string& genotype, int& a1, int& a2);

    // Calculates Hardy-Weinberg Expected Heterozygosity
    double calculateExpectedHet(int totalAlleles, int homRef, int homAlt, int het);

    // Calculates F-statistic for an individual
    double calculateF(int homozygous, int heterozygous, double expectedHet);
};

#endif // VCFX_INBREEDING_CALCULATOR_H
