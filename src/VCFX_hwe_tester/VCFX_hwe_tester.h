#ifndef VCFX_HWE_TESTER_H
#define VCFX_HWE_TESTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFX_hwe_tester: Header file for Hardy-Weinberg Equilibrium (HWE) testing tool
class VCFXHWETester {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Processes VCF input and performs HWE tests
    void performHWE(std::istream& in);

    // Parses genotype information to count genotypes
    bool parseGenotypes(const std::vector<std::string>& genotypes, int& homRef, int& het, int& homAlt);

    // Calculates HWE p-value using the exact test
    double calculateHWE(int homRef, int het, int homAlt);

    // Computes the probability of a given genotype count
    double genotypeProbability(int homRef, int het, int homAlt);

    // Computes log factorial to aid in probability calculations
    double logFactorial(int n);

    // Log-sum-exp trick for numerical stability
    double logSumExp(const std::vector<double>& logProbs);
};

#endif // VCFX_HWE_TESTER_H
