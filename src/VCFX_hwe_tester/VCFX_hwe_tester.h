#ifndef VCFX_HWE_TESTER_H
#define VCFX_HWE_TESTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFX_hwe_tester: Tool for Hardy-Weinberg Equilibrium (HWE) tests on a VCF.
class VCFXHWETester {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Reads a VCF from 'in' and writes results to stdout
    void performHWE(std::istream& in);

    // Checks if the line is strictly biallelic. If multi-allelic, skip.
    bool isBiallelic(const std::string &alt);

    // Parses genotype information to count 0/0, 0/1, 1/1 calls
    bool parseGenotypes(const std::vector<std::string>& genotypes, int& homRef, int& het, int& homAlt);

    // Calculate HWE p-value using an exact test
    double calculateHWE(int homRef, int het, int homAlt);

    // The “core” of the exact test probability
    double genotypeProbability(int homRef, int het, int homAlt);
};

#endif // VCFX_HWE_TESTER_H
