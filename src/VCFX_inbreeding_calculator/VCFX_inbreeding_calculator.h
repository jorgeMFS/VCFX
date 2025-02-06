#ifndef VCFX_INBREEDING_CALCULATOR_H
#define VCFX_INBREEDING_CALCULATOR_H

#include <iostream>
#include <string>
#include <vector>

// Represents a single VCF variant for biallelic analysis
struct InbreedingVariant {
    std::string chrom;
    int pos;
    // genotypeCodes[sampleIndex] in { -1,0,1,2 } => missing,0/0,0/1,1/1
    std::vector<int> genotypeCodes;
};

// VCFXInbreedingCalculator: calculates individual inbreeding coefficients
class VCFXInbreedingCalculator {
public:
    int run(int argc, char* argv[]);

private:
    // Print help
    void displayHelp();

    // Main function: read VCF => store biallelic variants => compute F
    void calculateInbreeding(std::istream& in, std::ostream& out);

    // Utility to parse a single genotype string => 0,1,2, or -1
    int parseGenotype(const std::string& s);

    // Helper to decide if ALT is biallelic
    bool isBiallelic(const std::string &alt);

    // Summation step: for each sample, we sum observed het and sum of expected
    //   expected is sum(2 p_excl(1 - p_excl)) across variants
};

#endif // VCFX_INBREEDING_CALCULATOR_H
