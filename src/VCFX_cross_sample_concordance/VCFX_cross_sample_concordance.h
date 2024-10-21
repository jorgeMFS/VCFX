#ifndef VCFX_CROSS_SAMPLE_CONCORDANCE_H
#define VCFX_CROSS_SAMPLE_CONCORDANCE_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// VCFXCrossSampleConcordance: Header file for Cross-Sample Variant Concordance Tool
class VCFXCrossSampleConcordance {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Processes the VCF input and calculates concordance
    void calculateConcordance(std::istream& in, std::ostream& out);

    // Structure to hold variant information
    struct Variant {
        std::string chrom;
        std::string pos;
        std::string ref;
        std::string alt;
        std::vector<std::string> genotypes;
    };

    // Stores all variants
    std::vector<Variant> variants;
};

#endif // VCFX_CROSS_SAMPLE_CONCORDANCE_H