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

// Frequency mode for computing inbreeding
enum class FrequencyMode {
    EXCLUDE_SAMPLE, // exclude the sample's own genotype from p
    GLOBAL          // use a single site-wide p for all samples
};

// VCFXInbreedingCalculator: calculates individual inbreeding coefficients
class VCFXInbreedingCalculator {
public:
    int run(int argc, char* argv[]);

private:
    // Command-line settings
    FrequencyMode freqMode_          = FrequencyMode::EXCLUDE_SAMPLE;
    bool skipBoundary_               = false;
    bool countBoundaryAsUsed_        = false;

    // Print help
    void displayHelp();

    // Main function: read VCF => store biallelic variants => compute F
    void calculateInbreeding(std::istream& in, std::ostream& out);

    // Utility to parse a single genotype string => 0,1,2, or -1
    int parseGenotype(const std::string& s);

    // Helper to decide if ALT is biallelic
    bool isBiallelic(const std::string &alt);

    // Command-line parser for freq-mode
    FrequencyMode parseFreqMode(const std::string &modeStr);
};

#endif // VCFX_INBREEDING_CALCULATOR_H
