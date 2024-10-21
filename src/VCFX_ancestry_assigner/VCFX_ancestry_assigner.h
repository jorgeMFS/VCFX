#ifndef VCFX_ANCESTRY_ASSIGNER_H
#define VCFX_ANCESTRY_ASSIGNER_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// VCFXAncestryAssigner: Header file for Ancestry Assignment Tool
class VCFXAncestryAssigner {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Parses a single line from the frequency file
    bool parseFrequencyLine(const std::string& line, std::string& chrom, int& pos, char& ref, char& alt, std::unordered_map<std::string, double>& freqMap, const std::vector<std::string>& populations);

    // Loads ancestral frequencies from the provided input stream
    bool loadAncestralFrequencies(std::istream& in);

    // Assigns ancestry to samples based on VCF input
    void assignAncestry(std::istream& vcfIn, std::ostream& out);

    // Stores populations
    std::vector<std::string> populations;

    // Stores variant frequencies: Key -> Population -> Frequency
    std::unordered_map<std::string, std::unordered_map<std::string, double>> variantFrequencies;
};

#endif // VCFX_ANCESTRY_ASSIGNER_H
