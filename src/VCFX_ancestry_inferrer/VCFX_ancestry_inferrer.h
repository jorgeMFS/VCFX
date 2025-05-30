#ifndef VCFX_ANCESTRY_INFERRER_H
#define VCFX_ANCESTRY_INFERRER_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// VCFXAncestryInferer: Header file for Ancestry Inference tool
class VCFXAncestryInferer {
  public:
    // Entry point for the tool
    int run(int argc, char *argv[]);

  private:
    // Displays the help message
    void displayHelp();

    // Loads population allele frequencies from a file
    bool loadPopulationFrequencies(const std::string &freqFilePath);

    // Infers ancestry for each sample based on allele frequencies
    void inferAncestry(std::istream &vcfInput, std::ostream &ancestryOutput);

    // Structure to hold allele frequencies per population
    struct PopulationFrequencies {
        std::string population;
        std::unordered_map<std::string, double> variantFrequencies; // Key: "chr:pos:ref:alt", Value: frequency
    };

    std::vector<PopulationFrequencies> populations;
};

#endif // VCFX_ANCESTRY_INFERRER_H
