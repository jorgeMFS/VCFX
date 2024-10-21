#ifndef VCFX_INDEL_NORMALIZER_H
#define VCFX_INDEL_NORMALIZER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXIndelNormalizer: Header file for Indel Normalization Tool
class VCFXIndelNormalizer {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Normalizes indels in the VCF input
    void normalizeIndels(std::istream& in, std::ostream& out);

    // Determines if a variant is an indel
    bool isIndel(const std::string& ref, const std::string& alt);

    // Normalizes a single indel variant
    bool normalizeVariant(std::string& chrom, std::string& pos, std::string& ref, std::string& alt);
};

#endif // VCFX_INDEL_NORMALIZER_H
