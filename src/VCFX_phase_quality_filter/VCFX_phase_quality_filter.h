#ifndef VCFX_PHASE_QUALITY_FILTER_H
#define VCFX_PHASE_QUALITY_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXPhaseQualityFilter: Header file for Variant Phasing Quality Filter Tool
class VCFXPhaseQualityFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input based on phasing quality scores
    void filterByPQ(std::istream& in, std::ostream& out, double threshold);

    // Parses the PQ score from the INFO field
    double parsePQScore(const std::string& infoField);
};

#endif // VCFX_PHASE_QUALITY_FILTER_H
