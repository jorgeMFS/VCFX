#ifndef VCFX_PROBABILITY_FILTER_H
#define VCFX_PROBABILITY_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXProbabilityFilter: Header file for Genotype Probability Filter tool
class VCFXProbabilityFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input based on the specified genotype probability condition
    void filterByProbability(std::istream& in, std::ostream& out, const std::string& condition);
};

#endif // VCFX_PROBABILITY_FILTER_H
