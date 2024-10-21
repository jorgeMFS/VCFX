#ifndef VCFX_POPULATION_FILTER_H
#define VCFX_POPULATION_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXPopulationFilter: Header file for Population Subset Filter tool
class VCFXPopulationFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input to include only samples from the specified population
    void filterPopulation(std::istream& in, std::ostream& out, const std::string& populationTag, const std::string& popMapFile);
};

#endif // VCFX_POPULATION_FILTER_H
