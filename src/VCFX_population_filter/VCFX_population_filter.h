#ifndef VCFX_POPULATION_FILTER_H
#define VCFX_POPULATION_FILTER_H

#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

class VCFXPopulationFilter {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();

    bool loadPopulationMap(const std::string &popMapFile, const std::string &popTag,
                           std::unordered_set<std::string> &samplesToInclude);
    void filterPopulation(std::istream &in, std::ostream &out, const std::unordered_set<std::string> &samplesToInclude,
                          const std::string &popTag);
};

#endif
