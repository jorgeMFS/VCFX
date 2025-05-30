#ifndef VCFX_INFO_AGGREGATOR_H
#define VCFX_INFO_AGGREGATOR_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

// VCFXInfoAggregator: A tool that reads a VCF, prints it unmodified, and at the end
// produces an aggregated summary of numeric INFO fields (like DP, AF, etc.)
class VCFXInfoAggregator {
  public:
    int run(int argc, char *argv[]);

  private:
    // Print usage
    void displayHelp();

    // The main function that scans the VCF from 'in', writes lines unchanged to 'out',
    // collecting numeric values from specified fields in 'infoFields'.
    // After reading the entire file, it appends a summary section.
    void aggregateInfo(std::istream &in, std::ostream &out, const std::vector<std::string> &infoFields);
};

#endif // VCFX_INFO_AGGREGATOR_H
