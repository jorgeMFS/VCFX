#ifndef VCFX_INFO_AGGREGATOR_H
#define VCFX_INFO_AGGREGATOR_H
#include <iostream>
#include <string>
#include <vector>
// VCFXInfoAggregator: Header file for INFO Field Aggregator tool
class VCFXInfoAggregator {
public:

    // Entry point for the tool
    int run(int argc, char* argv[]);

private:

    // Displays the help message
    void displayHelp();

    // Aggregates specified INFO fields across samples
    void aggregateInfo(std::istream& in, std::ostream& out, const std::vector<std::string>& infoFields);
};
#endif // VCFX_INFO_AGGREGATOR_H