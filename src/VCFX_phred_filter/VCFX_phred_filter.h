#ifndef VCFX_PHRED_FILTER_H
#define VCFX_PHRED_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFX_phred_filter: Header file for Phred score filtering tool
class VCFXPhredFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Processes VCF input from the given stream
    void processVCF(std::istream& in, double threshold);

    // Parses the QUAL field
    double parseQUAL(const std::string& qualStr);
};

#endif // VCFX_PHRED_FILTER_H
