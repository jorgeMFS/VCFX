#ifndef VCFX_NONREF_FILTER_H
#define VCFX_NONREF_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXNonRefFilter: Header file for Non-Reference Variant Filter tool
class VCFXNonRefFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input to exclude variants with all homozygous reference genotypes
    void filterNonRef(std::istream& in, std::ostream& out);
};

#endif // VCFX_NONREF_FILTER_H
