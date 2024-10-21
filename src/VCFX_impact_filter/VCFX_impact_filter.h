#ifndef VCFX_IMPACT_FILTER_H
#define VCFX_IMPACT_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXImpactFilter: Header file for Variant Impact Filter tool
class VCFXImpactFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input based on the specified impact level
    void filterByImpact(std::istream& in, std::ostream& out, const std::string& targetImpact);
};

#endif // VCFX_IMPACT_FILTER_H
