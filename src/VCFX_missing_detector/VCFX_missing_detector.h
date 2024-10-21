#ifndef VCFX_MISSING_DETECTOR_H
#define VCFX_MISSING_DETECTOR_H

#include <iostream>
#include <string>
#include <vector>

// VCFXMissingDetector: Header file for Missing Sample Detection Tool
class VCFXMissingDetector {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Detects missing genotypes in VCF input
    void detectMissingGenotypes(std::istream& in, std::ostream& out);
};

#endif // VCFX_MISSING_DETECTOR_H
