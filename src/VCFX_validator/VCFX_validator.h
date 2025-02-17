#ifndef VCFX_VALIDATOR_H
#define VCFX_VALIDATOR_H

#include <iostream>
#include <string>

class VCFXValidator {
public:
    int run(int argc, char* argv[]);

private:
    // If we add advanced checks for e.g. "strict" mode, we store a bool
    bool strictMode = false;

    // Show usage
    void displayHelp();

    // Main function that reads lines from in, does validation
    bool validateVCF(std::istream &in);

    // Validate a meta line "##" and returns true if it’s correct
    bool validateMetaLine(const std::string &line, int lineNumber);

    // Check #CHROM line
    bool validateChromHeader(const std::string &line, int lineNumber);

    // Validate a data line with at least 8 columns
    bool validateDataLine(const std::string &line, int lineNumber);
};

#endif
