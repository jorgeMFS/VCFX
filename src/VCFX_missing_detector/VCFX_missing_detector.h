#ifndef VCFX_MISSING_DETECTOR_H
#define VCFX_MISSING_DETECTOR_H

#include <iostream>
#include <string>
#include <vector>

/**
 * VCFXMissingDetector: A tool to detect variants with any missing sample genotypes
 * and flag them in the INFO field with MISSING_GENOTYPES=1
 */
class VCFXMissingDetector {
public:
    /**
     * Entry point for the tool
     */
    int run(int argc, char* argv[]);

private:
    /**
     * Displays the help message
     */
    void displayHelp();

    /**
     * Detects missing genotypes in VCF input from 'in', writes flagged lines to 'out'
     */
    void detectMissingGenotypes(std::istream& in, std::ostream& out);
};

#endif // VCFX_MISSING_DETECTOR_H
