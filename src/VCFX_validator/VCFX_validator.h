#ifndef VCFX_VALIDATOR_H
#define VCFX_VALIDATOR_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>

class VCFXValidator {
public:
    int run(int argc, char* argv[]);

private:
    // If we add advanced checks for e.g. "strict" mode, we store a bool
    bool strictMode = false;
    bool reportDuplicates = false;
    // Number of columns in the #CHROM header line
    int headerColumnCount = 0;
    // Whether the header includes FORMAT/sample columns
    bool headerHasFormat = false;
    // Number of sample columns
    int sampleCount = 0;

    struct FieldDef {
        std::string number;
        std::string type;
    };
    std::unordered_map<std::string, FieldDef> infoDefs;
    std::unordered_map<std::string, FieldDef> formatDefs;
    std::unordered_set<std::string> seenVariants;

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
