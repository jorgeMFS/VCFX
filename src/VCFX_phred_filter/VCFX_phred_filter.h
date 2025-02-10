#ifndef VCFX_PHRED_FILTER_H
#define VCFX_PHRED_FILTER_H

#include <iostream>
#include <string>

class VCFXPhredFilter {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();
    void processVCF(std::istream &in, double threshold, bool keepMissingAsPass);
    double parseQUAL(const std::string &qualStr, bool keepMissingAsPass);
};

#endif
