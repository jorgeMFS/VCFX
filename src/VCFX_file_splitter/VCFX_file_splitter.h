#ifndef VCFX_FILE_SPLITTER_H
#define VCFX_FILE_SPLITTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXFileSplitter: Header file for VCF File Splitter tool
class VCFXFileSplitter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Splits VCF input based on chromosome and outputs to separate files
    void splitVCFByChromosome(std::istream& in, const std::string& outputPrefix);
};

#endif // VCFX_FILE_SPLITTER_H
