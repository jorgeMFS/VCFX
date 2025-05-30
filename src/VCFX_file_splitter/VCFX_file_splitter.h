#ifndef VCFX_FILE_SPLITTER_H
#define VCFX_FILE_SPLITTER_H

#include <iostream>
#include <string>

// VCFXFileSplitter: Splits a VCF file by chromosome into multiple smaller VCFs.
class VCFXFileSplitter {
  public:
    // Entry point for the tool
    int run(int argc, char *argv[]);

  private:
    // Displays the help message
    void displayHelp();

    // Splits the input VCF by chromosome using the specified prefix
    void splitVCFByChromosome(std::istream &in, const std::string &outputPrefix);
};

#endif // VCFX_FILE_SPLITTER_H
