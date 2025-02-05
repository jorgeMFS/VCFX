#ifndef VCFX_FASTA_CONVERTER_H
#define VCFX_FASTA_CONVERTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXFastaConverter: Tool for converting a "variant-only" VCF into per-sample FASTA sequences
class VCFXFastaConverter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Converts VCF input to FASTA format
    void convertVCFtoFasta(std::istream& in, std::ostream& out);
};

#endif // VCFX_FASTA_CONVERTER_H
