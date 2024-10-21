#ifndef VCFX_REF_COMPARATOR_H
#define VCFX_REF_COMPARATOR_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// VCFXRefComparator: Header file for Reference Genome Comparator tool
class VCFXRefComparator {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Loads the reference genome from a FASTA file
    bool loadReference(const std::string& referencePath);

    // Compares VCF variants against the reference genome
    void compareWithReference(std::istream& vcfInput, std::ostream& vcfOutput);

    // Reference genome data: chromosome -> sequence
    std::unordered_map<std::string, std::string> referenceGenome;
};

#endif // VCFX_REF_COMPARATOR_H
