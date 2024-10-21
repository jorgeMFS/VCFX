#ifndef VCFX_ALIGNMENT_CHECKER_H
#define VCFX_ALIGNMENT_CHECKER_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// VCFXAlignmentChecker: Header file for Reference Alignment Discrepancy Finder Tool
class VCFXAlignmentChecker {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Loads the reference genome from a FASTA file
    bool loadReferenceGenome(std::istream& in);

    // Checks discrepancies between VCF variants and the reference genome
    void checkDiscrepancies(std::istream& vcfIn, std::istream& refIn, std::ostream& out);

    // Retrieves the reference base(s) from the reference genome at a specific position
    std::string getReferenceBases(const std::string& chrom, int pos, int length = 1);

    // Stores the reference genome sequences
    std::unordered_map<std::string, std::string> referenceGenome;

    // Helper function to convert chromosome names to consistent format
    std::string normalizeChromosome(const std::string& chrom);
};

#endif // VCFX_ALIGNMENT_CHECKER_H
