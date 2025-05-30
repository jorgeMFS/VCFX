#ifndef VCFX_ALIGNMENT_CHECKER_H
#define VCFX_ALIGNMENT_CHECKER_H

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// VCFXAlignmentChecker: Header file for Reference Alignment Discrepancy Finder Tool
class VCFXAlignmentChecker {
  public:
    // Entry point for the tool
    int run(int argc, char *argv[]);

  private:
    // Displays the help message
    void displayHelp();

    // Loads the reference genome from a FASTA file
    bool loadReferenceGenome(const std::string &path);

    // Checks discrepancies between VCF variants and the in-memory reference genome
    void checkDiscrepancies(std::istream &vcfIn, std::ostream &out);

    // Retrieves the reference base(s) from the reference genome at a specific position
    std::string getReferenceBases(const std::string &chrom, int pos, int length = 1);

    // Stores the reference genome sequences, keyed by normalized chromosome name
    struct FastaIndexEntry {
        std::streampos offset = 0;    // file offset to first base
        std::size_t length = 0;       // total bases in sequence
        std::size_t basesPerLine = 0; // number of bases per line in FASTA
        std::size_t bytesPerLine = 0; // bytes per line including newline
    };

    std::unordered_map<std::string, FastaIndexEntry> referenceIndex;
    std::ifstream referenceStream;
    std::string referencePath;

    // Helper function to convert chromosome names to a consistent format
    std::string normalizeChromosome(const std::string &chrom);
};

#endif // VCFX_ALIGNMENT_CHECKER_H
