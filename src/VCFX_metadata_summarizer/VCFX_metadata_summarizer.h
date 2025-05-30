#ifndef VCFX_METADATA_SUMMARIZER_H
#define VCFX_METADATA_SUMMARIZER_H

#include <iostream>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

// VCFX_metadata_summarizer: Summarizes key metadata from a VCF file.
class VCFXMetadataSummarizer {
  public:
    // Entry point for the tool
    int run(int argc, char *argv[]);

  private:
    // Displays the help message
    void displayHelp();

    // Processes VCF input and summarizes metadata
    void summarizeMetadata(std::istream &in);

    // Parses meta-information lines (##...) to extract contig/INFO/FILTER/FORMAT IDs
    void parseHeader(const std::string &line);

    // Prints the metadata summary
    void printSummary() const;

    // Data members that store the final results
    // We'll track unique contig IDs, info IDs, filter IDs, format IDs
    std::unordered_set<std::string> contigIDs;
    std::unordered_set<std::string> infoIDs;
    std::unordered_set<std::string> filterIDs;
    std::unordered_set<std::string> formatIDs;

    int numSamples = 0;
    int numVariants = 0;
};

#endif // VCFX_METADATA_SUMMARIZER_H
