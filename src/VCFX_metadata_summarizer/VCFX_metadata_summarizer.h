#ifndef VCFX_METADATA_SUMMARIZER_H
#define VCFX_METADATA_SUMMARIZER_H

#include <iostream>
#include <string>
#include <vector>
#include <map>

// VCFX_metadata_summarizer: Header file for VCF metadata summarization tool
class VCFXMetadataSummarizer {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Processes VCF input and summarizes metadata
    void summarizeMetadata(std::istream& in);

    // Parses the header lines to extract metadata
    void parseHeader(const std::string& line, std::map<std::string, int>& metadata);

    // Prints the metadata summary
    void printSummary(const std::map<std::string, int>& metadata);
};

#endif // VCFX_METADATA_SUMMARIZER_H
