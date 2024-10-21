#ifndef VCFX_MERGER_H
#define VCFX_MERGER_H

#include <iostream>
#include <string>
#include <vector>

// VCFX_merger: Header file for VCF file merging tool
class VCFXMerger {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Processes and merges VCF files
    void mergeVCF(const std::vector<std::string>& inputFiles, std::ostream& out);

    // Parses a VCF file and stores variants
    void parseVCF(const std::string& filename, std::vector<std::vector<std::string>>& variants, std::vector<std::string>& headerLines);

    // Compares variants based on chromosome and position
    bool compareVariants(const std::vector<std::string>& a, const std::vector<std::string>& b);
};

#endif // VCFX_MERGER_H
