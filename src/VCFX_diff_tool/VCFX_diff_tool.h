#ifndef VCFX_DIFF_TOOL_H
#define VCFX_DIFF_TOOL_H

#include <iostream>
#include <string>
#include <fstream>
#include <unordered_set>
#include <vector>

// VCFXDiffTool: Header file for VCF Diff Tool
class VCFXDiffTool {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Loads variants from a VCF file into a set
    bool loadVariants(const std::string& filePath, std::unordered_set<std::string>& variants);

    // Generates a unique key for a variant based on chromosome, position, ref, and alt
    std::string generateVariantKey(const std::string& chrom, const std::string& pos, const std::string& ref, const std::string& alt);
};

#endif // VCFX_DIFF_TOOL_H
