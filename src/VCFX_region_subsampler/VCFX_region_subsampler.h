#ifndef VCFX_REGION_SUBSAMPLER_H
#define VCFX_REGION_SUBSAMPLER_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// VCFXRegionSubsampler: Header file for Region-based Subsampling Tool
class VCFXRegionSubsampler {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Loads regions from a BED file into a map
    bool loadRegions(const std::string& bedFilePath, std::unordered_map<std::string, std::vector<std::pair<int, int>>>& regions);

    // Checks if a variant falls within any of the specified regions
    bool isVariantInRegions(const std::string& chrom, int pos, const std::unordered_map<std::string, std::vector<std::pair<int, int>>>& regions);

    // Subsamples variants from the VCF input based on regions
    void subsampleRegions(std::istream& in, std::ostream& out, const std::unordered_map<std::string, std::vector<std::pair<int, int>>>& regions);
};

#endif // VCFX_REGION_SUBSAMPLER_H
