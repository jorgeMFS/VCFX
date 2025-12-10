#ifndef VCFX_REGION_SUBSAMPLER_H
#define VCFX_REGION_SUBSAMPLER_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// Each region is (start, end), inclusive, 1-based
struct Region {
    int start;
    int end;
};

// VCFXRegionSubsampler:
// Reads a BED file with multiple lines => chromosome -> sorted intervals
// Then reads a VCF and keeps lines whose POS is within any interval for that CHROM.
class VCFXRegionSubsampler {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();
    bool loadRegions(const std::string &bedFile, std::unordered_map<std::string, std::vector<Region>> &chromRegions);
    // after loading, we sort the intervals for each chromosome, possibly merge them
    void sortAndMergeIntervals(std::unordered_map<std::string, std::vector<Region>> &chromRegions);

    bool isInAnyRegion(const std::string &chrom, int pos) const;

    // The map from chrom => sorted intervals
    std::unordered_map<std::string, std::vector<Region>> regions;

    void processVCF(std::istream &in, std::ostream &out);

    // Memory-mapped file processing (fast path)
    bool processVCFMmap(const char* filepath, std::ostream& out, bool quiet = false);
};

#endif
