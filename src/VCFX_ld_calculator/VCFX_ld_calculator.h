#ifndef VCFX_LD_CALCULATOR_H
#define VCFX_LD_CALCULATOR_H

#include <iostream>
#include <string>
#include <vector>

// Holds minimal data for a single variant
struct LDVariant {
    std::string chrom;
    int pos;
    std::vector<int> genotype; // 0 => homRef, 1 => het, 2 => homAlt, -1 => missing
};

class VCFXLDCalculator {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();
    // Parse the region string "chr1:10000-20000"
    // If none provided, regionChrom will be empty => use all
    bool parseRegion(const std::string &regionStr, std::string &regionChrom, int &regionStart, int &regionEnd);

    // The main logic: read VCF, store genotype codes for in-range variants, compute r^2, output
    void computeLD(std::istream &in, std::ostream &out, const std::string &regionChrom, int regionStart, int regionEnd);

    // parse a single genotype string => code
    int parseGenotype(const std::string &s);

    // compute r^2 for two variant’s genotype arrays
    double computeRsq(const std::vector<int> &g1, const std::vector<int> &g2);
};

#endif // VCFX_LD_CALCULATOR_H
