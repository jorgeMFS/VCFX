#ifndef VCFX_HAPLOTYPE_PHASER_H
#define VCFX_HAPLOTYPE_PHASER_H

#include <iostream>
#include <string>
#include <vector>

// A small struct to store a variant's key data: chromosome, position, plus the "allele sum" genotype
struct VariantData {
    std::string chrom;
    int pos;
    std::vector<int> genotype; // one element per sample (the sum of alleles, or -1 if missing)
};

// Stores the result of LD calculation
struct LDResult {
    double r;  // Correlation coefficient
    double r2; // Squared correlation coefficient
};

class VCFXHaplotypePhaser {
  public:
    // main runner
    int run(int argc, char *argv[]);

  private:
    // prints usage
    void displayHelp();

    // Main function that does phasing
    void phaseHaplotypes(std::istream &in, std::ostream &out, double ldThreshold);

    // Groups variants into haplotype blocks by naive r^2 threshold
    std::vector<std::vector<int>> groupVariants(const std::vector<VariantData> &variants, double ldThreshold);

    // calculates r^2 between two variants
    LDResult calculateLD(const VariantData &v1, const VariantData &v2);
};

#endif // VCFX_HAPLOTYPE_PHASER_H
