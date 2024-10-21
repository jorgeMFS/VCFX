#ifndef VCFX_HAPLOTYPE_PHASER_H
#define VCFX_HAPLOTYPE_PHASER_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// VCFXHaplotypePhaser: Header file for Haplotype-based Phasing Tool
class VCFXHaplotypePhaser {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Phases haplotypes in the VCF input
    void phaseHaplotypes(std::istream& in, std::ostream& out, double ldThreshold);

    // Groups variants into haplotype blocks based on linkage disequilibrium
    std::vector<std::vector<int>> groupVariants(const std::vector<std::vector<int>>& genotypeMatrix, double ldThreshold);

    // Calculates linkage disequilibrium (r²) between two variants
    double calculateLD(const std::vector<int>& variant1, const std::vector<int>& variant2);
};

#endif // VCFX_HAPLOTYPE_PHASER_H