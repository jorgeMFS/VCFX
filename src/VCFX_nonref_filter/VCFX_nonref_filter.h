#ifndef VCFX_NONREF_FILTER_H
#define VCFX_NONREF_FILTER_H

#include <iostream>
#include <string>
#include <vector>

/*
 * A fast tool for filtering out any variant whose genotypes are all
 * homozygous reference in every sample. Multi-ploid aware: a genotype
 * is considered homozygous reference if all alleles are '0'.
 */
class VCFXNonRefFilter {
public:
    // entry point with arguments
    int run(int argc, char* argv[]);

private:
    // prints usage
    void displayHelp();

    // does the actual filter: if every sample is hom-ref => skip
    void filterNonRef(std::istream& in, std::ostream& out);

    // parse genotype subfield; returns true if definitely hom-ref
    // returns false if any allele is not '0' or if missing
    bool isHomRef(const std::string &genotypeString) const;
};

#endif // VCFX_NONREF_FILTER_H
