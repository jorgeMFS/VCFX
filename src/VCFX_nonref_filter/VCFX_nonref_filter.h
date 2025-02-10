#ifndef VCFX_NONREF_FILTER_H
#define VCFX_NONREF_FILTER_H

#include <iostream>
#include <string>
#include <vector>

/*
 * VCFX_nonref_filter:
 *   Reads a VCF, discards any variant for which EVERY sample is hom-ref
 *   (all alleles = '0'). If a genotype is missing or partial, treat that
 *   as “not guaranteed hom-ref,” so we keep that variant.
 */
class VCFXNonRefFilter {
public:
    // main entry point
    int run(int argc, char* argv[]);

private:
    // prints usage
    void displayHelp();

    // checks the input line-by-line, printing lines that pass the filter
    void filterNonRef(std::istream& in, std::ostream& out);

    // parse a sample’s genotype from a field => returns true if definitely hom-ref
    bool isDefinitelyHomRef(const std::string &genotypeField) const;
};

#endif // VCFX_NONREF_FILTER_H
