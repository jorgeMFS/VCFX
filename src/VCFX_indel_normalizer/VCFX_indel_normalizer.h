#ifndef VCFX_INDEL_NORMALIZER_H
#define VCFX_INDEL_NORMALIZER_H

#include <iostream>
#include <string>
#include <vector>

// A tool for normalizing INDELs (and any variant) to a minimal left-aligned representation
// without requiring an external reference genome. This includes:
//   - Splitting multi-ALT lines into separate lines (one ALT each).
//   - Removing the largest possible shared leading prefix, adjusting POS accordingly.
//   - Removing the largest possible shared trailing suffix.
class VCFXIndelNormalizer {
public:
    int run(int argc, char* argv[]);

private:
    // Print usage
    void displayHelp();

    // The main function: read VCF from 'in', write normalized lines to 'out'
    void normalizeIndels(std::istream& in, std::ostream& out);

    // For each ALT allele, produce a separate line. Then do left+right trim
    bool normalizeVariant(std::string &chrom, int &posInt,
                          std::string &ref,
                          std::string &alt);

    // checks if line is a variant line (#CHROM line => we pass it as header)
};

#endif // VCFX_INDEL_NORMALIZER_H
