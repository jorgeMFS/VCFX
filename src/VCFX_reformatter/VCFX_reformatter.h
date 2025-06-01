#ifndef VCFX_REFORMATTER_H
#define VCFX_REFORMATTER_H

#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

class VCFXReformatter {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();

    // Reformat the VCF from 'in' to 'out' using the user-specified lists
    void reformatVCF(std::istream &in, std::ostream &out, const std::vector<std::string> &compressInfoFields,
                     const std::vector<std::string> &compressFormatFields, const std::vector<std::string> &reorderInfo,
                     const std::vector<std::string> &reorderFormat);

    // Remove user-specified fields from the semicolon-based INFO
    // (Used for compressing info keys)
    std::string compressInfo(const std::string &infoStr, const std::unordered_set<std::string> &keysToRemove);

    // Remove user-specified keys from colon-based FORMAT
    // and from each genotype subfield in that position
    // Returns the new format string, plus an index mapping
    std::string compressFormat(const std::string &formatStr, const std::unordered_set<std::string> &keysToRemove,
                               std::vector<int> &keepIndices);

    // Reorder a semicolon-based INFO
    std::string reorderInfo(const std::string &infoStr, const std::vector<std::string> &order);

    // Reorder a colon-based FORMAT; returns new format, plus old->new index mapping
    std::string reorderFormat(const std::string &fmtStr, const std::vector<std::string> &order,
                              std::vector<int> &oldToNew);

    // reorder genotype subfields for each sample based on oldToNew
    // or remove subfields if the new index is negative
    std::string applyFormatReorderToSample(const std::string &sampleStr, const std::vector<int> &oldToNew);
};

#endif
