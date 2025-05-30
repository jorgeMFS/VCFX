#ifndef VCFX_REF_COMPARATOR_H
#define VCFX_REF_COMPARATOR_H

#include <iostream>
#include <string>
#include <unordered_map>

class VCFXRefComparator {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();
    bool loadReference(const std::string &referenceFastaPath);
    void compareVCF(std::istream &vcfIn, std::ostream &vcfOut);

    // chromosome -> uppercase sequence
    std::unordered_map<std::string, std::string> referenceGenome;

    // store whether we have already inserted the INFO line
    bool infoHeaderInserted = false;
};

#endif
