#ifndef VCFX_HWE_TESTER_H
#define VCFX_HWE_TESTER_H

#include <iostream>
#include <string>
#include <vector>

class VCFXHWETester {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();
    void performHWE(std::istream &in);
    bool isBiallelic(const std::string &alt);
    bool parseGenotypes(const std::vector<std::string> &genotypes, int &homRef, int &het, int &homAlt);
    double calculateHWE(int homRef, int het, int homAlt);
    double genotypeProbability(int homRef, int het, int homAlt);
};

#endif // VCFX_HWE_TESTER_H
