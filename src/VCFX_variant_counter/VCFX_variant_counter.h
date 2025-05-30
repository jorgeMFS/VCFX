#ifndef VCFX_VARIANT_COUNTER_H
#define VCFX_VARIANT_COUNTER_H

#include <iostream>
#include <string>

class VCFXVariantCounter {
  public:
    int run(int argc, char *argv[]);

  private:
    // If true, any line with <8 columns is a fatal error
    bool strictMode = false;

    // Show usage
    void displayHelp();

    // The actual counting function
    int countVariants(std::istream &in);
    int countVariantsGzip(std::istream &in);
    bool processLine(const std::string &line, int lineNumber, int &count);
};

#endif
