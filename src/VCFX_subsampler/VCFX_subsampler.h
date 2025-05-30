#ifndef VCFX_SUBSAMPLER_H
#define VCFX_SUBSAMPLER_H

#include <iostream>
#include <string>
#include <vector>

class VCFXSubsampler {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();

    // Reservoir sampling method
    void subsampleLines(std::istream &in, std::ostream &out, int sampleSize, unsigned int seed);
};

#endif
