#ifndef VCFX_SAMPLE_EXTRACTOR_H
#define VCFX_SAMPLE_EXTRACTOR_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>

class VCFXSampleExtractor {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();

    // read the user’s --samples list => store them
    // read VCF from in => output a valid VCF with only those samples
    void extractSamples(std::istream &in, std::ostream &out,
                        const std::vector<std::string> &samples);

};

#endif
