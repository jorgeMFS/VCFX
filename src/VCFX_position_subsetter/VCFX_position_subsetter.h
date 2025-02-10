#ifndef VCFX_POSITION_SUBSETTER_H
#define VCFX_POSITION_SUBSETTER_H

#include <iostream>
#include <string>

class VCFXPositionSubsetter {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();

    // parse "chr1:10000-20000", store in regionChrom, regionStart, regionEnd
    bool parseRegion(const std::string &regionStr,
                     std::string &chrom,
                     int &start, int &end);

    // do the actual subsetting from in->out
    bool subsetVCFByPosition(std::istream &in, std::ostream &out,
                             const std::string &regionChrom,
                             int regionStart,
                             int regionEnd);
};

#endif // VCFX_POSITION_SUBSETTER_H
