#ifndef VCFX_INDEXER_H
#define VCFX_INDEXER_H

#include <iostream>
#include <string>

// A small tool that reads a VCF from stdin, and outputs a 3-column index:
//   CHROM   POS   FILE_OFFSET
// for each data line. FILE_OFFSET is a byte offset from the start of the file.
class VCFXIndexer {
public:
    int run(int argc, char* argv[]);
private:
    void displayHelp();

    // The core function that reads from 'in' and writes the index to 'out'
    void createVCFIndex(std::istream &in, std::ostream &out);
};

#endif // VCFX_INDEXER_H
