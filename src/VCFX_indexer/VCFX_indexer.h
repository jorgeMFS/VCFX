#ifndef VCFX_INDEXER_H
#define VCFX_INDEXER_H

#include <iostream>
#include <string>

class VCFXIndexer {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();
    void createVCFIndex(std::istream &in, std::ostream &out);
};

#endif // VCFX_INDEXER_H
