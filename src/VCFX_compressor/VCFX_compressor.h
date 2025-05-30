#ifndef VCFX_COMPRESSOR_H
#define VCFX_COMPRESSOR_H

#include <iostream>
#include <string>

// Function to display help message
void printHelp();

// Function to perform compression or decompression
bool compressDecompressVCF(std::istream &in, std::ostream &out, bool compress);

#endif // VCFX_COMPRESSOR_H
