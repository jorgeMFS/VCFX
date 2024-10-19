#ifndef VCFX_SUBSAMPLER_H
#define VCFX_SUBSAMPLER_H

#include <iostream>
#include <string>

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], int& sample_size);

// Function to perform reservoir sampling on VCF records
void subsampleVariants(std::istream& in, std::ostream& out, int sample_size);

#endif // VCFX_SUBSAMPLER_H
