#ifndef VCFX_HAPLOTYPE_EXTRACTOR_H
#define VCFX_HAPLOTYPE_EXTRACTOR_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// Function to display help message
void printHelp();

// Structure to represent a haplotype block
struct HaplotypeBlock {
    std::string chrom;
    int start;
    int end;
    std::vector<std::string> haplotypes; // One haplotype per sample
};

// Class to handle haplotype extraction
class HaplotypeExtractor {
public:
    HaplotypeExtractor();
    ~HaplotypeExtractor() = default;

    // Parses the VCF file and extracts haplotype blocks
    bool extractHaplotypes(std::istream& in, std::ostream& out);

private:
    std::vector<std::string> sampleNames;
    size_t numSamples = 0;

    // Parses the VCF header to extract sample names
    bool parseHeader(const std::string& headerLine);

    // Splits a string by a delimiter
    std::vector<std::string> splitString(const std::string& str, char delimiter);

    // Reconstructs haplotype blocks based on phased genotypes
    bool processVariant(const std::vector<std::string>& fields, std::vector<HaplotypeBlock>& haplotypeBlocks, int currentPos);

    // Validates if all samples are phased
    bool areAllSamplesPhased(const std::vector<std::string>& genotypeFields);

    // Updates haplotype blocks with new variant information
    void updateHaplotypeBlocks(const std::vector<std::string>& genotypeFields, std::vector<HaplotypeBlock>& haplotypeBlocks, int pos, const std::string& chrom);
};

#endif // VCFX_HAPLOTYPE_EXTRACTOR_H
