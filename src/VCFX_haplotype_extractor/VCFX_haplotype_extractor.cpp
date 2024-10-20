#include "VCFX_haplotype_extractor.h"
#include <sstream>
#include <algorithm>
#include <cctype>
#include <stdexcept>

// Constructor
HaplotypeExtractor::HaplotypeExtractor() {}

// Function to display help message
void printHelp() {
    std::cout << "VCFX_haplotype_extractor\n"
              << "Usage: VCFX_haplotype_extractor [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h                     Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Extracts phased haplotype blocks from genotype data in a VCF file. "
              << "It reconstructs haplotypes for each sample by analyzing phased genotype information, "
              << "grouping contiguous variants into haplotype blocks based on consistent phasing and proximity.\n\n"
              << "Examples:\n"
              << "  ./VCFX_haplotype_extractor < phased.vcf > haplotypes.tsv\n";
}

// Utility function to split a string by a delimiter
std::vector<std::string> HaplotypeExtractor::splitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string temp;
    while (std::getline(ss, temp, delimiter)) {
        tokens.push_back(temp);
    }
    return tokens;
}

// Parses the VCF header to extract sample names
bool HaplotypeExtractor::parseHeader(const std::string& headerLine) {
    std::vector<std::string> fields = splitString(headerLine, '\t');
    if (fields.size() <= 8) {
        std::cerr << "Error: VCF header does not contain sample columns.\n";
        return false;
    }

    // Extract sample names starting from the 9th column
    sampleNames.assign(fields.begin() + 9, fields.end());
    numSamples = sampleNames.size();
    return true;
}

// Validates if all samples are phased
bool HaplotypeExtractor::areAllSamplesPhased(const std::vector<std::string>& genotypeFields) {
    for (const auto& gt : genotypeFields) {
        // A phased genotype contains '|', e.g., 1
        if (gt.find('|') == std::string::npos) {
            return false;
        }
    }
    return true;
}

// Processes a single VCF variant and updates haplotype blocks
bool HaplotypeExtractor::processVariant(const std::vector<std::string>& fields, std::vector<HaplotypeBlock>& haplotypeBlocks, int currentPos) {
    if (fields.size() < 9) {
        std::cerr << "Warning: Skipping invalid VCF line with fewer than 9 fields.\n";
        return false;
    }

    std::string chrom = fields[0];
    int pos = 0;
    try {
        pos = std::stoi(fields[1]);
    } catch (...) {
        std::cerr << "Warning: Invalid POS value. Skipping line.\n";
        return false;
    }

    std::string format = fields[8];
    std::vector<std::string> formatFields = splitString(format, ':');
    int gtIndex = -1;
    for (size_t i = 0; i < formatFields.size(); ++i) {
        if (formatFields[i] == "GT") {
            gtIndex = static_cast<int>(i);
            break;
        }
    }

    if (gtIndex == -1) {
        std::cerr << "Warning: GT field not found in FORMAT column. Skipping line.\n";
        return false;
    }

    std::vector<std::string> genotypeFields;
    bool allPhased = true;
    for (size_t i = 9; i < fields.size(); ++i) {
        std::vector<std::string> sampleData = splitString(fields[i], ':');
        if (static_cast<size_t>(gtIndex) >= sampleData.size()) {
            genotypeFields.push_back("./.");
            allPhased = false;
            continue; // No GT data for this sample
        }
        std::string gt = sampleData[gtIndex];
        genotypeFields.push_back(gt);
        if (gt.find('|') == std::string::npos) {
            allPhased = false;
        }
    }

    if (!allPhased) {
        std::cerr << "Warning: Not all samples are phased at position " << pos << " on " << chrom << ". Skipping haplotype extraction for this variant.\n";
        return false;
    }

    // Update haplotype blocks
    updateHaplotypeBlocks(genotypeFields, haplotypeBlocks, pos, chrom);

    return true;
}

// Updates haplotype blocks with new variant information
void HaplotypeExtractor::updateHaplotypeBlocks(const std::vector<std::string>& genotypeFields, std::vector<HaplotypeBlock>& haplotypeBlocks, int pos, const std::string& chrom) {
    if (haplotypeBlocks.empty()) {
        // Initialize the first haplotype block
        HaplotypeBlock block;
        block.chrom = chrom;
        block.start = pos;
        block.end = pos;
        for (const auto& gt : genotypeFields) {
            block.haplotypes.push_back(gt);
        }
        haplotypeBlocks.push_back(block);
        return;
    }

    // Get the last haplotype block
    HaplotypeBlock& lastBlock = haplotypeBlocks.back();

    // Determine if the current variant should be part of the last block
    // Criteria:
    // 1. Same chromosome
    // 2. Distance between variants is below a threshold (e.g., 100 kb)
    // 3. Phasing is consistent (already ensured)

    const int distanceThreshold = 100000; // 100 kb
    if (chrom == lastBlock.chrom && (pos - lastBlock.end) <= distanceThreshold) {
        // Extend the last haplotype block
        lastBlock.end = pos;
        for (size_t i = 0; i < genotypeFields.size(); ++i) {
            lastBlock.haplotypes[i] += "|" + genotypeFields[i];
        }
    } else {
        // Start a new haplotype block
        HaplotypeBlock newBlock;
        newBlock.chrom = chrom;
        newBlock.start = pos;
        newBlock.end = pos;
        for (const auto& gt : genotypeFields) {
            newBlock.haplotypes.push_back(gt);
        }
        haplotypeBlocks.push_back(newBlock);
    }
}

// Parses the VCF file and extracts haplotype blocks
bool HaplotypeExtractor::extractHaplotypes(std::istream& in, std::ostream& out) {
    std::string line;
    bool header_found = false;
    std::vector<HaplotypeBlock> haplotypeBlocks;

    // Read the VCF file line by line
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                // Parse header to extract sample names
                if (!parseHeader(line)) {
                    return false;
                }
                header_found = true;
            }
            continue; // Skip header lines
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        // Split the VCF record into fields
        std::vector<std::string> fields = splitString(line, '\t');

        // Process the variant to extract haplotypes
        if (!processVariant(fields, haplotypeBlocks, 0)) {
            continue; // Skip variants that couldn't be processed
        }
    }

    // Output the haplotype blocks
    out << "CHROM\tSTART\tEND\t" ;
    for (const auto& sample : sampleNames) {
        out << sample << "\t";
    }
    out << "\n";

    for (const auto& block : haplotypeBlocks) {
        out << block.chrom << "\t" << block.start << "\t" << block.end << "\t";
        for (size_t i = 0; i < block.haplotypes.size(); ++i) {
            out << block.haplotypes[i];
            if (i != block.haplotypes.size() - 1) {
                out << "\t";
            }
        }
        out << "\n";
    }

    return true;
}

int main(int argc, char* argv[]) {
    // Argument parsing for help
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // Instantiate the haplotype extractor
    HaplotypeExtractor extractor;

    // Extract haplotype blocks from stdin and output to stdout
    bool success = extractor.extractHaplotypes(std::cin, std::cout);
    return success ? 0 : 1;
}
