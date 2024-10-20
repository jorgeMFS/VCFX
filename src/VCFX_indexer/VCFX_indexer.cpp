#include "VCFX_indexer.h"
#include <sstream>
#include <vector>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_indexer\n"
              << "Usage: VCFX_indexer [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h  Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Creates an index for a VCF file, mapping each variant's chromosome and position to its byte offset in the file.\n\n"
              << "Example:\n"
              << "  ./VCFX_indexer < input.vcf > index.tsv\n";
}

// Function to split a string by a delimiter and return a vector of tokens
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::stringstream ss(s);
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Function to create an index for VCF file
void createVCFIndex(std::istream& in, std::ostream& out) {
    std::string line;
    std::string header;
    long offset = 0;
    bool header_found = false;

    // Print header for index
    out << "CHROM\tPOS\tFILE_OFFSET\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            offset += (line.length() + 1); // +1 for newline character
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                header_found = true;
            }
            offset += (line.length() + 1);
            continue; // Skip header lines
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return;
        }

        std::vector<std::string> fields = split(line, '\t');
        if (fields.size() < 2) {
            std::cerr << "Warning: Skipping invalid VCF line with fewer than 2 fields.\n";
            offset += (line.length() + 1);
            continue;
        }

        std::string chrom = fields[0];
        int pos = 0;
        try {
            pos = std::stoi(fields[1]);
        } catch (...) {
            std::cerr << "Warning: Invalid POS value on line with CHROM " << chrom << ". Skipping line.\n";
            offset += (line.length() + 1);
            continue;
        }

        // Output the index entry
        out << chrom << "\t" << pos << "\t" << offset << "\n";

        // Update the offset
        offset += (line.length() + 1);
    }
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

    // Display help if no input is provided
    if (argc == 1) {
        printHelp();
        return 1;
    }

    // Create VCF index
    createVCFIndex(std::cin, std::cout);
    return 0;
}
