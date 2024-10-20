#include "VCFX_duplicate_remover.h"
#include <sstream>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_duplicate_remover\n"
              << "Usage: VCFX_duplicate_remover [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h               Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Removes duplicate variants from a VCF file based on chromosome, position, REF, and ALT alleles.\n\n"
              << "Examples:\n"
              << "  ./VCFX_duplicate_remover < input.vcf > unique_variants.vcf\n";
}

// Function to remove duplicate variants
bool removeDuplicates(std::istream& in, std::ostream& out) {
    std::string line;
    bool header_printed = false;
    std::unordered_set<VariantKey, VariantKeyHash> seen_variants;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            // Print header lines
            out << line << "\n";
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info)) {
            std::cerr << "Warning: Skipping invalid VCF line.\n";
            continue;
        }

        int pos_int = 0;
        try {
            pos_int = std::stoi(pos);
        } catch (...) {
            std::cerr << "Warning: Invalid POS value. Skipping line.\n";
            continue;
        }

        VariantKey key = {chrom, pos_int, ref, alt};
        if (seen_variants.find(key) == seen_variants.end()) {
            // Variant is unique
            seen_variants.insert(key);
            out << line << "\n";
        } else {
            // Duplicate variant found; skip
            // Optionally, you can log or count duplicates
        }
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

    // Remove duplicates
    bool success = removeDuplicates(std::cin, std::cout);
    return success ? 0 : 1;
}
