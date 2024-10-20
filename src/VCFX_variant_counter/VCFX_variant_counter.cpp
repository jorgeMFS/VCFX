#include "VCFX_variant_counter.h"
#include <string>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_variant_counter\n"
              << "Usage: VCFX_variant_counter [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h            Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Counts the total number of variants in a VCF file.\n\n"
              << "Example:\n"
              << "  ./VCFX_variant_counter < input.vcf > variant_count.txt\n";
}

// Function to count variants in VCF
int countVariants(std::istream& in) {
    int count = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip header lines
        }
        count++;
    }
    return count;
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

    // Count variants
    int total_variants = countVariants(std::cin);
    std::cout << "Total Variants: " << total_variants << std::endl;
    return 0;
}
