#include "VCFX_distance_calculator.h"
#include <sstream>
#include <unordered_map>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_distance_calculator\n"
              << "Usage: VCFX_distance_calculator [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h               Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Calculates the distance between consecutive variants along each chromosome in a VCF file.\n\n"
              << "Examples:\n"
              << "  ./VCFX_distance_calculator < input.vcf > variant_distances.tsv\n";
}

// Function to parse a VCF line into VCFVariant structure
bool parseVCFLine(const std::string& line, VCFVariant& variant) {
    if (line.empty() || line[0] == '#') return false;

    std::stringstream ss(line);
    std::string chrom, pos_str;

    if (!std::getline(ss, chrom, '\t') ||
        !std::getline(ss, pos_str, '\t')) {
        return false;
    }

    try {
        variant.chrom = chrom;
        variant.pos = std::stoi(pos_str);
    } catch (...) {
        return false;
    }

    return true;
}

// Function to calculate distances between consecutive variants
bool calculateDistances(std::istream& in, std::ostream& out) {
    std::string line;
    bool header_found = false;

    // Map to store the last position for each chromosome
    std::unordered_map<std::string, int> last_pos_map;

    // Output header
    out << "CHROM\tPOS\tPREV_POS\tDISTANCE\n";

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                header_found = true;
            }
            continue; // Skip header lines
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        VCFVariant variant;
        if (!parseVCFLine(line, variant)) {
            std::cerr << "Warning: Skipping malformed VCF line: " << line << "\n";
            continue;
        }

        if (last_pos_map.find(variant.chrom) != last_pos_map.end()) {
            int distance = variant.pos - last_pos_map[variant.chrom];
            out << variant.chrom << "\t" << variant.pos << "\t" << last_pos_map[variant.chrom] << "\t" << distance << "\n";
        } else {
            out << variant.chrom << "\t" << variant.pos << "\tNA\tNA\n";
        }

        // Update last position
        last_pos_map[variant.chrom] = variant.pos;
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

    // Calculate variant distances
    bool success = calculateDistances(std::cin, std::cout);
    return success ? 0 : 1;
}
