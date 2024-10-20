#include "VCFX_haplotype_extractor.h"
#include <sstream>
#include <vector>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_haplotype_extractor\n"
              << "Usage: VCFX_haplotype_extractor [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h                     Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Extracts haplotype blocks from phased genotype data in a VCF file.\n\n"
              << "Example:\n"
              << "  ./VCFX_haplotype_extractor < phased.vcf > haplotypes.tsv\n";
}

// Helper function to split a string by a delimiter
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::stringstream ss(s);
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Function to extract haplotype blocks from VCF
bool extractHaplotypeBlocks(std::istream& in, std::ostream& out) {
    std::string line;
    bool header_found = false;
    std::vector<std::string> sample_names;

    // Output header for haplotype blocks
    out << "CHROM\tSTART\tEND\tHAPLOTYPES\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                // Parse sample names
                std::vector<std::string> fields = split(line, '\t');
                if (fields.size() > 9) {
                    sample_names.assign(fields.begin() + 9, fields.end());
                }
                header_found = true;
            }
            continue; // Skip header lines
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        std::vector<std::string> fields = split(line, '\t');
        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping invalid VCF line with fewer than 9 fields.\n";
            continue;
        }

        std::string chrom = fields[0];
        int pos = 0;
        try {
            pos = std::stoi(fields[1]);
        } catch (...) {
            std::cerr << "Warning: Invalid POS value. Skipping line.\n";
            continue;
        }

        std::string format = fields[8];
        std::vector<std::string> format_fields = split(format, ':');
        int gt_index = -1;
        for (size_t i = 0; i < format_fields.size(); ++i) {
            if (format_fields[i] == "GT") {
                gt_index = static_cast<int>(i);
                break;
            }
        }

        if (gt_index == -1) {
            std::cerr << "Warning: GT field not found in FORMAT column. Skipping line.\n";
            continue;
        }

        std::vector<std::string> haplotypes;
        for (size_t i = 9; i < fields.size(); ++i) {
            std::vector<std::string> sample_data = split(fields[i], ':');
            if (static_cast<size_t>(gt_index) >= sample_data.size()) {
                haplotypes.push_back(".");
                continue; // No GT data for this sample
            }
            haplotypes.push_back(sample_data[gt_index]);
        }

        // For simplicity, consider each variant as a separate haplotype block
        // In practice, you might want to group contiguous variants into blocks
        out << chrom << "\t" << pos << "\t" << pos << "\t";

        // Concatenate haplotypes separated by commas
        std::string haplo_str;
        for (size_t i = 0; i < haplotypes.size(); ++i) {
            haplo_str += haplotypes[i];
            if (i != haplotypes.size() - 1) {
                haplo_str += ",";
            }
        }
        out << haplo_str << "\n";
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

    // Extract haplotype blocks
    bool success = extractHaplotypeBlocks(std::cin, std::cout);
    return success ? 0 : 1;
}
