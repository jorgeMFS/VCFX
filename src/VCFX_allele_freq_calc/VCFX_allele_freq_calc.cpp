#include "VCFX_allele_freq_calc.h"
#include <sstream>
#include <vector>
#include <iomanip>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_allele_freq_calc\n"
              << "Usage: VCFX_allele_freq_calc [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h  Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Calculates allele frequencies from VCF records and outputs them in a tab-separated format.\n\n"
              << "Example:\n"
              << "  ./VCFX_allele_freq_calc < input.vcf > allele_frequencies.tsv\n";
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

// Function to parse genotype string and count allele occurrences
void parseGenotype(const std::string& genotype, int& alt_count, int& total_count) {
    std::vector<std::string> alleles = split(genotype, '/');
    for (const auto& allele : alleles) {
        if (allele == "0") {
            total_count += 1;
        } else if (allele == "1") {
            alt_count += 1;
            total_count += 1;
        }
        // Extend this if there are more alleles (e.g., "2", "3", ...)
    }
}

// Function to perform allele frequency calculation on VCF records
void calculateAlleleFrequency(std::istream& in, std::ostream& out) {
    std::string line;
    std::string header;
    int sample_count = 0;

    // Print header for output
    out << "CHROM\tPOS\tID\tREF\tALT\tAllele_Frequency\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                // Determine the number of samples
                std::vector<std::string> fields = split(line, '\t');
                if (fields.size() > 9) {
                    sample_count = static_cast<int>(fields.size()) - 9;
                }
            }
            continue; // Skip header lines
        }

        std::vector<std::string> fields = split(line, '\t');
        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping invalid VCF line with fewer than 9 fields.\n";
            continue;
        }

        std::string chrom = fields[0];
        std::string pos = fields[1];
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        std::string format = fields[8];

        // Find the GT index in the FORMAT column
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

        int alt_count = 0;
        int total_count = 0;

        // Iterate over each sample to count alleles
        for (size_t i = 9; i < fields.size(); ++i) {
            std::vector<std::string> sample_data = split(fields[i], ':');
            if (static_cast<size_t>(gt_index) >= sample_data.size()) {
                continue; // No GT data for this sample
            }
            parseGenotype(sample_data[gt_index], alt_count, total_count);
        }

        // Calculate allele frequency
        double allele_freq = (total_count > 0) ? static_cast<double>(alt_count) / static_cast<double>(total_count) : 0.0;

        // Output the result with fixed decimal precision
        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t"
            << std::fixed << std::setprecision(4) << allele_freq << "\n";
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

    // Perform allele frequency calculation
    calculateAlleleFrequency(std::cin, std::cout);
    return 0;
}
