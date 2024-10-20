#include "VCFX_concordance_checker.h"
#include <sstream>
#include <algorithm>
#include <cctype>
#include <unordered_map>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_concordance_checker\n"
              << "Usage: VCFX_concordance_checker [OPTIONS] < input.vcf > concordance_report.tsv\n\n"
              << "Options:\n"
              << "  --samples, -s \"Sample1 Sample2\"  Specify the two sample names to compare.\n"
              << "  --help, -h                        Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Compares genotypes between two specified samples in a VCF file and calculates concordance metrics.\n\n"
              << "Examples:\n"
              << "  ./VCFX_concordance_checker --samples \"SampleA SampleB\" < input.vcf > concordance_report.tsv\n";
}

// Utility function to split a string by a delimiter
std::vector<std::string> splitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string temp;
    while (std::getline(ss, temp, delimiter)) {
        tokens.push_back(temp);
    }
    return tokens;
}

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], ConcordanceArguments& args) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--samples" || arg == "-s") && i + 1 < argc) {
            std::string samples_str = argv[++i];
            std::vector<std::string> samples = splitString(samples_str, ' ');
            if (samples.size() != 2) {
                std::cerr << "Error: Please specify exactly two sample names.\n";
                return false;
            }
            args.sample1 = samples[0];
            args.sample2 = samples[1];
        }
        else if (arg == "--help" || arg == "-h") {
            printHelp();
            exit(0);
        }
        else {
            std::cerr << "Warning: Unrecognized argument '" << arg << "'.\n";
        }
    }

    if (args.sample1.empty() || args.sample2.empty()) {
        std::cerr << "Error: Two sample names must be specified using --samples or -s option.\n";
        std::cerr << "Use --help for usage information.\n";
        return false;
    }

    return true;
}

// Function to extract genotype from genotype string
std::string extractGenotype(const std::string& genotype_str) {
    if (genotype_str.empty()) return "./.";

    std::vector<std::string> alleles = splitString(genotype_str, ':');
    if (!alleles.empty()) {
        return alleles[0];
    }
    return "./.";
}

// Function to process VCF and calculate concordance
bool calculateConcordance(std::istream& in, std::ostream& out, const ConcordanceArguments& args) {
    std::string line;
    std::vector<std::string> header_fields;
    bool header_found = false;
    int sample1_index = -1;
    int sample2_index = -1;

    // Concordance metrics
    int total_variants = 0;
    int concordant = 0;
    int discordant = 0;

    // Header for output
    out << "Sample1\tSample2\tConcordance\n";

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                header_fields = splitString(line, '\t');
                // Identify sample columns
                std::unordered_map<std::string, int> sample_map;
                for (size_t i = 9; i < header_fields.size(); ++i) {
                    sample_map[header_fields[i]] = static_cast<int>(i);
                }

                if (sample_map.find(args.sample1) == sample_map.end()) {
                    std::cerr << "Error: Sample '" << args.sample1 << "' not found in VCF header.\n";
                    return false;
                }
                if (sample_map.find(args.sample2) == sample_map.end()) {
                    std::cerr << "Error: Sample '" << args.sample2 << "' not found in VCF header.\n";
                    return false;
                }

                sample1_index = sample_map[args.sample1];
                sample2_index = sample_map[args.sample2];
                header_found = true;
            }
            continue;
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        std::vector<std::string> fields = splitString(line, '\t');
        if (fields.size() <= std::max(sample1_index, sample2_index)) {
            std::cerr << "Warning: Skipping line with insufficient columns.\n";
            continue;
        }

        std::string gt1 = extractGenotype(fields[sample1_index]);
        std::string gt2 = extractGenotype(fields[sample2_index]);

        if (gt1.empty() || gt2.empty()) {
            // Skip variants with missing genotypes in either sample
            continue;
        }

        total_variants++;

        // Compare genotypes
        if (gt1 == gt2) {
            concordant++;
        }
        else {
            discordant++;
        }

        // Output per-variant concordance
        std::string concordance_status = (gt1 == gt2) ? "Concordant" : "Discordant";
        out << args.sample1 << "\t" << args.sample2 << "\t" << concordance_status << "\n";
    }

    // Summary statistics
    std::cerr << "Total Variants Compared: " << total_variants << "\n";
    std::cerr << "Concordant Genotypes: " << concordant << "\n";
    std::cerr << "Discordant Genotypes: " << discordant << "\n";

    return true;
}

int main(int argc, char* argv[]) {
    ConcordanceArguments args;
    if (!parseArguments(argc, argv, args)) {
        return 1;
    }

    bool success = calculateConcordance(std::cin, std::cout, args);
    return success ? 0 : 1;
}
