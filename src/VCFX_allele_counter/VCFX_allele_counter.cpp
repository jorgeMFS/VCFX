#include "VCFX_allele_counter.h"
#include <sstream>
#include <algorithm>
#include <cctype>
#include <unordered_map>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_allele_counter\n"
              << "Usage: VCFX_allele_counter [OPTIONS] < input.vcf > allele_counts.tsv\n\n"
              << "Options:\n"
              << "  --samples, -s \"Sample1 Sample2\"   Specify the sample names to include in allele counting. If not specified, all samples are processed.\n"
              << "  --help, -h                        Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Counts the occurrences of reference and alternate alleles across specified samples for each variant in a VCF file.\n\n"
              << "Examples:\n"
              << "  Count alleles for SampleA and SampleB:\n"
              << "    ./VCFX_allele_counter --samples \"SampleA SampleB\" < input.vcf > allele_counts.tsv\n\n"
              << "  Count alleles for all samples:\n"
              << "    ./VCFX_allele_counter < input.vcf > allele_counts_all.tsv\n";
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
bool parseArguments(int argc, char* argv[], AlleleCounterArguments& args) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--samples" || arg == "-s") && i + 1 < argc) {
            std::string samples_str = argv[++i];
            args.samples = splitString(samples_str, ' ');
            // Trim whitespace from each sample name
            for (auto& sample : args.samples) {
                sample.erase(0, sample.find_first_not_of(" \t\n\r\f\v"));
                sample.erase(sample.find_last_not_of(" \t\n\r\f\v") + 1);
            }
        }
        else if (arg == "--help" || arg == "-h") {
            printHelp();
            exit(0);
        }
        else {
            std::cerr << "Warning: Unrecognized argument '" << arg << "'.\n";
        }
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

// Function to process VCF and count alleles
bool countAlleles(std::istream& in, std::ostream& out, const AlleleCounterArguments& args) {
    std::string line;
    std::vector<std::string> header_fields;
    bool header_found = false;
    std::vector<int> sample_indices;

    // Output header
    out << "CHROM\tPOS\tID\tREF\tALT\tSample\tRef_Count\tAlt_Count\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                header_fields = splitString(line, '\t');
                // Identify sample columns
                std::unordered_map<std::string, int> sample_map;
                for (size_t i = 9; i < header_fields.size(); ++i) {
                    sample_map[header_fields[i]] = static_cast<int>(i);
                }

                if (!args.samples.empty()) {
                    for (const auto& sample : args.samples) {
                        if (sample_map.find(sample) == sample_map.end()) {
                            std::cerr << "Error: Sample '" << sample << "' not found in VCF header.\n";
                            return false;
                        }
                        sample_indices.push_back(sample_map[sample]);
                    }
                }
                else {
                    // If no samples specified, include all
                    for (size_t i = 9; i < header_fields.size(); ++i) {
                        sample_indices.push_back(static_cast<int>(i));
                    }
                }
                header_found = true;
            }
            out << line << "\n"; // Retain header in output if needed
            continue;
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        std::vector<std::string> fields = splitString(line, '\t');
        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping invalid VCF line with fewer than 9 fields.\n";
            continue;
        }

        std::string chrom = fields[0];
        std::string pos = fields[1];
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];

        // Iterate over specified samples
        for (const auto& idx : sample_indices) {
            if (idx >= static_cast<int>(fields.size())) {
                std::cerr << "Warning: Sample index out of range. Skipping sample.\n";
                continue;
            }

            std::string sample_data = fields[idx];
            std::vector<std::string> genotype_fields = splitString(sample_data, ':');
            if (genotype_fields.empty()) {
                continue;
            }

            std::string genotype = extractGenotype(genotype_fields[0]);
            if (genotype.empty() || genotype == "./." || genotype == ".|.") {
                // Skip missing genotypes
                continue;
            }

            // Split genotype by '/' or '|'
            char sep = '/';
            size_t sep_pos = genotype.find('/');
            if (sep_pos == std::string::npos) {
                sep = '|';
                sep_pos = genotype.find('|');
            }

            std::vector<std::string> alleles;
            if (sep_pos != std::string::npos) {
                std::string allele1 = genotype.substr(0, sep_pos);
                std::string allele2 = genotype.substr(sep_pos + 1);
                alleles.push_back(allele1);
                alleles.push_back(allele2);
            }
            else {
                alleles.push_back(genotype);
            }

            int ref_count = 0;
            int alt_count = 0;

            for (const auto& allele : alleles) {
                if (allele == "0") {
                    ref_count++;
                }
                else if (!allele.empty() && allele != ".") {
                    alt_count++;
                }
            }

            // Retrieve sample name
            std::string sample_name = header_fields[idx];

            out << chrom << "\t"
                << pos << "\t"
                << id << "\t"
                << ref << "\t"
                << alt << "\t"
                << sample_name << "\t"
                << ref_count << "\t"
                << alt_count << "\n";
        }
    }

    return true;
}

int main(int argc, char* argv[]) {
    AlleleCounterArguments args;
    parseArguments(argc, argv, args);

    if (!args.samples.empty()) {
        std::cerr << "Info: Counting alleles for specified samples.\n";
    }
    else {
        std::cerr << "Info: Counting alleles for all samples.\n";
    }

    bool success = countAlleles(std::cin, std::cout, args);
    return success ? 0 : 1;
}
