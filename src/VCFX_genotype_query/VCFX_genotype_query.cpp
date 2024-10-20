#include "VCFX_genotype_query.h"
#include <sstream>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_genotype_query\n"
              << "Usage: VCFX_genotype_query [OPTIONS]\n\n"
              << "Options:\n"
              << "  --genotype-query, -g \"GENOTYPE\"  Specify the genotype to query (e.g., \"0/1\", \"1/1\").\n"
              << "  --help, -h                        Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Filters VCF records based on the specified genotype query. Only records matching the genotype "
              << "criteria will be outputted.\n\n"
              << "Example:\n"
              << "  ./VCFX_genotype_query --genotype-query \"0/1\" < input.vcf > output.vcf\n";
}

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], std::string& genotype_query) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--genotype-query" || arg == "-g") && i + 1 < argc) {
            genotype_query = argv[++i];
            return true;
        } else if (arg.find("--genotype-query=") == 0) {
            genotype_query = arg.substr(17);
            return true;
        } else if (arg == "--help" || arg == "-h") {
            printHelp();
            exit(0);
        }
    }
    return false;
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

// Function to perform genotype query on VCF records
void genotypeQuery(std::istream& in, std::ostream& out, const std::string& genotype_query) {
    std::string line;
    std::string header;
    std::vector<std::string> sample_names;
    bool header_found = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                header = line;
                std::vector<std::string> fields = split(line, '\t');
                // Samples start from the 10th column (index 9)
                if (fields.size() > 9) {
                    sample_names.assign(fields.begin() + 9, fields.end());
                }
                header_found = true;
            }
            out << line << "\n"; // Print header lines
            continue;
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return;
        }

        std::vector<std::string> fields = split(line, '\t');
        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping invalid VCF line with fewer than 9 fields.\n";
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

        bool match_found = false;
        for (size_t i = 9; i < fields.size(); ++i) {
            std::vector<std::string> sample_data = split(fields[i], ':');
            if (static_cast<size_t>(gt_index) >= sample_data.size()) {
                continue; // No GT data for this sample
            }
            if (sample_data[gt_index] == genotype_query) {
                match_found = true;
                break;
            }
        }

        if (match_found) {
            out << line << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    std::string genotype_query;
    if (!parseArguments(argc, argv, genotype_query)) {
        std::cerr << "Usage: " << argv[0] << " --genotype-query \"0/1\" < input.vcf > output.vcf\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    genotypeQuery(std::cin, std::cout, genotype_query);
    return 0;
}
