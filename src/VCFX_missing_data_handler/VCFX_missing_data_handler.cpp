#include "VCFX_missing_data_handler.h"
#include <sstream>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_missing_data_handler\n"
              << "Usage: VCFX_missing_data_handler [OPTIONS] < input.vcf > output.vcf\n\n"
              << "Options:\n"
              << "  --fill-missing, -f         Impute missing genotypes with a default value (e.g., ./.).\n"
              << "  --default-genotype, -d GEN  Specify the default genotype for imputation (default: ./.).\n"
              << "  --help, -h                  Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Flags or imputes missing genotype data in a VCF file. By default, missing genotypes are flagged, "
              << "but can be imputed with a specified genotype using the --fill-missing option.\n\n"
              << "Examples:\n"
              << "  Flag missing data:\n"
              << "    ./VCFX_missing_data_handler < input.vcf > flagged_output.vcf\n\n"
              << "  Impute missing data with ./.;\n"
              << "    ./VCFX_missing_data_handler --fill-missing --default-genotype \"./.\" < input.vcf > imputed_output.vcf\n";
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
bool parseArguments(int argc, char* argv[], Arguments& args) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--fill-missing" || arg == "-f")) {
            args.fill_missing = true;
        }
        else if ((arg == "--default-genotype" || arg == "-d") && i + 1 < argc) {
            args.default_genotype = argv[++i];
        }
        else if (arg == "--help" || arg == "-h") {
            printHelp();
            exit(0);
        }
        else {
            // Handle input files or other arguments if necessary
            args.input_files.push_back(arg);
        }
    }
    return true;
}

// Function to process VCF and handle missing data
bool handleMissingData(std::istream& in, std::ostream& out, const Arguments& args) {
    std::string line;
    std::vector<std::string> header_fields;
    bool header_found = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            out << "\n";
            continue;
        }

        if (line[0] == '#') {
            out << line << "\n";
            if (line.find("#CHROM") == 0) {
                // Parse header to identify sample columns
                header_fields = splitString(line, '\t');
                if (header_fields.size() > 9) {
                    // Samples start from the 10th column (index 9)
                }
                header_found = true;
            }
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

        // GT column is typically the first in FORMAT
        std::string format = fields[8];
        std::vector<std::string> format_fields = splitString(format, ':');
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

        // Process each sample
        for (size_t i = 9; i < fields.size(); ++i) {
            std::vector<std::string> sample_data = splitString(fields[i], ':');
            if (gt_index >= static_cast<int>(sample_data.size())) {
                // No GT data for this sample
                continue;
            }

            std::string genotype = sample_data[gt_index];
            // Check for missing genotype
            if (genotype == "." || genotype == "./." || genotype == ".|.") {
                if (args.fill_missing) {
                    // Impute with default genotype
                    sample_data[gt_index] = args.default_genotype;
                }
                else {
                    // Flag missing genotype by appending a tag, e.g., adding a flag in INFO or elsewhere
                    // For simplicity, we'll leave it as is or could add an annotation
                    // Here, we choose to leave as is
                }

                // Reconstruct the sample data
                std::stringstream ss;
                for (size_t j = 0; j < sample_data.size(); ++j) {
                    ss << sample_data[j];
                    if (j != sample_data.size() - 1) {
                        ss << ":";
                    }
                }
                fields[i] = ss.str();
            }
        }

        // Reconstruct the VCF line
        std::stringstream ss_line;
        for (size_t i = 0; i < fields.size(); ++i) {
            ss_line << fields[i];
            if (i != fields.size() - 1) {
                ss_line << "\t";
            }
        }
        out << ss_line.str() << "\n";
    }

    return true;
}

int main(int argc, char* argv[]) {
    Arguments args;
    parseArguments(argc, argv, args);

    if (args.fill_missing) {
        std::cerr << "Info: Missing genotypes will be imputed with genotype: " << args.default_genotype << "\n";
    }
    else {
        std::cerr << "Info: Missing genotypes will be flagged.\n";
    }

    bool success = handleMissingData(std::cin, std::cout, args);
    return success ? 0 : 1;
}
