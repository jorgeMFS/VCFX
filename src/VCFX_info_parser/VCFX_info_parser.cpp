#include "VCFX_info_parser.h"
#include <sstream>
#include <algorithm>
#include <unordered_map>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_info_parser\n"
              << "Usage: VCFX_info_parser [OPTIONS]\n\n"
              << "Options:\n"
              << "  --info, -i \"FIELD1,FIELD2\"   Specify the INFO fields to display (e.g., \"DP,AF\").\n"
              << "  --help, -h                    Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Parses the INFO field of a VCF file and displays the selected INFO fields in a user-friendly format.\n\n"
              << "Examples:\n"
              << "  ./VCFX_info_parser --info \"DP,AF\" < input.vcf > output_info.tsv\n";
}

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], std::vector<std::string>& info_fields) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--info" || arg == "-i") && i + 1 < argc) {
            std::string fields_str = argv[++i];
            // Split the fields by comma
            std::stringstream ss(fields_str);
            std::string field;
            while (std::getline(ss, field, ',')) {
                // Trim whitespace from field
                field.erase(0, field.find_first_not_of(" \t\n\r\f\v"));
                field.erase(field.find_last_not_of(" \t\n\r\f\v") + 1);
                if (!field.empty()) {
                    info_fields.push_back(field);
                }
            }
            return true;
        } else if (arg.find("--info=") == 0) {
            std::string fields_str = arg.substr(7);
            std::stringstream ss(fields_str);
            std::string field;
            while (std::getline(ss, field, ',')) {
                // Trim whitespace from field
                field.erase(0, field.find_first_not_of(" \t\n\r\f\v"));
                field.erase(field.find_last_not_of(" \t\n\r\f\v") + 1);
                if (!field.empty()) {
                    info_fields.push_back(field);
                }
            }
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

// Function to parse the INFO field and display selected fields
bool parseInfoFields(std::istream& in, std::ostream& out, const std::vector<std::string>& info_fields) {
    std::string line;
    bool header_printed = false;

    // Print header
    if (!info_fields.empty()) {
        out << "CHROM\tPOS\tID\tREF\tALT";
        for (const auto& field : info_fields) {
            out << "\t" << field;
        }
        out << "\n";
    }

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            // Skip header lines except #CHROM
            if (line.find("#CHROM") == 0) {
                // Optionally, you could include parsing and displaying header information
            }
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info)) {
            std::cerr << "Warning: Skipping invalid VCF line.\n";
            continue;
        }

        // Parse INFO field into key-value pairs
        std::vector<std::string> info_entries = split(info, ';');
        std::unordered_map<std::string, std::string> info_map;
        for (const auto& entry : info_entries) {
            size_t eq = entry.find('=');
            if (eq != std::string::npos) {
                std::string key = entry.substr(0, eq);
                std::string value = entry.substr(eq + 1);
                info_map[key] = value;
            } else {
                // Flag without a value
                info_map[entry] = "";
            }
        }

        // Output selected INFO fields
        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt;
        for (const auto& field : info_fields) {
            auto it = info_map.find(field);
            if (it != info_map.end()) {
                out << "\t" << it->second;
            } else {
                out << "\t.";
            }
        }
        out << "\n";
    }

    return true;
}

int main(int argc, char* argv[]) {
    std::vector<std::string> info_fields;

    // Argument parsing
    if (!parseArguments(argc, argv, info_fields)) {
        std::cerr << "Error: INFO fields not specified.\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    // Parse and display INFO fields
    bool success = parseInfoFields(std::cin, std::cout, info_fields);
    return success ? 0 : 1;
}
