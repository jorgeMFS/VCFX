#include "VCFX_validator.h"
#include <sstream>
#include <vector>
#include <cctype>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_validator\n"
              << "Usage: VCFX_validator [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h            Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Validates the integrity and format of a VCF file.\n\n"
              << "Example:\n"
              << "  ./VCFX_validator < input.vcf\n";
}

// Function to trim whitespace from both ends of a string
static inline std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");
    if (start == std::string::npos)
        return "";
    return s.substr(start, end - start + 1);
}

// Function to validate VCF meta-information headers
bool validateVCFHeader(const std::string& line) {
    // Basic validation: check if header starts with ##
    if (line.size() < 2 || line[0] != '#' || line[1] != '#') {
        return false;
    }

    // Further validations can be implemented as per VCF specifications
    // For simplicity, we'll assume headers starting with ## are valid
    return true;
}

// Function to validate a single VCF record
bool validateVCFRecord(const std::string& line, int line_number) {
    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> fields;

    // Split the line by tabs
    while (std::getline(ss, field, '\t')) {
        fields.push_back(trim(field));
    }

    // VCF records should have at least 8 fields
    if (fields.size() < 8) {
        std::cerr << "Error: Line " << line_number << " has fewer than 8 fields.\n";
        return false;
    }

    // Validate CHROM: non-empty
    if (fields[0].empty()) {
        std::cerr << "Error: Line " << line_number << " has empty CHROM field.\n";
        return false;
    }

    // Validate POS: positive integer
    try {
        int pos = std::stoi(fields[1]);
        if (pos <= 0) {
            std::cerr << "Error: Line " << line_number << " has invalid POS value.\n";
            return false;
        }
    } catch (...) {
        std::cerr << "Error: Line " << line_number << " has invalid POS value.\n";
        return false;
    }

    // REF: non-empty
    if (fields[3].empty()) {
        std::cerr << "Error: Line " << line_number << " has empty REF field.\n";
        return false;
    }

    // ALT: non-empty
    if (fields[4].empty()) {
        std::cerr << "Error: Line " << line_number << " has empty ALT field.\n";
        return false;
    }

    // QUAL: Either '.' or a non-negative floating-point number
    if (fields[5] != ".") {
        try {
            float qual = std::stof(fields[5]);
            if (qual < 0) {
                std::cerr << "Error: Line " << line_number << " has negative QUAL value.\n";
                return false;
            }
        } catch (...) {
            std::cerr << "Error: Line " << line_number << " has invalid QUAL value.\n";
            return false;
        }
    }

    // FILTER: non-empty (can be "PASS" or specific filter names)
    if (fields[6].empty()) {
        std::cerr << "Error: Line " << line_number << " has empty FILTER field.\n";
        return false;
    }

    // INFO: can be empty but should follow key=value format separated by semicolons
    // Basic check: if not empty, contains at least one '=' or contains flags (no '=')
    if (!fields[7].empty()) {
        bool has_key_value = false;
        bool has_flag = false;
        std::stringstream info_ss(fields[7]);
        std::string info_field;
        while (std::getline(info_ss, info_field, ';')) {
            if (info_field.find('=') != std::string::npos) {
                has_key_value = true;
            } else {
                has_flag = true;
            }
        }
        if (!has_key_value && !has_flag) {
            std::cerr << "Error: Line " << line_number << " has invalid INFO field.\n";
            return false;
        }
    }

    // FORMAT and sample fields can be further validated if necessary
    // This implementation focuses on basic validations

    return true;
}

// Function to validate the entire VCF file
bool validateVCF(std::istream& in, std::ostream& out) {
    std::string line;
    int line_number = 0;
    bool header_found = false;
   
    while (std::getline(in, line)) {
        line_number++;
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                header_found = true;
                continue; // Header line found, continue to records
            } else {
                // Validate meta-information headers
                if (!validateVCFHeader(line)) {
                    std::cerr << "Error: Invalid VCF meta-information header at line " << line_number << ".\n";
                    return false;
                }
            }
            continue; // Continue processing headers
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records at line " << line_number << ".\n";
            return false;
        }

        // Validate the VCF record
        if (!validateVCFRecord(line, line_number)) {
            // Detailed error already printed in validateVCFRecord
            return false;
        }
    }

    if (!header_found) {
        std::cerr << "Error: VCF header (#CHROM) not found in the file.\n";
        return false;
    }

    std::cout << "VCF file is valid.\n";
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

    // Validate VCF
    bool is_valid = validateVCF(std::cin, std::cout);
    return is_valid ? 0 : 1;
}