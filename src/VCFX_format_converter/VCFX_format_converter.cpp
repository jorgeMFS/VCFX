#include "VCFX_format_converter.h"
#include <sstream>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_format_converter\n"
              << "Usage: VCFX_format_converter [OPTIONS]\n\n"
              << "Options:\n"
              << "  --to-bed             Convert VCF to BED format.\n"
              << "  --to-csv             Convert VCF to CSV format.\n"
              << "  --help, -h           Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Converts VCF files to specified formats (BED or CSV).\n\n"
              << "Example:\n"
              << "  ./VCFX_format_converter --to-bed < input.vcf > output.bed\n"
              << "  ./VCFX_format_converter --to-csv < input.vcf > output.csv\n";
}

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], OutputFormat& format) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--to-bed") {
            format = OutputFormat::BED;
            return true;
        } else if (arg == "--to-csv") {
            format = OutputFormat::CSV;
            return true;
        }
    }
    format = OutputFormat::UNKNOWN;
    return false;
}

// Function to convert VCF to BED
void convertVCFtoBED(std::istream& in, std::ostream& out) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip header lines
        }

        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string field;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() < 5) {
            continue; // Insufficient fields
        }

        std::string chrom = fields[0];
        int pos = std::stoi(fields[1]);
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        int end = pos + ref.length(); // Simple end position

        out << chrom << "\t" << pos - 1 << "\t" << end << "\t" << id << "\n";
    }
}

// Function to convert VCF to CSV
void convertVCFtoCSV(std::istream& in, std::ostream& out) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip header lines or handle accordingly
        }

        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string field;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        // Join fields with commas
        for (size_t i = 0; i < fields.size(); ++i) {
            // Handle fields with commas by enclosing in quotes
            if (fields[i].find(',') != std::string::npos) {
                out << "\"" << fields[i] << "\"";
            } else {
                out << fields[i];
            }

            if (i != fields.size() - 1) {
                out << ",";
            }
        }
        out << "\n";
    }
}

int main(int argc, char* argv[]) {
    // Argument parsing
    OutputFormat format;
    if (!parseArguments(argc, argv, format)) {
        std::cerr << "No valid output format specified.\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    // Handle help and invalid formats
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    switch (format) {
        case OutputFormat::BED:
            convertVCFtoBED(std::cin, std::cout);
            break;
        case OutputFormat::CSV:
            convertVCFtoCSV(std::cin, std::cout);
            break;
        default:
            std::cerr << "Unsupported output format.\n";
            std::cerr << "Use --help for usage information.\n";
            return 1;
    }

    return 0;
}
