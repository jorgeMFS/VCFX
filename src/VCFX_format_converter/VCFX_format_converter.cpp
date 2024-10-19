#include "VCFX_format_converter.h"
#include <sstream>
#include <algorithm>

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
    bool header_written = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (!header_written) {
                // Convert VCF header to CSV header
                std::cout << line << "\n"; // Optionally, skip or transform header
            }
            continue; // Skip header lines
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
    OutputFormat format;
    if (!parseArguments(argc, argv, format)) {
        std::cerr << "Usage: " << argv[0] << " --to-bed|--to-csv < input.vcf > output" << std::endl;
        return 1;
    }

    switch (format) {
        case OutputFormat::BED:
            convertVCFtoBED(std::cin, std::cout);
            break;
        case OutputFormat::CSV:
            convertVCFtoCSV(std::cin, std::cout);
            break;
        default:
            std::cerr << "Unsupported output format." << std::endl;
            return 1;
    }

    return 0;
}
