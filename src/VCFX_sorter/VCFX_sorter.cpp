#include "VCFX_sorter.h"
#include <sstream>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_sorter\n"
              << "Usage: VCFX_sorter [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h            Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Sorts VCF records based on chromosome and position.\n\n"
              << "Example:\n"
              << "  ./VCFX_sorter < unsorted.vcf > sorted.vcf\n";
}

// Implement comparator based on chromosome and position
bool VCFRecord::operator<(const VCFRecord& other) const {
    if (chrom != other.chrom) {
        return chrom < other.chrom;
    }
    return pos < other.pos;
}

// Function to parse a VCF line into a VCFRecord
bool parseVCFLine(const std::string& line, VCFRecord& record) {
    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> fields;

    while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
    }

    if (fields.size() < 8) {
        return false; // Invalid VCF line
    }

    record.chrom = fields[0];
    try {
        record.pos = std::stoi(fields[1]);
    } catch (...) {
        return false; // Invalid position
    }
    record.id = fields[2];
    record.ref = fields[3];
    record.alt = fields[4];
    record.qual = fields[5];
    record.filter = fields[6];
    record.info = fields[7];

    // Parse sample fields if present
    if (fields.size() > 8) {
        record.samples.assign(fields.begin() + 8, fields.end());
    }

    return true;
}

// Function to sort VCF records
void sortVCFRecords(std::vector<VCFRecord>& records) {
    std::sort(records.begin(), records.end());
}

// Function to print sorted VCF records
void printSortedVCF(const std::vector<VCFRecord>& records, const std::string& header) {
    std::cout << header << "\n";
    for (const auto& record : records) {
        std::cout << record.chrom << "\t"
                  << record.pos << "\t"
                  << record.id << "\t"
                  << record.ref << "\t"
                  << record.alt << "\t"
                  << record.qual << "\t"
                  << record.filter << "\t"
                  << record.info;

        // Print sample fields if present
        if (!record.samples.empty()) {
            for (const auto& sample : record.samples) {
                std::cout << "\t" << sample;
            }
        }
        std::cout << "\n";
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

    std::string line;
    std::string header;
    std::vector<VCFRecord> records;

    while (std::getline(std::cin, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            header = line; // Save header
            continue;
        }

        VCFRecord record;
        if (parseVCFLine(line, record)) {
            records.push_back(record);
        } else {
            std::cerr << "Warning: Skipping invalid VCF line.\n";
        }
    }

    sortVCFRecords(records);
    printSortedVCF(records, header);

    return 0;
}
