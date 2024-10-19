#include "VCFX_sorter.h"
#include <sstream>
#include <algorithm>

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
    record.pos = std::stoi(fields[1]);
    record.id = fields[2];
    record.ref = fields[3];
    record.alt = fields[4];
    record.qual = fields[5];
    record.filter = fields[6];
    record.info = fields[7];

    // Extract sample columns if present
    record.samples.clear();
    if (fields.size() > 8) {
        for (size_t i = 8; i < fields.size(); ++i) {
            record.samples.push_back(fields[i]);
        }
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

        for (const auto& sample : record.samples) {
            std::cout << "\t" << sample;
        }
        std::cout << "\n";
    }
}

int main(int argc, char* argv[]) {
    std::vector<VCFRecord> records;
    std::string line;
    std::string header;

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
