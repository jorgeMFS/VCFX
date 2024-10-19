#ifndef VCFX_SORTER_H
#define VCFX_SORTER_H

#include <iostream>
#include <string>
#include <vector>

// Structure to represent a VCF record
struct VCFRecord {
    std::string chrom;
    int pos;
    std::string id;
    std::string ref;
    std::string alt;
    std::string qual;
    std::string filter;
    std::string info;
    std::vector<std::string> samples;

    // Comparator for sorting
    bool operator<(const VCFRecord& other) const;
};

// Function to parse a VCF line into a VCFRecord
bool parseVCFLine(const std::string& line, VCFRecord& record);

// Function to sort VCF records
void sortVCFRecords(std::vector<VCFRecord>& records);

// Function to print sorted VCF records
void printSortedVCF(const std::vector<VCFRecord>& records, const std::string& header);

#endif // VCFX_SORTER_H
