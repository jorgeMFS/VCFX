#ifndef VCFX_RECORD_FILTER_H
#define VCFX_RECORD_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// Enumeration for comparison operators
enum class Operator {
    GREATER_THAN,    // >
    LESS_THAN,       // <
    GREATER_EQUAL,   // >=
    LESS_EQUAL,      // <=
    EQUAL            // ==
};

// Structure to hold individual filter criteria
struct FilterCriterion {
    std::string field;
    Operator op;
    double value;
};

// Function to parse and populate filter criteria from a string
bool parseCriteria(const std::string& criteria_str, std::vector<FilterCriterion>& criteria);

// Function to apply all filter criteria to a single VCF record
bool applyFilters(const std::string& record, const std::vector<FilterCriterion>& criteria);

// Function to process and filter VCF records based on the provided criteria
void processRecords(std::istream& in, std::ostream& out, const std::vector<FilterCriterion>& criteria);

#endif // VCFX_RECORD_FILTER_H
