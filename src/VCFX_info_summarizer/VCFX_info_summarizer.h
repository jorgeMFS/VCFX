#ifndef VCFX_INFO_SUMMARIZER_H
#define VCFX_INFO_SUMMARIZER_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// Function to display help message
void printHelp();

// Structure to hold statistical summaries
struct StatSummary {
    double mean;
    double median;
    double mode;
};

// Function to parse command-line arguments
bool parseArguments(int argc, char *argv[], std::vector<std::string> &info_fields);

// Function to calculate mean
double calculateMean(const std::vector<double> &data);

// Function to calculate median
double calculateMedian(std::vector<double> data);

// Function to calculate mode
double calculateMode(const std::vector<double> &data);

// Function to parse the INFO field and collect specified fields
bool summarizeInfoFields(std::istream &in, std::ostream &out, const std::vector<std::string> &info_fields);

// Memory-mapped file processing (fast path)
bool summarizeInfoFieldsMmap(const char* filepath, std::ostream& out,
                              const std::vector<std::string>& info_fields, bool quiet = false);

#endif // VCFX_INFO_SUMMARIZER_H
