#include "VCFX_info_summarizer.h"
#include <sstream>
#include <algorithm>
#include <cctype>
#include <unordered_set>
#include <iomanip>
#include <map>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_info_summarizer\n"
              << "Usage: VCFX_info_summarizer [OPTIONS]\n\n"
              << "Options:\n"
              << "  --info, -i \"FIELD1,FIELD2\"   Specify the INFO fields to summarize (e.g., \"DP,AF\").\n"
              << "  --help, -h                    Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Summarizes numeric fields in the INFO column of a VCF file by calculating statistics such as mean, median, and mode.\n\n"
              << "Examples:\n"
              << "  ./VCFX_info_summarizer --info \"DP,AF\" < input.vcf > summary_stats.tsv\n";
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
    std::cerr << "Error: INFO fields not specified.\n";
    std::cerr << "Use --help for usage information.\n";
    return false;
}

// Function to calculate mean
double calculateMean(const std::vector<double>& data) {
    if (data.empty()) return 0.0;
    double sum = 0.0;
    for (const auto& val : data) {
        sum += val;
    }
    return sum / static_cast<double>(data.size());
}

// Function to calculate median
double calculateMedian(std::vector<double> data) {
    if (data.empty()) return 0.0;
    std::sort(data.begin(), data.end());
    size_t n = data.size();
    if (n % 2 == 0) {
        return (data[n / 2 - 1] + data[n / 2]) / 2.0;
    } else {
        return data[n / 2];
    }
}

// Function to calculate mode
double calculateMode(const std::vector<double>& data) {
    if (data.empty()) return 0.0;
    std::unordered_map<double, int> frequency;
    int max_freq = 0;
    double mode = data[0];
    for (const auto& val : data) {
        frequency[val]++;
        if (frequency[val] > max_freq) {
            max_freq = frequency[val];
            mode = val;
        }
    }
    return mode;
}

// Function to parse the INFO field and collect specified fields
bool summarizeInfoFields(std::istream& in, std::ostream& out, const std::vector<std::string>& info_fields) {
    std::string line;
    bool header_found = false;

    // Map to store vectors of values for each INFO field
    std::map<std::string, std::vector<double>> info_data;

    // Initialize map keys
    for (const auto& field : info_fields) {
        info_data[field] = std::vector<double>();
    }

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                header_found = true;
            }
            continue; // Skip header lines
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info;
        if (!std::getline(ss, chrom, '\t') ||
            !std::getline(ss, pos, '\t') ||
            !std::getline(ss, id, '\t') ||
            !std::getline(ss, ref, '\t') ||
            !std::getline(ss, alt, '\t') ||
            !std::getline(ss, qual, '\t') ||
            !std::getline(ss, filter, '\t') ||
            !std::getline(ss, info, '\t')) {
            std::cerr << "Warning: Skipping malformed VCF line: " << line << "\n";
            continue;
        }

        // Parse INFO field into key-value pairs
        std::stringstream info_ss(info);
        std::string kv;
        std::unordered_map<std::string, std::string> info_map;
        while (std::getline(info_ss, kv, ';')) {
            size_t eq = kv.find('=');
            if (eq != std::string::npos) {
                std::string key = kv.substr(0, eq);
                std::string value = kv.substr(eq + 1);
                info_map[key] = value;
            } else {
                // Flags without values are set to "1"
                info_map[kv] = "1";
            }
        }

        // Extract and store specified INFO fields
        for (const auto& field : info_fields) {
            if (info_map.find(field) != info_map.end()) {
                std::string value_str = info_map[field];
                // Handle multiple values separated by commas
                std::stringstream val_ss(value_str);
                std::string val;
                while (std::getline(val_ss, val, ',')) {
                    try {
                        double val_num = std::stod(val);
                        info_data[field].push_back(val_num);
                    } catch (const std::invalid_argument& e) {
                        // Non-numeric value, skip
                        std::cerr << "Warning: Non-numeric value for field " << field << " in line: " << line << "\n";
                        continue;
                    }
                }
            }
        }
    }

    // Output summary statistics
    out << "INFO_Field\tMean\tMedian\tMode\n";
    for (const auto& field : info_fields) {
        const auto& data = info_data[field];
        if (data.empty()) {
            out << field << "\tNA\tNA\tNA\n";
            continue;
        }
        double mean = calculateMean(data);
        double median = calculateMedian(data);
        double mode = calculateMode(data);
        out << field << "\t" << std::fixed << std::setprecision(4) << mean << "\t" << median << "\t" << mode << "\n";
    }

    return true;
}

int main(int argc, char* argv[]) {
    std::vector<std::string> info_fields;

    // Parse command-line arguments
    if (!parseArguments(argc, argv, info_fields)) {
        return 1;
    }

    // Summarize INFO fields
    bool success = summarizeInfoFields(std::cin, std::cout, info_fields);
    return success ? 0 : 1;
}
