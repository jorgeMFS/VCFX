#include "VCFX_gl_filter.h"
#include <getopt.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <string>
// Implementation of VCFXGLFilter
int VCFXGLFilter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string filterCondition;

    static struct option long_options[] = {
        {"help",    no_argument,       0, 'h'},
        {"filter",  required_argument, 0, 'f'},
        {0,         0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hf:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'f':
                filterCondition = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || filterCondition.empty()) {
        displayHelp();
        return 1;
    }

    // Filter VCF based on genotype likelihood
    filterByGL(std::cin, std::cout, filterCondition);

    return 0;
}

void VCFXGLFilter::displayHelp() {
    std::cout << "VCFX_gl_filter: Filter VCF based on genotype likelihood scores (e.g., GQ > 20).\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_gl_filter --filter \"<CONDITION>\" [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                Display this help message and exit\n";
    std::cout << "  -f, --filter <CONDITION>  Specify the genotype likelihood filter condition (e.g., GQ>20)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_gl_filter --filter \"GQ>20\" < input.vcf > filtered.vcf\n";
}

void VCFXGLFilter::filterByGL(std::istream& in, std::ostream& out, const std::string& filterCondition) {
    // Parse the filter condition using regex (e.g., "GQ>20")
    std::regex conditionRegex(R"((\w+)\s*(>=|<=|>|<|==|!=)\s*(\d+))");
    std::smatch matches;
    if (!std::regex_match(filterCondition, matches, conditionRegex)) {
        std::cerr << "Error: Invalid filter condition format. Expected format like \"GQ>20\".\n";
        return;
    }

    std::string field = matches[1];
    std::string op = matches[2];
    double threshold = std::stod(matches[3]);

    // Process VCF
    bool headerParsed = false;
    std::vector<std::string> headerFields;
    std::vector<int> filterIndices;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            out << "\n";
            continue;
        }

        if (line[0] == '#') {
            if (line.substr(0, 6) == "#CHROM") {
                std::stringstream ss(line);
                std::string fieldStr;
                // Extract header fields
                while (std::getline(ss, fieldStr, '\t')) {
                    headerFields.push_back(fieldStr);
                }

                // Write header as-is
                out << line << "\n";
                headerParsed = true;
            } else {
                // Write other header lines as-is
                out << line << "\n";
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header line with #CHROM not found.\n";
            return;
        }

        // Parse VCF data lines
        std::stringstream ss(line);
        std::string fieldVal;
        std::vector<std::string> fieldsVec;
        while (std::getline(ss, fieldVal, '\t')) {
            fieldsVec.push_back(fieldVal);
        }

        if (fieldsVec.size() < 9) {
            std::cerr << "Warning: Invalid VCF line with fewer than 9 fields: " << line << "\n";
            continue;
        }

        std::string format = fieldsVec[8];
        std::stringstream fmt_ss(format);
        std::string fmt_field;
        std::vector<std::string> formatFields;
        while (std::getline(fmt_ss, fmt_field, ':')) {
            formatFields.push_back(fmt_field);
        }

        // Find the index of the specified field in FORMAT
        int fieldIndex = -1;
        for (size_t i = 0; i < formatFields.size(); ++i) {
            if (formatFields[i] == field) {
                fieldIndex = static_cast<int>(i);
                break;
            }
        }

        if (fieldIndex == -1) {
            std::cerr << "Warning: Field \"" << field << "\" not found in FORMAT column.\n";
            // Optionally, skip this variant or include it
            out << line << "\n";
            continue;
        }

        bool pass = true;
        // Iterate over each sample to check the condition
        for (size_t i = 9; i < fieldsVec.size(); ++i) {
            std::string sample = fieldsVec[i];
            std::stringstream samp_ss(sample);
            std::string samp_field;
            std::vector<std::string> sampleFields;
            while (std::getline(samp_ss, samp_field, ':')) {
                sampleFields.push_back(samp_field);
            }

            if (fieldIndex >= static_cast<int>(sampleFields.size())) {
                std::cerr << "Warning: Field index out of range in sample fields.\n";
                pass = false;
                break;
            }

            std::string valueStr = sampleFields[fieldIndex];
            if (valueStr.empty() || valueStr == ".") {
                pass = false;
                break;
            }

            double value;
            try {
                value = std::stod(valueStr);
            } catch (...) {
                std::cerr << "Warning: Unable to convert value \"" << valueStr << "\" to number.\n";
                pass = false;
                break;
            }

            // Apply the filter condition
            if (op == ">") {
                if (!(value > threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "<") {
                if (!(value < threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == ">=") {
                if (!(value >= threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "<=") {
                if (!(value <= threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "==") {
                if (!(value == threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "!=") {
                if (!(value != threshold)) {
                    pass = false;
                    break;
                }
            }
        }

        if (pass) {
            out << line << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXGLFilter glFilter;
    return glFilter.run(argc, argv);
}