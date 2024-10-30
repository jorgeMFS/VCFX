#include "VCFX_reformatter.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cctype>

// Implementation

int VCFXReformatter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::vector<std::string> compressFields;
    std::vector<std::string> reorderInfoFields;
    std::vector<std::string> reorderFormatFields;

    static struct option long_options[] = {
        {"help",            no_argument,       0, 'h'},
        {"compress-info",   required_argument, 0, 'c'},
        {"compress-format", required_argument, 0, 'f'},
        {"reorder-info",    required_argument, 0, 'i'},
        {"reorder-format",  required_argument, 0, 'o'},
        {0,                 0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hc:f:i:o:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'c':
                {
                    // Compress INFO fields separated by commas
                    std::string fields(optarg);
                    std::stringstream ss(fields);
                    std::string field;
                    while (std::getline(ss, field, ',')) {
                        // Trim whitespace
                        field.erase(field.find_last_not_of(" \n\r\t")+1);
                        field.erase(0, field.find_first_not_of(" \n\r\t"));
                        if (!field.empty()) {
                            compressFields.push_back(field);
                        }
                    }
                }
                break;
            case 'f':
                {
                    // Compress FORMAT fields separated by commas
                    std::string fields(optarg);
                    std::stringstream ss(fields);
                    std::string field;
                    while (std::getline(ss, field, ',')) {
                        // Trim whitespace
                        field.erase(field.find_last_not_of(" \n\r\t")+1);
                        field.erase(0, field.find_first_not_of(" \n\r\t"));
                        if (!field.empty()) {
                            compressFields.push_back(field);
                        }
                    }
                }
                break;
            case 'i':
                {
                    // Reorder INFO fields, order separated by commas
                    std::string fields(optarg);
                    std::stringstream ss(fields);
                    std::string field;
                    while (std::getline(ss, field, ',')) {
                        // Trim whitespace
                        field.erase(field.find_last_not_of(" \n\r\t")+1);
                        field.erase(0, field.find_first_not_of(" \n\r\t"));
                        if (!field.empty()) {
                            reorderInfoFields.push_back(field);
                        }
                    }
                }
                break;
            case 'o':
                {
                    // Reorder FORMAT fields, order separated by commas
                    std::string fields(optarg);
                    std::stringstream ss(fields);
                    std::string field;
                    while (std::getline(ss, field, ',')) {
                        // Trim whitespace
                        field.erase(field.find_last_not_of(" \n\r\t")+1);
                        field.erase(0, field.find_first_not_of(" \n\r\t"));
                        if (!field.empty()) {
                            reorderFormatFields.push_back(field);
                        }
                    }
                }
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    // Perform VCF reformatting
    reformatVCF(std::cin, std::cout, compressFields, reorderInfoFields, reorderFormatFields);

    return 0;
}

void VCFXReformatter::displayHelp() {
    std::cout << "VCFX_reformatter: Reformat VCF fields (e.g., compressing or reordering INFO/FORMAT fields).\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_reformatter [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                      Display this help message and exit\n";
    std::cout << "  -c, --compress-info <FIELDS>    Compress specified INFO fields (comma-separated)\n";
    std::cout << "  -f, --compress-format <FIELDS>  Compress specified FORMAT fields (comma-separated)\n";
    std::cout << "  -i, --reorder-info <ORDER>      Reorder INFO fields as per specified order (comma-separated)\n";
    std::cout << "  -o, --reorder-format <ORDER>    Reorder FORMAT fields as per specified order (comma-separated)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_reformatter --compress-info AF,DP --reorder-info AF,DP,INFO < input.vcf > reformatted.vcf\n";
}

void VCFXReformatter::reformatVCF(std::istream& in, std::ostream& out,
                                 const std::vector<std::string>& compressFields,
                                 const std::vector<std::string>& reorderInfoFields,
                                 const std::vector<std::string>& reorderFormatFields) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            // Header lines are printed as-is
            out << line << "\n";
            continue;
        }

        // Split the line into fields
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() < 8) {
            std::cerr << "Warning: Skipping invalid VCF line (less than 8 fields): " << line << "\n";
            continue;
        }

        std::string chrom = fields[0];
        std::string pos = fields[1];
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        std::string qual = fields[5];
        std::string filter = fields[6];
        std::string info = fields[7];

        // If there are genotype fields, they start from index 8
        std::vector<std::string> genotypeFields;
        std::string formatField;
        if (fields.size() > 8) {
            formatField = fields[8];
            genotypeFields.assign(fields.begin() + 9, fields.end());
        }

        // Compress specified INFO fields
        if (!compressFields.empty()) {
            info = compressFieldsFunction(info, compressFields);
        }

        // Reorder INFO fields
        if (!reorderInfoFields.empty()) {
            info = reorderInfo(info, reorderInfoFields);
        }

        // Compress specified FORMAT fields
        if (!compressFields.empty() && !formatField.empty()) {
            formatField = compressFieldsFunction(formatField, compressFields);
        }

        // Reorder FORMAT fields
        if (!reorderFormatFields.empty() && !formatField.empty()) {
            formatField = reorderFormat(formatField, reorderFormatFields);
        }

        // Reconstruct the VCF line
        std::ostringstream oss;
        oss << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt
            << "\t" << qual << "\t" << filter << "\t" << info;

        if (!formatField.empty()) {
            oss << "\t" << formatField;
            for (const auto& gt : genotypeFields) {
                oss << "\t" << gt;
            }
        }

        oss << "\n";
        out << oss.str();
    }
}

std::string VCFXReformatter::compressFieldsFunction(const std::string& fieldValue, const std::vector<std::string>& fieldsToCompress) {
    // Example compression: Remove specified key-value pairs from INFO or FORMAT fields
    std::stringstream ss(fieldValue);
    std::string token;
    std::vector<std::string> compressedTokens;

    while (std::getline(ss, token, ';')) {
        bool toCompress = false;
        for (const auto& key : fieldsToCompress) {
            if (token.find(key + "=") == 0) { // Starts with key=
                toCompress = true;
                break;
            }
        }
        if (!toCompress) {
            compressedTokens.push_back(token);
        }
    }

    // Reconstruct the field
    std::string compressedField;
    for (size_t i = 0; i < compressedTokens.size(); ++i) {
        compressedField += compressedTokens[i];
        if (i != compressedTokens.size() - 1) {
            compressedField += ";";
        }
    }

    return compressedField;
}

std::string VCFXReformatter::reorderInfo(const std::string& infoField, const std::vector<std::string>& reorderOrder) {
    // Split INFO field into key-value pairs
    std::stringstream ss(infoField);
    std::string token;
    std::unordered_map<std::string, std::string> infoMap;
    std::vector<std::string> keys;

    while (std::getline(ss, token, ';')) {
        size_t eqPos = token.find('=');
        if (eqPos != std::string::npos) {
            std::string key = token.substr(0, eqPos);
            std::string value = token.substr(eqPos + 1);
            infoMap[key] = value;
            keys.push_back(key);
        } else {
            // Flag without value
            infoMap[token] = "";
            keys.push_back(token);
        }
    }

    // Reorder based on reorderOrder
    std::ostringstream oss;
    for (size_t i = 0; i < reorderOrder.size(); ++i) {
        const std::string& key = reorderOrder[i];
        if (infoMap.find(key) != infoMap.end()) {
            if (!infoMap[key].empty()) {
                oss << key << "=" << infoMap[key];
            } else {
                oss << key;
            }
            if (i != reorderOrder.size() - 1) {
                oss << ";";
            }
            infoMap.erase(key);
        }
    }

    // Append remaining keys in original order
    for (const auto& key : keys) {
        if (infoMap.find(key) != infoMap.end()) {
            if (!oss.str().empty() && oss.str().back() != ';') {
                oss << ";";
            }
            if (!infoMap[key].empty()) {
                oss << key << "=" << infoMap[key];
            } else {
                oss << key;
            }
            infoMap.erase(key);
        }
    }

    return oss.str();
}

std::string VCFXReformatter::reorderFormat(const std::string& formatField, const std::vector<std::string>& reorderOrder) {
    // Split FORMAT field into keys
    std::stringstream ss(formatField);
    std::string token;
    std::vector<std::string> formatKeys;
    std::vector<std::string> remainingKeys;

    while (std::getline(ss, token, ':')) {
        formatKeys.push_back(token);
    }

    // Reorder based on reorderOrder
    std::vector<std::string> reorderedFormat;
    for (const auto& key : reorderOrder) {
        auto it = std::find(formatKeys.begin(), formatKeys.end(), key);
        if (it != formatKeys.end()) {
            reorderedFormat.push_back(key);
            formatKeys.erase(it);
        }
    }

    // Append remaining keys in original order
    for (const auto& key : formatKeys) {
        reorderedFormat.push_back(key);
    }

    // Reconstruct the FORMAT field
    std::string reorderedFormatStr;
    for (size_t i = 0; i < reorderedFormat.size(); ++i) {
        reorderedFormatStr += reorderedFormat[i];
        if (i != reorderedFormat.size() - 1) {
            reorderedFormatStr += ":";
        }
    }

    return reorderedFormatStr;
}

int main(int argc, char* argv[]) {
    VCFXReformatter reformatter;
    return reformatter.run(argc, argv);
}