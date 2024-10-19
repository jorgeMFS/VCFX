#include "VCFX_field_extractor.h"
#include <sstream>

std::vector<std::string> parseFields(const std::string& record, const std::vector<std::string>& fields) {
    std::vector<std::string> extracted;
    std::stringstream ss(record);
    std::string token;
    std::vector<std::string> vcf_fields;

    // VCF standard fields: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, etc.
    while (std::getline(ss, token, '\t')) {
        vcf_fields.push_back(token);
    }

    // Assuming the first 8 fields are standard
    std::vector<std::string> standard_fields = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"};

    for (const auto& field : fields) {
        auto it = std::find(standard_fields.begin(), standard_fields.end(), field);
        if (it != standard_fields.end()) {
            size_t index = std::distance(standard_fields.begin(), it);
            if (index < vcf_fields.size()) {
                extracted.push_back(vcf_fields[index]);
            } else {
                extracted.push_back(".");
            }
        }
        // Handle additional INFO fields if needed
        else {
            // Parse INFO field
            auto info_it = std::find(standard_fields.begin(), standard_fields.end(), "INFO");
            if (info_it != standard_fields.end()) {
                size_t info_index = std::distance(standard_fields.begin(), info_it);
                if (info_index < vcf_fields.size()) {
                    std::string info = vcf_fields[info_index];
                    std::stringstream info_ss(info);
                    std::string info_token;
                    while (std::getline(info_ss, info_token, ';')) {
                        size_t eq = info_token.find('=');
                        if (eq != std::string::npos) {
                            std::string key = info_token.substr(0, eq);
                            std::string value = info_token.substr(eq + 1);
                            if (key == field) {
                                extracted.push_back(value);
                                break;
                            }
                        } else {
                            if (info_token == field) {
                                extracted.push_back("1");
                                break;
                            }
                        }
                    }
                } else {
                    extracted.push_back(".");
                }
            } else {
                extracted.push_back(".");
            }
        }
    }

    return extracted;
}

void extractFields(std::istream& in, std::ostream& out, const std::vector<std::string>& fields) {
    std::string line;
    // Print header
    for (size_t i = 0; i < fields.size(); ++i) {
        out << fields[i];
        if (i != fields.size() - 1) {
            out << "\t";
        }
    }
    out << "\n";

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip header lines
        }
        std::vector<std::string> extracted = parseFields(line, fields);
        for (size_t i = 0; i < extracted.size(); ++i) {
            out << extracted[i];
            if (i != extracted.size() - 1) {
                out << "\t";
            }
        }
        out << "\n";
    }
}

int main(int argc, char* argv[]) {
    // Simple command-line argument parsing
    std::vector<std::string> fields;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.find("--fields") == 0) {
            size_t eq = arg.find('=');
            if (eq != std::string::npos) {
                std::string fields_str = arg.substr(eq + 1);
                std::stringstream ss(fields_str);
                std::string field;
                while (std::getline(ss, field, ',')) {
                    fields.push_back(field);
                }
            } else if (i + 1 < argc) {
                std::string fields_str = argv[++i];
                std::stringstream ss(fields_str);
                std::string field;
                while (std::getline(ss, field, ',')) {
                    fields.push_back(field);
                }
            }
        }
    }

    if (fields.empty()) {
        std::cerr << "No fields specified. Usage: " << argv[0] << " --fields \"CHROM,POS,ID,QUAL\"" << std::endl;
        return 1;
    }

    extractFields(std::cin, std::cout, fields);
    return 0;
}
