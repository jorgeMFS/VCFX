#include "VCFX_record_filter.h"
#include <sstream>
#include <vector>
#include <algorithm>

// Helper function to trim whitespace from both ends of a string
static inline std::string trim(const std::string& s) {
    auto start = s.begin();
    while (start != s.end() && isspace(*start)) {
        start++;
    }

    auto end = s.end();
    do {
        end--;
    } while (distance(start, end) > 0 && isspace(*end));

    return std::string(start, end + 1);
}

// Function to parse individual filter criterion
bool parseSingleCriterion(const std::string& token, FilterCriterion& criterion) {
    size_t pos = token.find_first_of("><=!");
    if (pos == std::string::npos) {
        return false;
    }

    // Determine the operator
    Operator op;
    std::string op_str;
    if (token[pos] == '>') {
        if (pos + 1 < token.size() && token[pos + 1] == '=') {
            op = Operator::GREATER_EQUAL;
            op_str = ">=";
            pos += 1;
        } else {
            op = Operator::GREATER_THAN;
            op_str = ">";
        }
    } else if (token[pos] == '<') {
        if (pos + 1 < token.size() && token[pos + 1] == '=') {
            op = Operator::LESS_EQUAL;
            op_str = "<=";
            pos += 1;
        } else {
            op = Operator::LESS_THAN;
            op_str = "<";
        }
    } else if (token[pos] == '=') {
        if (pos + 1 < token.size() && token[pos + 1] == '=') {
            op = Operator::EQUAL;
            op_str = "==";
            pos += 1;
        } else {
            // Single '=' treated as '=='
            op = Operator::EQUAL;
            op_str = "==";
        }
    } else {
        // Unsupported operator
        return false;
    }

    // Extract field and value
    std::string field = trim(token.substr(0, pos));
    std::string value_str = trim(token.substr(pos + 1));

    if (field.empty() || value_str.empty()) {
        return false;
    }

    try {
        double value = std::stod(value_str);
        criterion.field = field;
        criterion.op = op;
        criterion.value = value;
    } catch (const std::invalid_argument&) {
        return false;
    }

    return true;
}

// Function to parse filter criteria string into structured criteria
bool parseCriteria(const std::string& criteria_str, std::vector<FilterCriterion>& criteria) {
    std::stringstream ss(criteria_str);
    std::string token;

    while (std::getline(ss, token, ';')) {
        token = trim(token);
        if (token.empty()) {
            continue;
        }

        FilterCriterion criterion;
        if (!parseSingleCriterion(token, criterion)) {
            std::cerr << "Failed to parse filter criterion: '" << token << "'" << std::endl;
            return false;
        }

        criteria.push_back(criterion);
    }

    if (criteria.empty()) {
        std::cerr << "No valid filter criteria found." << std::endl;
        return false;
    }

    return true;
}

// Function to split a string by a delimiter
static inline std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::stringstream ss(s);
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Function to retrieve the value of a specific field from a VCF record
bool getFieldValue(const std::vector<std::string>& fields, const std::string& field_name, double& value) {
    static const std::vector<std::string> standard_fields = {
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
    };

    auto it = std::find(standard_fields.begin(), standard_fields.end(), field_name);
    if (it != standard_fields.end()) {
        size_t index = std::distance(standard_fields.begin(), it);
        if (index >= fields.size()) {
            return false;
        }

        if (field_name == "POS") {
            try {
                value = std::stod(fields[index]);
            } catch (...) {
                return false;
            }
            return true;
        } else if (field_name == "QUAL") {
            if (fields[index] == ".") {
                value = 0.0;
            } else {
                try {
                    value = std::stod(fields[index]);
                } catch (...) {
                    return false;
                }
            }
            return true;
        } else {
            // Currently, only POS and QUAL are supported as numeric fields
            return false;
        }
    } else {
        // Attempt to extract from INFO field
        if (standard_fields.size() < 8) {
            return false;
        }
        std::string info_field = fields[7];
        std::vector<std::string> info_entries = split(info_field, ';');
        for (const auto& entry : info_entries) {
            size_t eq = entry.find('=');
            if (eq != std::string::npos) {
                std::string key = trim(entry.substr(0, eq));
                std::string val_str = trim(entry.substr(eq + 1));
                if (key == field_name) {
                    try {
                        value = std::stod(val_str);
                        return true;
                    } catch (...) {
                        return false;
                    }
                }
            } else {
                if (entry == field_name) {
                    // Flag present, treat as true (1.0)
                    value = 1.0;
                    return true;
                }
            }
        }
    }

    return false;
}

// Function to apply all filters to a single record
bool applyFilters(const std::string& record, const std::vector<FilterCriterion>& criteria) {
    std::vector<std::string> fields = split(record, '\t');

    for (const auto& criterion : criteria) {
        double field_value = 0.0;
        if (!getFieldValue(fields, criterion.field, field_value)) {
            // If the field is not found or not numeric, skip this criterion
            return false;
        }

        bool condition_met = false;
        switch (criterion.op) {
            case Operator::GREATER_THAN:
                condition_met = (field_value > criterion.value);
                break;
            case Operator::LESS_THAN:
                condition_met = (field_value < criterion.value);
                break;
            case Operator::GREATER_EQUAL:
                condition_met = (field_value >= criterion.value);
                break;
            case Operator::LESS_EQUAL:
                condition_met = (field_value <= criterion.value);
                break;
            case Operator::EQUAL:
                condition_met = (field_value == criterion.value);
                break;
        }

        if (!condition_met) {
            return false; // If any criterion is not met, reject the record
        }
    }

    return true; // All criteria met
}

// Function to process and filter records
void processRecords(std::istream& in, std::ostream& out, const std::vector<FilterCriterion>& criteria) {
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            continue; // Skip header lines
        }

        if (applyFilters(line, criteria)) {
            out << line << "\n";
        }
    }
}

void printHelp() {
    std::cout << "VCFX_record_filter\n"
              << "Usage: VCFX_record_filter --filter \"CRITERIA\" [OPTIONS]\n\n"
              << "Options:\n"
              << "  --filter, -f          Specify filter criteria (e.g., \"QUAL>30;DP<100\").\n"
              << "  --help, -h            Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Filters VCF records based on specified criteria.\n\n"
              << "Example:\n"
              << "  ./VCFX_record_filter --filter \"QUAL>30;DP<100\" < input.vcf > filtered.vcf\n";
}

int main(int argc, char* argv[]) {
    // Argument parsing
    std::string criteria_str;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
        if (arg == "--filter" || arg == "-f") {
            if (i + 1 < argc) {
                criteria_str = argv[++i];
            } else {
                std::cerr << "Error: --filter option requires an argument.\n";
                return 1;
            }
        } else if (arg.find("--filter=") == 0) {
            criteria_str = arg.substr(9);
        }
    }

    if (criteria_str.empty()) {
        std::cerr << "No filter criteria provided.\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    std::vector<FilterCriterion> criteria;
    if (!parseCriteria(criteria_str, criteria)) {
        std::cerr << "Failed to parse filter criteria.\n";
        return 1;
    }

    processRecords(std::cin, std::cout, criteria);
    return 0;
}
