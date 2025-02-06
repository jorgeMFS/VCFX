./VCFX_format_converter/VCFX_format_converter.cpp
#include "VCFX_format_converter.h"
#include <sstream>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_format_converter\n"
              << "Usage: VCFX_format_converter [OPTIONS]\n\n"
              << "Options:\n"
              << "  --to-bed             Convert VCF to BED format.\n"
              << "  --to-csv             Convert VCF to CSV format.\n"
              << "  --help, -h           Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Converts VCF files to specified formats (BED or CSV).\n\n"
              << "Example:\n"
              << "  ./VCFX_format_converter --to-bed < input.vcf > output.bed\n"
              << "  ./VCFX_format_converter --to-csv < input.vcf > output.csv\n";
}

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], OutputFormat& format) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--to-bed") {
            format = OutputFormat::BED;
            return true;
        } else if (arg == "--to-csv") {
            format = OutputFormat::CSV;
            return true;
        }
    }
    format = OutputFormat::UNKNOWN;
    return false;
}

// Function to convert VCF to BED
void convertVCFtoBED(std::istream& in, std::ostream& out) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip header lines
        }

        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string field;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() < 5) {
            continue; // Insufficient fields
        }

        std::string chrom = fields[0];
        int pos = std::stoi(fields[1]);
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        int end = pos + ref.length(); // Simple end position

        out << chrom << "\t" << pos - 1 << "\t" << end << "\t" << id << "\n";
    }
}

// Function to convert VCF to CSV
void convertVCFtoCSV(std::istream& in, std::ostream& out) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip header lines or handle accordingly
        }

        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string field;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        // Join fields with commas
        for (size_t i = 0; i < fields.size(); ++i) {
            // Handle fields with commas by enclosing in quotes
            if (fields[i].find(',') != std::string::npos) {
                out << "\"" << fields[i] << "\"";
            } else {
                out << fields[i];
            }

            if (i != fields.size() - 1) {
                out << ",";
            }
        }
        out << "\n";
    }
}

int main(int argc, char* argv[]) {
    // Argument parsing
    OutputFormat format;
    if (!parseArguments(argc, argv, format)) {
        std::cerr << "No valid output format specified.\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    // Handle help and invalid formats
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    switch (format) {
        case OutputFormat::BED:
            convertVCFtoBED(std::cin, std::cout);
            break;
        case OutputFormat::CSV:
            convertVCFtoCSV(std::cin, std::cout);
            break;
        default:
            std::cerr << "Unsupported output format.\n";
            std::cerr << "Use --help for usage information.\n";
            return 1;
    }

    return 0;
}


./VCFX_format_converter/VCFX_format_converter.h
#ifndef VCFX_FORMAT_CONVERTER_H
#define VCFX_FORMAT_CONVERTER_H

#include <iostream>
#include <string>
#include <vector>

// Enumeration for output formats
enum class OutputFormat {
    BED,
    CSV,
    UNKNOWN
};

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], OutputFormat& format);

// Function to convert VCF to BED
void convertVCFtoBED(std::istream& in, std::ostream& out);

// Function to convert VCF to CSV
void convertVCFtoCSV(std::istream& in, std::ostream& out);

// Function to display help message
void printHelp();

#endif // VCFX_FORMAT_CONVERTER_H


./VCFX_genotype_query/VCFX_genotype_query.cpp
#include "VCFX_genotype_query.h"
#include <sstream>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_genotype_query\n"
              << "Usage: VCFX_genotype_query [OPTIONS]\n\n"
              << "Options:\n"
              << "  --genotype-query, -g \"GENOTYPE\"  Specify the genotype to query (e.g., \"0/1\", \"1/1\").\n"
              << "  --help, -h                        Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Filters VCF records based on the specified genotype query. Only records matching the genotype "
              << "criteria will be outputted.\n\n"
              << "Example:\n"
              << "  ./VCFX_genotype_query --genotype-query \"0/1\" < input.vcf > output.vcf\n";
}

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], std::string& genotype_query) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--genotype-query" || arg == "-g") && i + 1 < argc) {
            genotype_query = argv[++i];
            return true;
        } else if (arg.find("--genotype-query=") == 0) {
            genotype_query = arg.substr(17);
            return true;
        } else if (arg == "--help" || arg == "-h") {
            printHelp();
            exit(0);
        }
    }
    return false;
}

// Function to split a string by a delimiter and return a vector of tokens
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::stringstream ss(s);
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Function to perform genotype query on VCF records
void genotypeQuery(std::istream& in, std::ostream& out, const std::string& genotype_query) {
    std::string line;
    std::string header;
    std::vector<std::string> sample_names;
    bool header_found = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                header = line;
                std::vector<std::string> fields = split(line, '\t');
                // Samples start from the 10th column (index 9)
                if (fields.size() > 9) {
                    sample_names.assign(fields.begin() + 9, fields.end());
                }
                header_found = true;
            }
            out << line << "\n"; // Print header lines
            continue;
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return;
        }

        std::vector<std::string> fields = split(line, '\t');
        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping invalid VCF line with fewer than 9 fields.\n";
            continue;
        }

        std::string format = fields[8];
        std::vector<std::string> format_fields = split(format, ':');
        int gt_index = -1;
        for (size_t i = 0; i < format_fields.size(); ++i) {
            if (format_fields[i] == "GT") {
                gt_index = static_cast<int>(i);
                break;
            }
        }

        if (gt_index == -1) {
            std::cerr << "Warning: GT field not found in FORMAT column. Skipping line.\n";
            continue;
        }

        bool match_found = false;
        for (size_t i = 9; i < fields.size(); ++i) {
            std::vector<std::string> sample_data = split(fields[i], ':');
            if (static_cast<size_t>(gt_index) >= sample_data.size()) {
                continue; // No GT data for this sample
            }
            if (sample_data[gt_index] == genotype_query) {
                match_found = true;
                break;
            }
        }

        if (match_found) {
            out << line << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    std::string genotype_query;
    if (!parseArguments(argc, argv, genotype_query)) {
        std::cerr << "Usage: " << argv[0] << " --genotype-query \"0/1\" < input.vcf > output.vcf\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    genotypeQuery(std::cin, std::cout, genotype_query);
    return 0;
}


./VCFX_genotype_query/VCFX_genotype_query.h
#ifndef VCFX_GENOTYPE_QUERY_H
#define VCFX_GENOTYPE_QUERY_H

#include <iostream>
#include <string>
#include <vector>

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], std::string& genotype_query);

// Function to display help message
void printHelp();

// Function to perform genotype query on VCF records
void genotypeQuery(std::istream& in, std::ostream& out, const std::string& genotype_query);

#endif // VCFX_GENOTYPE_QUERY_H


./VCFX_gl_filter/VCFX_gl_filter.cpp
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

./VCFX_gl_filter/VCFX_gl_filter.h
#ifndef VCFX_GL_FILTER_H
#define VCFX_GL_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXGLFilter: Header file for Genotype Likelihood Filter tool
class VCFXGLFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input based on genotype likelihood scores (e.g., GQ > threshold)
    void filterByGL(std::istream& in, std::ostream& out, const std::string& filterCondition);
};

#endif // VCFX_GL_FILTER_H


./VCFX_haplotype_extractor/VCFX_haplotype_extractor.cpp
#include "VCFX_haplotype_extractor.h"
#include <sstream>
#include <algorithm>
#include <cctype>
#include <stdexcept>

// Constructor
HaplotypeExtractor::HaplotypeExtractor() {}

// Function to display help message
void printHelp() {
    std::cout << "VCFX_haplotype_extractor\n"
              << "Usage: VCFX_haplotype_extractor [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h                     Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Extracts phased haplotype blocks from genotype data in a VCF file. "
              << "It reconstructs haplotypes for each sample by analyzing phased genotype information, "
              << "grouping contiguous variants into haplotype blocks based on consistent phasing and proximity.\n\n"
              << "Examples:\n"
              << "  ./VCFX_haplotype_extractor < phased.vcf > haplotypes.tsv\n";
}

// Utility function to split a string by a delimiter
std::vector<std::string> HaplotypeExtractor::splitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string temp;
    while (std::getline(ss, temp, delimiter)) {
        tokens.push_back(temp);
    }
    return tokens;
}

// Parses the VCF header to extract sample names
bool HaplotypeExtractor::parseHeader(const std::string& headerLine) {
    std::vector<std::string> fields = splitString(headerLine, '\t');
    if (fields.size() <= 8) {
        std::cerr << "Error: VCF header does not contain sample columns.\n";
        return false;
    }

    // Extract sample names starting from the 9th column
    sampleNames.assign(fields.begin() + 9, fields.end());
    numSamples = sampleNames.size();
    return true;
}

// Validates if all samples are phased
bool HaplotypeExtractor::areAllSamplesPhased(const std::vector<std::string>& genotypeFields) {
    for (const auto& gt : genotypeFields) {
        // A phased genotype contains '|', e.g., 1
        if (gt.find('|') == std::string::npos) {
            return false;
        }
    }
    return true;
}

// Processes a single VCF variant and updates haplotype blocks
bool HaplotypeExtractor::processVariant(const std::vector<std::string>& fields, std::vector<HaplotypeBlock>& haplotypeBlocks, int currentPos) {
    if (fields.size() < 9) {
        std::cerr << "Warning: Skipping invalid VCF line with fewer than 9 fields.\n";
        return false;
    }

    std::string chrom = fields[0];
    int pos = 0;
    try {
        pos = std::stoi(fields[1]);
    } catch (...) {
        std::cerr << "Warning: Invalid POS value. Skipping line.\n";
        return false;
    }

    std::string format = fields[8];
    std::vector<std::string> formatFields = splitString(format, ':');
    int gtIndex = -1;
    for (size_t i = 0; i < formatFields.size(); ++i) {
        if (formatFields[i] == "GT") {
            gtIndex = static_cast<int>(i);
            break;
        }
    }

    if (gtIndex == -1) {
        std::cerr << "Warning: GT field not found in FORMAT column. Skipping line.\n";
        return false;
    }

    std::vector<std::string> genotypeFields;
    bool allPhased = true;
    for (size_t i = 9; i < fields.size(); ++i) {
        std::vector<std::string> sampleData = splitString(fields[i], ':');
        if (static_cast<size_t>(gtIndex) >= sampleData.size()) {
            genotypeFields.push_back("./.");
            allPhased = false;
            continue; // No GT data for this sample
        }
        std::string gt = sampleData[gtIndex];
        genotypeFields.push_back(gt);
        if (gt.find('|') == std::string::npos) {
            allPhased = false;
        }
    }

    if (!allPhased) {
        std::cerr << "Warning: Not all samples are phased at position " << pos << " on " << chrom << ". Skipping haplotype extraction for this variant.\n";
        return false;
    }

    // Update haplotype blocks
    updateHaplotypeBlocks(genotypeFields, haplotypeBlocks, pos, chrom);

    return true;
}

// Updates haplotype blocks with new variant information
void HaplotypeExtractor::updateHaplotypeBlocks(const std::vector<std::string>& genotypeFields, std::vector<HaplotypeBlock>& haplotypeBlocks, int pos, const std::string& chrom) {
    if (haplotypeBlocks.empty()) {
        // Initialize the first haplotype block
        HaplotypeBlock block;
        block.chrom = chrom;
        block.start = pos;
        block.end = pos;
        for (const auto& gt : genotypeFields) {
            block.haplotypes.push_back(gt);
        }
        haplotypeBlocks.push_back(block);
        return;
    }

    // Get the last haplotype block
    HaplotypeBlock& lastBlock = haplotypeBlocks.back();

    // Determine if the current variant should be part of the last block
    // Criteria:
    // 1. Same chromosome
    // 2. Distance between variants is below a threshold (e.g., 100 kb)
    // 3. Phasing is consistent (already ensured)

    const int distanceThreshold = 100000; // 100 kb
    if (chrom == lastBlock.chrom && (pos - lastBlock.end) <= distanceThreshold) {
        // Extend the last haplotype block
        lastBlock.end = pos;
        for (size_t i = 0; i < genotypeFields.size(); ++i) {
            lastBlock.haplotypes[i] += "|" + genotypeFields[i];
        }
    } else {
        // Start a new haplotype block
        HaplotypeBlock newBlock;
        newBlock.chrom = chrom;
        newBlock.start = pos;
        newBlock.end = pos;
        for (const auto& gt : genotypeFields) {
            newBlock.haplotypes.push_back(gt);
        }
        haplotypeBlocks.push_back(newBlock);
    }
}

// Parses the VCF file and extracts haplotype blocks
bool HaplotypeExtractor::extractHaplotypes(std::istream& in, std::ostream& out) {
    std::string line;
    bool header_found = false;
    std::vector<HaplotypeBlock> haplotypeBlocks;

    // Read the VCF file line by line
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                // Parse header to extract sample names
                if (!parseHeader(line)) {
                    return false;
                }
                header_found = true;
            }
            continue; // Skip header lines
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        // Split the VCF record into fields
        std::vector<std::string> fields = splitString(line, '\t');

        // Process the variant to extract haplotypes
        if (!processVariant(fields, haplotypeBlocks, 0)) {
            continue; // Skip variants that couldn't be processed
        }
    }

    // Output the haplotype blocks
    out << "CHROM\tSTART\tEND\t" ;
    for (const auto& sample : sampleNames) {
        out << sample << "\t";
    }
    out << "\n";

    for (const auto& block : haplotypeBlocks) {
        out << block.chrom << "\t" << block.start << "\t" << block.end << "\t";
        for (size_t i = 0; i < block.haplotypes.size(); ++i) {
            out << block.haplotypes[i];
            if (i != block.haplotypes.size() - 1) {
                out << "\t";
            }
        }
        out << "\n";
    }

    return true;
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

    // Instantiate the haplotype extractor
    HaplotypeExtractor extractor;

    // Extract haplotype blocks from stdin and output to stdout
    bool success = extractor.extractHaplotypes(std::cin, std::cout);
    return success ? 0 : 1;
}


./VCFX_haplotype_extractor/VCFX_haplotype_extractor.h
#ifndef VCFX_HAPLOTYPE_EXTRACTOR_H
#define VCFX_HAPLOTYPE_EXTRACTOR_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// Function to display help message
void printHelp();

// Structure to represent a haplotype block
struct HaplotypeBlock {
    std::string chrom;
    int start;
    int end;
    std::vector<std::string> haplotypes; // One haplotype per sample
};

// Class to handle haplotype extraction
class HaplotypeExtractor {
public:
    HaplotypeExtractor();
    ~HaplotypeExtractor() = default;

    // Parses the VCF file and extracts haplotype blocks
    bool extractHaplotypes(std::istream& in, std::ostream& out);

private:
    std::vector<std::string> sampleNames;
    size_t numSamples = 0;

    // Parses the VCF header to extract sample names
    bool parseHeader(const std::string& headerLine);

    // Splits a string by a delimiter
    std::vector<std::string> splitString(const std::string& str, char delimiter);

    // Reconstructs haplotype blocks based on phased genotypes
    bool processVariant(const std::vector<std::string>& fields, std::vector<HaplotypeBlock>& haplotypeBlocks, int currentPos);

    // Validates if all samples are phased
    bool areAllSamplesPhased(const std::vector<std::string>& genotypeFields);

    // Updates haplotype blocks with new variant information
    void updateHaplotypeBlocks(const std::vector<std::string>& genotypeFields, std::vector<HaplotypeBlock>& haplotypeBlocks, int pos, const std::string& chrom);
};

#endif // VCFX_HAPLOTYPE_EXTRACTOR_H


./VCFX_haplotype_phaser/VCFX_haplotype_phaser.cpp
#include "VCFX_haplotype_phaser.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <numeric>

int VCFXHaplotypePhaser::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    double ldThreshold = 0.8; // Default LD threshold

    static struct option long_options[] = {
        {"help",          no_argument,       0, 'h'},
        {"ld-threshold", required_argument, 0, 'l'},
        {0,               0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hl:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'l':
                try {
                    ldThreshold = std::stod(optarg);
                } catch (const std::invalid_argument&) {
                    std::cerr << "Error: Invalid LD threshold value.\n";
                    displayHelp();
                    return 1;
                }
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || ldThreshold < 0.0 || ldThreshold > 1.0) {
        displayHelp();
        return 1;
    }

    // Perform haplotype phasing
    phaseHaplotypes(std::cin, std::cout, ldThreshold);

    return 0;
}

void VCFXHaplotypePhaser::displayHelp() {
    std::cout << "VCFX_haplotype_phaser: Group variants into phased haplotype blocks based on linkage disequilibrium.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_haplotype_phaser [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                    Display this help message and exit\n";
    std::cout << "  -l, --ld-threshold <VALUE>    Specify the LD threshold for grouping (0.0 - 1.0, default: 0.8)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_haplotype_phaser --ld-threshold 0.9 < input.vcf > phased_blocks.txt\n";
}

void VCFXHaplotypePhaser::phaseHaplotypes(std::istream& in, std::ostream& out, double ldThreshold) {
    std::string line;
    std::vector<std::string> sampleNames;
    bool headerParsed = false;
    std::vector<std::vector<int>> genotypeMatrix; // Rows: Variants, Columns: Samples

    // Read header to get sample names
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line.substr(0, 6) == "#CHROM") {
            std::stringstream ss(line);
            std::string field;
            std::vector<std::string> headers;
            while (std::getline(ss, field, '\t')) {
                headers.push_back(field);
            }
            // Sample names start from the 10th column
            for (size_t i = 9; i < headers.size(); ++i) {
                sampleNames.push_back(headers[i]);
            }
            headerParsed = true;
            // Output header
            out << line << "\n";
            break;
        }
    }

    if (!headerParsed) {
        std::cerr << "Error: VCF header line with #CHROM not found.\n";
        return;
    }

    // Read and store all variants
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info, format;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        std::vector<int> genotypes;
        std::string genotype;
        while (ss >> genotype) {
            // Parse genotype (assuming GT field first)
            size_t sep = genotype.find_first_of("/|");
            if (sep == std::string::npos) {
                genotypes.push_back(-1); // Missing or unknown genotype
                continue;
            }
            std::string gt1 = genotype.substr(0, sep);
            std::string gt2 = genotype.substr(sep + 1);
            int allele1 = (gt1.empty()) ? -1 : std::stoi(gt1);
            int allele2 = (gt2.empty()) ? -1 : std::stoi(gt2);
            // For simplicity, store count of reference alleles
            if (allele1 == -1 || allele2 == -1) {
                genotypes.push_back(-1); // Missing genotype
            } else {
                genotypes.push_back(allele1 + allele2);
            }
        }

        genotypeMatrix.push_back(genotypes);
    }

    if (genotypeMatrix.empty()) {
        std::cerr << "Error: No variant data found in VCF.\n";
        return;
    }

    // Group variants into haplotype blocks based on LD
    std::vector<std::vector<int>> haplotypeBlocks = groupVariants(genotypeMatrix, ldThreshold);

    // Output haplotype blocks
    out << "#HAPLOTYPE_BLOCKS_START\n";
    for (size_t i = 0; i < haplotypeBlocks.size(); ++i) {
        out << "Block " << (i + 1) << ": ";
        for (size_t j = 0; j < haplotypeBlocks[i].size(); ++j) {
            out << haplotypeBlocks[i][j];
            if (j != haplotypeBlocks[i].size() - 1) {
                out << ", ";
            }
        }
        out << "\n";
    }
    out << "#HAPLOTYPE_BLOCKS_END\n";
}

std::vector<std::vector<int>> VCFXHaplotypePhaser::groupVariants(const std::vector<std::vector<int>>& genotypeMatrix, double ldThreshold) {
    std::vector<std::vector<int>> blocks;
    std::vector<int> currentBlock;

    for (size_t i = 0; i < genotypeMatrix.size(); ++i) {
        if (currentBlock.empty()) {
            currentBlock.push_back(i);
        } else {
            // Calculate LD between the last variant in the block and the current variant
            int lastVariantIndex = currentBlock.back();
            double ld = calculateLD(genotypeMatrix[lastVariantIndex], genotypeMatrix[i]);

            if (ld >= ldThreshold) {
                currentBlock.push_back(i);
            } else {
                blocks.push_back(currentBlock);
                currentBlock.clear();
                currentBlock.push_back(i);
            }
        }
    }

    if (!currentBlock.empty()) {
        blocks.push_back(currentBlock);
    }

    return blocks;
}

double VCFXHaplotypePhaser::calculateLD(const std::vector<int>& variant1, const std::vector<int>& variant2) {
    // Calculate r² linkage disequilibrium between two variants
    int n = 0;
    int sum_X = 0, sum_Y = 0;
    int sum_XY = 0;
    int sum_X2 = 0, sum_Y2 = 0;

    for (size_t i = 0; i < variant1.size(); ++i) {
        int x = variant1[i];
        int y = variant2[i];
        if (x == -1 || y == -1) {
            continue; // Skip missing genotypes
        }

        n++;
        sum_X += x;
        sum_Y += y;
        sum_XY += x * y;
        sum_X2 += x * x;
        sum_Y2 += y * y;
    }

    if (n == 0) {
        return 0.0; // No data to calculate LD
    }

    double mean_X = static_cast<double>(sum_X) / n;
    double mean_Y = static_cast<double>(sum_Y) / n;

    double cov = (static_cast<double>(sum_XY) / n) - (mean_X * mean_Y);
    double var_X = (static_cast<double>(sum_X2) / n) - (mean_X * mean_X);
    double var_Y = (static_cast<double>(sum_Y2) / n) - (mean_Y * mean_Y);

    if (var_X == 0.0 || var_Y == 0.0) {
        return 0.0; // Avoid division by zero
    }

    double r = cov / (std::sqrt(var_X) * std::sqrt(var_Y));
    double r_squared = r * r;

    return r_squared;
}

int main(int argc, char* argv[]) {
    VCFXHaplotypePhaser haplotypePhaser;
    return haplotypePhaser.run(argc, argv);
}

./VCFX_haplotype_phaser/VCFX_haplotype_phaser.h
#ifndef VCFX_HAPLOTYPE_PHASER_H
#define VCFX_HAPLOTYPE_PHASER_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// VCFXHaplotypePhaser: Header file for Haplotype-based Phasing Tool
class VCFXHaplotypePhaser {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Phases haplotypes in the VCF input
    void phaseHaplotypes(std::istream& in, std::ostream& out, double ldThreshold);

    // Groups variants into haplotype blocks based on linkage disequilibrium
    std::vector<std::vector<int>> groupVariants(const std::vector<std::vector<int>>& genotypeMatrix, double ldThreshold);

    // Calculates linkage disequilibrium (r²) between two variants
    double calculateLD(const std::vector<int>& variant1, const std::vector<int>& variant2);
};

#endif // VCFX_HAPLOTYPE_PHASER_H

./VCFX_header_parser/VCFX_header_parser.cpp
#include "VCFX_header_parser.h"
#include <iostream>
#include <sstream>

void printHelp() {
    std::cout << "VCFX_header_parser\n"
              << "Usage: VCFX_header_parser [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h            Display this help message and exit.\n"
              << "\n"
              << "Description:\n"
              << "  Extracts and displays the header lines from a VCF file.\n\n"
              << "Example:\n"
              << "  ./VCFX_header_parser < input.vcf > header.txt\n";
}

void processHeader(std::istream& in, std::ostream& out) {
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '#') {
            out << line << std::endl;
        } else {
            break; // Stop reading after header
        }
    }
}

int main(int argc, char* argv[]) {
    // Simple argument parsing
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // Process header
    processHeader(std::cin, std::cout);
    return 0;
}


./VCFX_header_parser/VCFX_header_parser.h
#ifndef VCFX_HEADER_PARSER_H
#define VCFX_HEADER_PARSER_H

#include <iostream>
#include <string>

// Function to process and extract header lines from VCF
void processHeader(std::istream& in, std::ostream& out);

#endif // VCFX_HEADER_PARSER_H


./VCFX_hwe_tester/VCFX_hwe_tester.cpp
#include "VCFX_hwe_tester.h"
#include <getopt.h>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <vector>

// Implementation of VCFX_hwe_tester
int VCFXHWETester::run(int argc, char* argv[]) {
    // Parse command-line arguments (if any in future extensions)
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {0,      0,           0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Perform HWE tests on stdin
    performHWE(std::cin);

    return 0;
}

void VCFXHWETester::displayHelp() {
    std::cout << "VCFX_hwe_tester: Perform Hardy-Weinberg Equilibrium (HWE) tests on a VCF file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_hwe_tester [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_hwe_tester < input.vcf\n";
}

void VCFXHWETester::performHWE(std::istream& in) {
    std::string line;
    // Print header for output
    std::cout << "CHROM\tPOS\tID\tREF\tALT\tHWE_p-value\n";

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue; // Skip header lines

        // Split the VCF line into fields
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() < 10) {
            std::cerr << "Invalid VCF line with fewer than 10 fields.\n";
            continue;
        }

        std::string chrom = fields[0];
        std::string pos = fields[1];
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        std::string format = fields[8];

        // Find the index of the GT field
        std::vector<std::string> format_fields;
        std::stringstream fmt_ss(format);
        std::string fmt_field;
        while (std::getline(fmt_ss, fmt_field, ':')) {
            format_fields.push_back(fmt_field);
        }

        int gt_index = -1;
        for (size_t i = 0; i < format_fields.size(); ++i) {
            if (format_fields[i] == "GT") {
                gt_index = static_cast<int>(i);
                break;
            }
        }

        if (gt_index == -1) {
            std::cerr << "GT field not found in FORMAT column.\n";
            continue;
        }

        // Collect all genotype information
        std::vector<std::string> genotypes;
        for (size_t i = 9; i < fields.size(); ++i) {
            std::vector<std::string> sample_fields;
            std::stringstream samp_ss(fields[i]);
            std::string samp_field;
            while (std::getline(samp_ss, samp_field, ':')) {
                sample_fields.push_back(samp_field);
            }

            if (gt_index >= static_cast<int>(sample_fields.size())) {
                std::cerr << "GT index out of range in sample fields.\n";
                continue;
            }

            std::string genotype = sample_fields[gt_index];
            genotypes.push_back(genotype);
        }

        // Parse genotypes to count homozygotes and heterozygotes
        int homRef = 0;
        int het = 0;
        int homAlt = 0;
        bool valid = parseGenotypes(genotypes, homRef, het, homAlt);
        if (!valid) {
            std::cerr << "Error parsing genotypes for variant at position " << chrom << ":" << pos << "\n";
            continue;
        }

        // Calculate HWE p-value
        double pValue = calculateHWE(homRef, het, homAlt);

        // Output the result
        std::cout << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" 
                  << std::fixed << std::setprecision(5) << pValue << "\n";
    }
}

bool VCFXHWETester::parseGenotypes(const std::vector<std::string>& genotypes, int& homRef, int& het, int& homAlt) {
    homRef = 0;
    het = 0;
    homAlt = 0;

    for (const auto& gt : genotypes) {
        if (gt.empty() || gt == "." || gt == "./." || gt == ".|.") continue;

        std::string gt_clean = gt;
        // Replace '|' with '/' for consistency
        std::replace(gt_clean.begin(), gt_clean.end(), '|', '/');

        // Split genotype into alleles
        std::vector<std::string> alleles;
        std::stringstream ss(gt_clean);
        std::string allele;
        while (std::getline(ss, allele, '/')) {
            alleles.push_back(allele);
        }

        if (alleles.size() != 2) {
            // Not a diploid genotype
            return false;
        }

        if (alleles[0] == "0" && alleles[1] == "0") {
            homRef++;
        } else if (alleles[0] == "1" && alleles[1] == "1") {
            homAlt++;
        } else if ((alleles[0] == "0" && alleles[1] == "1") || (alleles[0] == "1" && alleles[1] == "0")) {
            het++;
        } else {
            // Invalid allele
            return false;
        }
    }

    return true;
}

double VCFXHWETester::calculateHWE(int homRef, int het, int homAlt) {
    int N = homRef + het + homAlt;
    int x = 2 * homAlt + het; // Count of alternate alleles
    int y = 2 * homRef + het; // Count of reference alleles

    // Ensure the site is biallelic
    if (x + y != 2 * N) {
        return 1.0; // Return 1.0 for non-biallelic sites
    }

    return genotypeProbability(homRef, het, homAlt);
}

double VCFXHWETester::genotypeProbability(int homRef, int het, int homAlt) {
    int N = homRef + het + homAlt;
    int x = 2 * homAlt + het; // Count of alternate alleles
    int y = 2 * homRef + het; // Count of reference alleles

    // Initialize variables for exact test
    int n = y; // Number of ref alleles
    int N_total = 2 * N;

    // Initialize log factorials
    // Precompute log factorials up to 2N
    std::vector<double> logFac(2 * N + 1, 0.0);
    for (int i = 1; i <= 2 * N; ++i) {
        logFac[i] = logFac[i - 1] + std::log(static_cast<double>(i));
    }

    // Function to calculate log(n!), using precomputed values
    auto logFactorial = [&](int n) -> double {
        if (n < 0 || n > 2 * N) return 0.0; // Undefined
        return logFac[n];
    };

    // Function to calculate the log probability of a genotype configuration
    auto logProb = [&](int a) -> double {
        // a = number of heterozygotes
        // number of homozygous refs = (y - a) / 2
        // number of homozygous alts = (x - a) / 2
        if ((y - a) < 0 || (x - a) < 0) return -INFINITY;
        if ((y - a) % 2 != 0 || (x - a) % 2 != 0) return -INFINITY;

        int homRef = (y - a) / 2;
        int homAlt = (x - a) / 2;

        // Multinomial coefficient: (N)! / (n_AA! * n_Aa! * n_Aa!)
        double logCoef = logFactorial(N) - logFactorial(homRef) - logFactorial(homAlt) - logFactorial(a);

        // Probability under HWE: (p^2)^homRef * (2pq)^het * (q^2)^homAlt
        // Taking log:
        // 2*log(p) * homRef + log(2) + log(p) + log(q) * het + 2*log(q) * homAlt
        // Simplified as:
        // 2*homRef*log(p) + het*log(2*p*q) + 2*homAlt*log(q)
        // Where p = y / (y + x), q = x / (y + x)
        double p = static_cast<double>(y) / (y + x);
        double q = static_cast<double>(x) / (y + x);

        if (p == 0.0 || q == 0.0) return -INFINITY;

        double logP = 2.0 * homRef * std::log(p) + a * std::log(2.0 * p * q) + 2.0 * homAlt * std::log(q);

        return logCoef + logP;
    };

    // Calculate observed heterozygotes
    int observedHet = het;

    // Calculate log probability of observed genotype
    double logProbObs = logProb(observedHet);

    // Find minimum log probability to avoid underflow
    double minLogProb = logProbObs;
    for (int a = 0; a <= std::min(y, x); ++a) {
        double lp = logProb(a);
        if (lp < minLogProb && lp != -INFINITY) {
            minLogProb = lp;
        }
    }

    // Calculate the exponents relative to minLogProb for numerical stability
    double sum = 0.0;
    double probObs = 0.0;

    for (int a = 0; a <= std::min(y, x); ++a) {
        double lp = logProb(a);
        if (lp == -INFINITY) continue;

        double relative = lp - minLogProb;
        double prob = std::exp(relative);
        sum += prob;

        if (a == observedHet) {
            probObs = prob;
        }
    }

    // Normalize
    double probNormalized = probObs / sum;

    // Calculate p-value: sum of probabilities <= probObs
    double pValue = 0.0;
    for (int a = 0; a <= std::min(y, x); ++a) {
        double lp = logProb(a);
        if (lp == -INFINITY) continue;

        double relative = lp - minLogProb;
        double prob = std::exp(relative);

        if (prob <= probObs + 1e-8) { // Added a small epsilon for floating point precision
            pValue += prob;
        }
    }

    // Normalize p-value
    pValue /= sum;

    // Ensure p-value is between 0 and 1
    if (pValue > 1.0) pValue = 1.0;
    if (pValue < 0.0) pValue = 0.0;

    return pValue;
}

int main(int argc, char* argv[]) {
    VCFXHWETester tester;
    return tester.run(argc, argv);
}


./VCFX_hwe_tester/VCFX_hwe_tester.h
#ifndef VCFX_HWE_TESTER_H
#define VCFX_HWE_TESTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFX_hwe_tester: Header file for Hardy-Weinberg Equilibrium (HWE) testing tool
class VCFXHWETester {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Processes VCF input and performs HWE tests
    void performHWE(std::istream& in);

    // Parses genotype information to count genotypes
    bool parseGenotypes(const std::vector<std::string>& genotypes, int& homRef, int& het, int& homAlt);

    // Calculates HWE p-value using the exact test
    double calculateHWE(int homRef, int het, int homAlt);

    // Computes the probability of a given genotype count
    double genotypeProbability(int homRef, int het, int homAlt);

    // Computes log factorial to aid in probability calculations
    double logFactorial(int n);

    // Log-sum-exp trick for numerical stability
    double logSumExp(const std::vector<double>& logProbs);
};

#endif // VCFX_HWE_TESTER_H


./VCFX_impact_filter/VCFX_impact_filter.cpp
#include "VCFX_impact_filter.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <regex>

// Implementation of VCFXImpactFilter
int VCFXImpactFilter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string targetImpact;

    static struct option long_options[] = {
        {"help",        no_argument,       0, 'h'},
        {"filter-impact", required_argument, 0, 'i'},
        {0,             0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hi:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'i':
                targetImpact = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || targetImpact.empty()) {
        displayHelp();
        return 1;
    }

    // Perform impact filtering on stdin and output to stdout
    filterByImpact(std::cin, std::cout, targetImpact);

    return 0;
}

void VCFXImpactFilter::displayHelp() {
    std::cout << "VCFX_impact_filter: Filter VCF variants based on predicted impact from annotations.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_impact_filter --filter-impact \"<IMPACT_LEVEL>\" [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                 Display this help message and exit\n";
    std::cout << "  -i, --filter-impact <level> Specify the impact level to filter (e.g., HIGH, MODERATE)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_impact_filter --filter-impact \"HIGH\" < input.vcf > filtered.vcf\n";
}

void VCFXImpactFilter::filterByImpact(std::istream& in, std::ostream& out, const std::string& targetImpact) {
    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerFields;
    size_t impactIndex = std::string::npos;

    // Define possible impact levels in hierarchy
    enum ImpactLevel { UNKNOWN, MODIFIER, LOW, MODERATE, HIGH };
    ImpactLevel targetLevel;

    if (targetImpact == "HIGH") {
        targetLevel = HIGH;
    } else if (targetImpact == "MODERATE") {
        targetLevel = MODERATE;
    } else if (targetImpact == "LOW") {
        targetLevel = LOW;
    } else if (targetImpact == "MODIFIER") {
        targetLevel = MODIFIER;
    } else {
        std::cerr << "Error: Invalid impact level \"" << targetImpact << "\". Choose from HIGH, MODERATE, LOW, MODIFIER.\n";
        return;
    }

    // Regex to extract Impact=LEVEL from INFO field
    std::regex impactRegex(R"(Impact=([A-Z]+))");

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            // Handle header lines
            if (line.substr(0, 6) == "#CHROM") {
                // Parse header to find the index of 'Impact' in INFO field if necessary
                out << line << "\tREF_IMPACT\n"; // Add a new INFO field for filtered impact
                headerParsed = true;
            } else {
                // Other header lines
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
        std::string field;
        std::vector<std::string> fieldsVec;

        while (std::getline(ss, field, '\t')) {
            fieldsVec.push_back(field);
        }

        if (fieldsVec.size() < 8) {
            std::cerr << "Warning: Invalid VCF line with fewer than 8 fields: " << line << "\n";
            continue;
        }

        std::string infoField = fieldsVec[7];
        std::smatch matches;
        std::string impactValue = "UNKNOWN";

        if (std::regex_search(infoField, matches, impactRegex)) {
            impactValue = matches[1];
        }

        // Determine the impact level
        ImpactLevel variantLevel;
        if (impactValue == "HIGH") {
            variantLevel = HIGH;
        } else if (impactValue == "MODERATE") {
            variantLevel = MODERATE;
        } else if (impactValue == "LOW") {
            variantLevel = LOW;
        } else if (impactValue == "MODIFIER") {
            variantLevel = MODIFIER;
        } else {
            variantLevel = UNKNOWN;
        }

        // Define hierarchy: HIGH > MODERATE > LOW > MODIFIER
        // Include variants >= targetImpact level
        bool includeVariant = false;
        switch (targetLevel) {
            case HIGH:
                if (variantLevel == HIGH) includeVariant = true;
                break;
            case MODERATE:
                if (variantLevel == HIGH || variantLevel == MODERATE) includeVariant = true;
                break;
            case LOW:
                if (variantLevel == HIGH || variantLevel == MODERATE || variantLevel == LOW) includeVariant = true;
                break;
            case MODIFIER:
                includeVariant = true; // Include all
                break;
            default:
                includeVariant = false;
        }

        if (includeVariant) {
            out << line << "\t" << impactValue << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXImpactFilter impactFilter;
    return impactFilter.run(argc, argv);
}   

./VCFX_impact_filter/VCFX_impact_filter.h
#ifndef VCFX_IMPACT_FILTER_H
#define VCFX_IMPACT_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXImpactFilter: Header file for Variant Impact Filter tool
class VCFXImpactFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input based on the specified impact level
    void filterByImpact(std::istream& in, std::ostream& out, const std::string& targetImpact);
};

#endif // VCFX_IMPACT_FILTER_H


./VCFX_inbreeding_calculator/VCFX_inbreeding_calculator.cpp
#include "VCFX_inbreeding_calculator.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>

// Implementation

int VCFXInbreedingCalculator::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {
        {"help", no_argument,       0, 'h'},
        {0,      0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    // Calculate inbreeding coefficients
    calculateInbreedingCoefficients(std::cin, std::cout);

    return 0;
}

void VCFXInbreedingCalculator::displayHelp() {
    std::cout << "VCFX_inbreeding_calculator: Calculate inbreeding coefficients (F-statistics) for each individual in a population based on VCF genotypes.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_inbreeding_calculator [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help               Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_inbreeding_calculator < input.vcf > inbreeding_coefficients.txt\n";
}

bool VCFXInbreedingCalculator::parseGenotype(const std::string& genotype, int& a1, int& a2) {
    // Parses genotype in the form GT:... e.g., 0/1, 1/1, 0/0, ./.
    if (genotype.empty()) return false;

    size_t sep = genotype.find(':');
    std::string gt = (sep != std::string::npos) ? genotype.substr(0, sep) : genotype;

    // Handle missing data
    if (gt.find('.') != std::string::npos) {
        return false;
    }

    // Split alleles
    size_t slash = gt.find('/');
    size_t pipe = gt.find('|');
    if (slash != std::string::npos) {
        a1 = std::stoi(gt.substr(0, slash));
        a2 = std::stoi(gt.substr(slash + 1));
        return true;
    }
    else if (pipe != std::string::npos) {
        a1 = std::stoi(gt.substr(0, pipe));
        a2 = std::stoi(gt.substr(pipe + 1));
        return true;
    }

    return false;
}

double VCFXInbreedingCalculator::calculateExpectedHet(int totalAlleles, int homRef, int homAlt, int het) {
    int refAlleles = (homRef * 2) + het;
    int altAlleles = (homAlt * 2) + het;
    double p = static_cast<double>(refAlleles) / totalAlleles;
    double q = static_cast<double>(altAlleles) / totalAlleles;
    return 2.0 * p * q;
}

double VCFXInbreedingCalculator::calculateF(int homozygous, int heterozygous, double expectedHet) {
    if (expectedHet == 0.0) return 0.0;
    return 1.0 - static_cast<double>(heterozygous) / expectedHet;
}

void VCFXInbreedingCalculator::calculateInbreedingCoefficients(std::istream& in, std::ostream& out) {
    std::string line;
    std::vector<std::string> sampleNames;
    bool headerParsed = false;

    // Initialize counts
    std::unordered_map<std::string, int> homozygousCounts;
    std::unordered_map<std::string, int> heterozygousCounts;
    std::unordered_map<std::string, int> totalCounts;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line.substr(0, 6) == "#CHROM") {
            std::stringstream ss(line);
            std::string field;
            // Parse header fields
            std::vector<std::string> headers;
            while (std::getline(ss, field, '\t')) {
                headers.push_back(field);
            }
            // Sample names start from the 10th column
            for (size_t i = 9; i < headers.size(); ++i) {
                sampleNames.push_back(headers[i]);
                homozygousCounts[headers[i]] = 0;
                heterozygousCounts[headers[i]] = 0;
                totalCounts[headers[i]] = 0;
            }
            headerParsed = true;
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header line with #CHROM not found.\n";
            return;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info, format;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        // Read genotype fields if present
        std::vector<std::string> genotypes;
        std::string genotype;
        while (ss >> genotype) {
            genotypes.push_back(genotype);
        }

        for (size_t i = 0; i < genotypes.size() && i < sampleNames.size(); ++i) {
            int a1, a2;
            if (parseGenotype(genotypes[i], a1, a2)) {
                if (a1 == a2) {
                    homozygousCounts[sampleNames[i]] += 1;
                }
                else {
                    heterozygousCounts[sampleNames[i]] += 1;
                }
                totalCounts[sampleNames[i]] += 1;
            }
        }
    }

    // Calculate F-statistics
    out << "Sample\tInbreeding_Coefficient(F)\n";
    for (const auto& sample : sampleNames) {
        int homo = homozygousCounts[sample];
        int het = heterozygousCounts[sample];
        int total = totalCounts[sample];
        if (total == 0) {
            out << sample << "\tNA\n";
            continue;
        }
        int totalAlleles = total * 2;
        int homRef = homo;
        int homAlt = homo; // Assuming homozygous can be either homRef or homAlt; lacking information
        double expectedHet = calculateExpectedHet(totalAlleles, homRef, homAlt, het);
        double F = calculateF(homo, het, expectedHet);
        out << sample << "\t" << F << "\n";
    }
}

int main(int argc, char* argv[]) {
    VCFXInbreedingCalculator inbreedingCalculator;
    return inbreedingCalculator.run(argc, argv);
}

./VCFX_inbreeding_calculator/VCFX_inbreeding_calculator.h
#ifndef VCFX_INBREEDING_CALCULATOR_H
#define VCFX_INBREEDING_CALCULATOR_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// VCFXInbreedingCalculator: Header file for Inbreeding Coefficient Calculator Tool
class VCFXInbreedingCalculator {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Calculates inbreeding coefficients (F-statistics) from VCF input
    void calculateInbreedingCoefficients(std::istream& in, std::ostream& out);

    // Parses genotype and returns allele counts
    bool parseGenotype(const std::string& genotype, int& a1, int& a2);

    // Calculates Hardy-Weinberg Expected Heterozygosity
    double calculateExpectedHet(int totalAlleles, int homRef, int homAlt, int het);

    // Calculates F-statistic for an individual
    double calculateF(int homozygous, int heterozygous, double expectedHet);
};

#endif // VCFX_INBREEDING_CALCULATOR_H


./VCFX_indel_normalizer/VCFX_indel_normalizer.cpp
#include "VCFX_indel_normalizer.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>

// Implementation of VCFXIndelNormalizer
int VCFXIndelNormalizer::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {
        {"help",            no_argument,       0, 'h'},
        {0,                 0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    // Perform indel normalization on stdin and output to stdout
    normalizeIndels(std::cin, std::cout);

    return 0;
}

void VCFXIndelNormalizer::displayHelp() {
    std::cout << "VCFX_indel_normalizer: Normalize indels to their left-most representation.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_indel_normalizer [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_indel_normalizer < input.vcf > normalized.vcf\n";
}

bool VCFXIndelNormalizer::isIndel(const std::string& ref, const std::string& alt) {
    return (ref.length() != alt.length());
}

bool VCFXIndelNormalizer::normalizeVariant(std::string& chrom, std::string& pos, std::string& ref, std::string& alt) {
    // Left-align the indel as per VCF specifications
    size_t start = 0;
    while (start < ref.length() && start < alt.length() && ref[ref.length() - 1 - start] == alt[alt.length() - 1 - start]) {
        start++;
    }

    if (start == 0) {
        return false; // No normalization needed
    }

    ref = ref.substr(0, ref.length() - start);
    alt = alt.substr(0, alt.length() - start);
    
    // Ensure there is at least one base before the indel
    if (ref.empty() || alt.empty()) {
        return false;
    }

    return true;
}

void VCFXIndelNormalizer::normalizeIndels(std::istream& in, std::ostream& out) {
    std::string line;
    bool headerPassed = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            out << line << "\n";
            if (line.substr(0, 6) == "#CHROM") {
                headerPassed = true;
            }
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info, format;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        if (isIndel(ref, alt)) {
            if (normalizeVariant(chrom, pos, ref, alt)) {
                // Reconstruct the VCF line with normalized indel
                std::string rest_of_line;
                getline(ss, rest_of_line);
                out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t"
                    << qual << "\t" << filter << "\t" << info << "\t" << format << rest_of_line << "\n";
                continue;
            }
        }

        // Output the original line if no normalization is needed
        out << line << "\n";
    }
}

int main(int argc, char* argv[]) {
    VCFXIndelNormalizer indelNormalizer;
    return indelNormalizer.run(argc, argv);
}

./VCFX_indel_normalizer/VCFX_indel_normalizer.h
#ifndef VCFX_INDEL_NORMALIZER_H
#define VCFX_INDEL_NORMALIZER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXIndelNormalizer: Header file for Indel Normalization Tool
class VCFXIndelNormalizer {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Normalizes indels in the VCF input
    void normalizeIndels(std::istream& in, std::ostream& out);

    // Determines if a variant is an indel
    bool isIndel(const std::string& ref, const std::string& alt);

    // Normalizes a single indel variant
    bool normalizeVariant(std::string& chrom, std::string& pos, std::string& ref, std::string& alt);
};

#endif // VCFX_INDEL_NORMALIZER_H


./VCFX_indexer/VCFX_indexer.cpp
#include "VCFX_indexer.h"
#include <sstream>
#include <vector>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_indexer\n"
              << "Usage: VCFX_indexer [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h  Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Creates an index for a VCF file, mapping each variant's chromosome and position to its byte offset in the file.\n\n"
              << "Example:\n"
              << "  ./VCFX_indexer < input.vcf > index.tsv\n";
}

// Function to split a string by a delimiter and return a vector of tokens
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::stringstream ss(s);
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Function to create an index for VCF file
void createVCFIndex(std::istream& in, std::ostream& out) {
    std::string line;
    std::string header;
    long offset = 0;
    bool header_found = false;

    // Print header for index
    out << "CHROM\tPOS\tFILE_OFFSET\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            offset += (line.length() + 1); // +1 for newline character
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                header_found = true;
            }
            offset += (line.length() + 1);
            continue; // Skip header lines
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return;
        }

        std::vector<std::string> fields = split(line, '\t');
        if (fields.size() < 2) {
            std::cerr << "Warning: Skipping invalid VCF line with fewer than 2 fields.\n";
            offset += (line.length() + 1);
            continue;
        }

        std::string chrom = fields[0];
        int pos = 0;
        try {
            pos = std::stoi(fields[1]);
        } catch (...) {
            std::cerr << "Warning: Invalid POS value on line with CHROM " << chrom << ". Skipping line.\n";
            offset += (line.length() + 1);
            continue;
        }

        // Output the index entry
        out << chrom << "\t" << pos << "\t" << offset << "\n";

        // Update the offset
        offset += (line.length() + 1);
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

    // Display help if no input is provided
    if (argc == 1) {
        printHelp();
        return 1;
    }

    // Create VCF index
    createVCFIndex(std::cin, std::cout);
    return 0;
}


./VCFX_indexer/VCFX_indexer.h
#ifndef VCFX_INDEXER_H
#define VCFX_INDEXER_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// Structure to hold variant information for indexing
struct VariantIndex {
    std::string chrom;
    int pos;
    long file_offset; // Byte offset in the file
};

// Function to display help message
void printHelp();

// Function to create an index for VCF file
void createVCFIndex(std::istream& in, std::ostream& out);

#endif // VCFX_INDEXER_H


./VCFX_info_aggregator/VCFX_info_aggregator.cpp

#include "VCFX_info_aggregator.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <map>
#include <vector>

// Implementation of VCFXInfoAggregator
int VCFXInfoAggregator::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string infoFieldsStr;

    static struct option long_options[] = {
        {"help",          no_argument,       0, 'h'},
        {"aggregate-info", required_argument, 0, 'a'},
        {0,               0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "ha:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'a':
                infoFieldsStr = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || infoFieldsStr.empty()) {
        displayHelp();
        return 1;
    }

    // Split infoFieldsStr by comma to get individual INFO fields
    std::vector<std::string> infoFields;
    std::stringstream ss(infoFieldsStr);
    std::string field;
    while (std::getline(ss, field, ',')) {
        // Trim whitespace
        field.erase(field.find_last_not_of(" \n\r\t")+1);
        field.erase(0, field.find_first_not_of(" \n\r\t"));
        if (!field.empty()) {
            infoFields.push_back(field);
        }
    }

    if (infoFields.empty()) {
        std::cerr << "Error: No valid INFO fields specified for aggregation.\n";
        return 1;
    }

    // Perform INFO field aggregation on stdin and output to stdout
    aggregateInfo(std::cin, std::cout, infoFields);

    return 0;
}

void VCFXInfoAggregator::displayHelp() {
    std::cout << "VCFX_info_aggregator: Aggregate numeric values in the INFO field across samples.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_info_aggregator --aggregate-info \"<INFO_FIELDS>\" [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                       Display this help message and exit\n";
    std::cout << "  -a, --aggregate-info <fields>    Comma-separated list of INFO fields to aggregate (e.g., DP,AF)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_info_aggregator --aggregate-info \"DP,AF\" < input.vcf > aggregated_info.txt\n";
}

void VCFXInfoAggregator::aggregateInfo(std::istream& in, std::ostream& out, const std::vector<std::string>& infoFields) {
    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerFields;
    std::map<std::string, std::vector<double>> aggregates; // Key: INFO field, Value: list of numeric values

    // Initialize aggregates map
    for (const auto& field : infoFields) {
        aggregates[field] = std::vector<double>();
    }

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            // Handle header lines
            if (line.substr(0, 6) == "#CHROM") {
                // Output header with additional aggregated fields
                out << line;
                for (const auto& field : infoFields) {
                    out << "\t" << "AGG_" << field;
                }
                out << "\n";
                headerParsed = true;
            } else {
                // Other header lines
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
        std::string fieldEntry;
        std::vector<std::string> fieldsVec;

        while (std::getline(ss, fieldEntry, '\t')) {
            fieldsVec.push_back(fieldEntry);
        }

        if (fieldsVec.size() < 8) {
            std::cerr << "Warning: Invalid VCF line with fewer than 8 fields: " << line << "\n";
            continue;
        }

        std::string infoField = fieldsVec[7];
        std::map<std::string, double> currentVariantAggregates;

        // Extract specified INFO fields
        for (const auto& aggField : infoFields) {
            size_t pos = infoField.find(aggField + "=");
            if (pos != std::string::npos) {
                size_t start = pos + aggField.length() + 1;
                size_t end = infoField.find(';', start);
                std::string valueStr = (end != std::string::npos) ? infoField.substr(start, end - start) : infoField.substr(start);
                try {
                    double value = std::stod(valueStr);
                    aggregates[aggField].push_back(value);
                    currentVariantAggregates[aggField] = value;
                } catch (...) {
                    std::cerr << "Warning: Non-numeric value for INFO field \"" << aggField << "\": " << valueStr << "\n";
                    currentVariantAggregates[aggField] = 0.0;
                }
            } else {
                // INFO field not present; assign default value
                aggregates[aggField].push_back(0.0);
                currentVariantAggregates[aggField] = 0.0;
            }
        }

        // Append aggregated values to the VCF line
        out << line;
        for (const auto& aggField : infoFields) {
            out << "\t" << currentVariantAggregates[aggField];
        }
        out << "\n";
    }

    std::cout << "Aggregated INFO Fields:\n";
    for (const auto& agg : aggregates) {
        double sum = 0.0;
        for (const auto& val : agg.second) {
            sum += val;
        }
        double average = (agg.second.empty()) ? 0.0 : sum / agg.second.size();
        std::cout << agg.first << ": Sum = " << sum << ", Average = " << average << "\n";
    }
}

int main(int argc, char* argv[]) {
    VCFXInfoAggregator infoAggregator;
    return infoAggregator.run(argc, argv);
}

./VCFX_info_aggregator/VCFX_info_aggregator.h
#ifndef VCFX_INFO_AGGREGATOR_H
#define VCFX_INFO_AGGREGATOR_H
#include <iostream>
#include <string>
#include <vector>
// VCFXInfoAggregator: Header file for INFO Field Aggregator tool
class VCFXInfoAggregator {
public:

    // Entry point for the tool
    int run(int argc, char* argv[]);

private:

    // Displays the help message
    void displayHelp();

    // Aggregates specified INFO fields across samples
    void aggregateInfo(std::istream& in, std::ostream& out, const std::vector<std::string>& infoFields);
};
#endif // VCFX_INFO_AGGREGATOR_H

./VCFX_info_parser/VCFX_info_parser.cpp
#include "VCFX_info_parser.h"
#include <sstream>
#include <algorithm>
#include <unordered_map>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_info_parser\n"
              << "Usage: VCFX_info_parser [OPTIONS]\n\n"
              << "Options:\n"
              << "  --info, -i \"FIELD1,FIELD2\"   Specify the INFO fields to display (e.g., \"DP,AF\").\n"
              << "  --help, -h                    Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Parses the INFO field of a VCF file and displays the selected INFO fields in a user-friendly format.\n\n"
              << "Examples:\n"
              << "  ./VCFX_info_parser --info \"DP,AF\" < input.vcf > output_info.tsv\n";
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
    return false;
}

// Function to split a string by a delimiter and return a vector of tokens
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::stringstream ss(s);
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Function to parse the INFO field and display selected fields
bool parseInfoFields(std::istream& in, std::ostream& out, const std::vector<std::string>& info_fields) {
    std::string line;
    bool header_printed = false;

    // Print header
    if (!info_fields.empty()) {
        out << "CHROM\tPOS\tID\tREF\tALT";
        for (const auto& field : info_fields) {
            out << "\t" << field;
        }
        out << "\n";
    }

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            // Skip header lines except #CHROM
            if (line.find("#CHROM") == 0) {
                // Optionally, you could include parsing and displaying header information
            }
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info)) {
            std::cerr << "Warning: Skipping invalid VCF line.\n";
            continue;
        }

        // Parse INFO field into key-value pairs
        std::vector<std::string> info_entries = split(info, ';');
        std::unordered_map<std::string, std::string> info_map;
        for (const auto& entry : info_entries) {
            size_t eq = entry.find('=');
            if (eq != std::string::npos) {
                std::string key = entry.substr(0, eq);
                std::string value = entry.substr(eq + 1);
                info_map[key] = value;
            } else {
                // Flag without a value
                info_map[entry] = "";
            }
        }

        // Output selected INFO fields
        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt;
        for (const auto& field : info_fields) {
            auto it = info_map.find(field);
            if (it != info_map.end()) {
                out << "\t" << it->second;
            } else {
                out << "\t.";
            }
        }
        out << "\n";
    }

    return true;
}

int main(int argc, char* argv[]) {
    std::vector<std::string> info_fields;

    // Argument parsing
    if (!parseArguments(argc, argv, info_fields)) {
        std::cerr << "Error: INFO fields not specified.\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    // Parse and display INFO fields
    bool success = parseInfoFields(std::cin, std::cout, info_fields);
    return success ? 0 : 1;
}


./VCFX_info_parser/VCFX_info_parser.h
#ifndef VCFX_INFO_PARSER_H
#define VCFX_INFO_PARSER_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], std::vector<std::string>& info_fields);

// Function to parse the INFO field and display selected fields
bool parseInfoFields(std::istream& in, std::ostream& out, const std::vector<std::string>& info_fields);

#endif // VCFX_INFO_PARSER_H


./VCFX_info_summarizer/VCFX_info_summarizer.cpp
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


./VCFX_info_summarizer/VCFX_info_summarizer.h
#ifndef VCFX_INFO_SUMMARIZER_H
#define VCFX_INFO_SUMMARIZER_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// Function to display help message
void printHelp();

// Structure to hold statistical summaries
struct StatSummary {
    double mean;
    double median;
    double mode;
};

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], std::vector<std::string>& info_fields);

// Function to calculate mean
double calculateMean(const std::vector<double>& data);

// Function to calculate median
double calculateMedian(std::vector<double> data);

// Function to calculate mode
double calculateMode(const std::vector<double>& data);

// Function to parse the INFO field and collect specified fields
bool summarizeInfoFields(std::istream& in, std::ostream& out, const std::vector<std::string>& info_fields);

#endif // VCFX_INFO_SUMMARIZER_H


./VCFX_ld_calculator/VCFX_ld_calculator.cpp
#include "VCFX_ld_calculator.h"

int main(int argc, char* argv[]) {
    // TODO: Implement VCFX_ld_calculator functionality
    return 0;
}


./VCFX_ld_calculator/VCFX_ld_calculator.h
#ifndef VCFX_LD_CALCULATOR_H
#define VCFX_LD_CALCULATOR_H

// Declarations for VCFX_ld_calculator

#endif // VCFX_LD_CALCULATOR_H


./VCFX_merger/VCFX_merger.cpp
#include "VCFX_merger.h"
#include <getopt.h>
#include <fstream>
#include <algorithm>
#include <map>

// Implementation of VCFX_merger
int VCFXMerger::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::vector<std::string> inputFiles;

    static struct option long_options[] = {
        {"merge", required_argument, 0, 'm'},
        {"help",  no_argument,       0, 'h'},
        {0,       0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "m:h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'm':
                {
                    // Split comma-separated file names
                    std::string files = optarg;
                    size_t pos = 0;
                    while ((pos = files.find(',')) != std::string::npos) {
                        inputFiles.emplace_back(files.substr(0, pos));
                        files.erase(0, pos + 1);
                    }
                    inputFiles.emplace_back(files);
                }
                break;
            case 'h':
            default:
                showHelp = true;
        }
    }

    if (showHelp || inputFiles.empty()) {
        displayHelp();
        return 0;
    }

    // Merge VCF files and output to stdout
    mergeVCF(inputFiles, std::cout);

    return 0;
}

void VCFXMerger::displayHelp() {
    std::cout << "VCFX_merger: Merge multiple VCF files by variant position.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_merger --merge file1.vcf,file2.vcf,... [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -m, --merge    Comma-separated list of VCF files to merge\n";
    std::cout << "  -h, --help     Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_merger --merge sample1.vcf,sample2.vcf > merged.vcf\n";
}

void VCFXMerger::mergeVCF(const std::vector<std::string>& inputFiles, std::ostream& out) {
    std::vector<std::vector<std::string>> allVariants;
    std::vector<std::string> allHeaders;

    for (const auto& file : inputFiles) {
        std::vector<std::vector<std::string>> variants;
        std::vector<std::string> headerLines;
        parseVCF(file, variants, headerLines);

        // Store headers (assuming all headers are identical)
        if (allHeaders.empty()) {
            allHeaders = headerLines;
        }

        allVariants.insert(allVariants.end(), variants.begin(), variants.end());
    }

    // Sort all variants by chromosome and position
    std::sort(allVariants.begin(), allVariants.end(), 
        [this](const std::vector<std::string>& a, const std::vector<std::string>& b) -> bool {
            if (a[0] == b[0]) {
                return std::stoi(a[1]) < std::stoi(b[1]);
            }
            return a[0] < b[0];
        }
    );

    // Output headers
    for (const auto& header : allHeaders) {
        out << header << "\n";
    }

    // Output merged variants
    for (const auto& variant : allVariants) {
        for (size_t i = 0; i < variant.size(); ++i) {
            out << variant[i];
            if (i < variant.size() - 1) {
                out << "\t";
            }
        }
        out << "\n";
    }
}

void VCFXMerger::parseVCF(const std::string& filename, std::vector<std::vector<std::string>>& variants, std::vector<std::string>& headerLines) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            headerLines.push_back(line);
            continue;
        }

        // Split the VCF line into fields
        std::vector<std::string> fields;
        std::string field;
        size_t pos = 0;
        while ((pos = line.find('\t')) != std::string::npos) {
            field = line.substr(0, pos);
            fields.push_back(field);
            line.erase(0, pos + 1);
        }
        fields.push_back(line); // Add the last field

        variants.push_back(fields);
    }

    infile.close();
}

bool VCFXMerger::compareVariants(const std::vector<std::string>& a, const std::vector<std::string>& b) {
    if (a[0] != b[0]) {
        return a[0] < b[0];
    }
    return std::stoi(a[1]) < std::stoi(b[1]);
}

int main(int argc, char* argv[]) {
    VCFXMerger merger;
    return merger.run(argc, argv);
}


./VCFX_merger/VCFX_merger.h
#ifndef VCFX_MERGER_H
#define VCFX_MERGER_H

#include <iostream>
#include <string>
#include <vector>

// VCFX_merger: Header file for VCF file merging tool
class VCFXMerger {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Processes and merges VCF files
    void mergeVCF(const std::vector<std::string>& inputFiles, std::ostream& out);

    // Parses a VCF file and stores variants
    void parseVCF(const std::string& filename, std::vector<std::vector<std::string>>& variants, std::vector<std::string>& headerLines);

    // Compares variants based on chromosome and position
    bool compareVariants(const std::vector<std::string>& a, const std::vector<std::string>& b);
};

#endif // VCFX_MERGER_H


./VCFX_metadata_summarizer/VCFX_metadata_summarizer.cpp
#include "VCFX_metadata_summarizer.h"
#include <getopt.h>
#include <sstream>

// Implementation of VCFX_metadata_summarizer
int VCFXMetadataSummarizer::run(int argc, char* argv[]) {
    // Parse command-line arguments (if any in future extensions)
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {0,      0,           0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Summarize metadata from stdin
    summarizeMetadata(std::cin);

    return 0;
}

void VCFXMetadataSummarizer::displayHelp() {
    std::cout << "VCFX_metadata_summarizer: Summarize key metadata from a VCF file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_metadata_summarizer [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_metadata_summarizer < input.vcf\n";
}

void VCFXMetadataSummarizer::summarizeMetadata(std::istream& in) {
    std::string line;
    std::map<std::string, int> metadata;
    int numSamples = 0;
    int numVariants = 0;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (line.find("##") == 0) {
                // Meta-information lines
                parseHeader(line, metadata);
            } else if (line.find("#CHROM") == 0) {
                // Header line with column names and sample information
                std::vector<std::string> fields;
                std::stringstream ss(line);
                std::string field;
                while (std::getline(ss, field, '\t')) {
                    fields.push_back(field);
                }
                numSamples = static_cast<int>(fields.size()) - 9; // Standard VCF has 9 fixed columns
                metadata["Number of Samples"] = numSamples;
            }
            continue;
        }

        // Count variants
        numVariants++;
    }

    metadata["Number of Variants"] = numVariants;

    // Additional summaries can be added here (e.g., contigs, variant types)
    
    printSummary(metadata);
}

void VCFXMetadataSummarizer::parseHeader(const std::string& line, std::map<std::string, int>& metadata) {
    // Example meta-information:
    // ##fileformat=VCFv4.2
    // ##source=...
    // ##contig=<ID=1,length=248956422>
    if (line.find("##contig=") != std::string::npos) {
        metadata["Number of Contigs"]++;
    }
    if (line.find("##INFO=") != std::string::npos) {
        metadata["Number of INFO Fields"]++;
    }
    if (line.find("##FILTER=") != std::string::npos) {
        metadata["Number of FILTER Fields"]++;
    }
    if (line.find("##FORMAT=") != std::string::npos) {
        metadata["Number of FORMAT Fields"]++;
    }
}

void VCFXMetadataSummarizer::printSummary(const std::map<std::string, int>& metadata) {
    std::cout << "VCF Metadata Summary:\n";
    std::cout << "---------------------\n";
    for (const auto& [key, value] : metadata) {
        std::cout << key << ": " << value << "\n";
    }
}

int main(int argc, char* argv[]) {
    VCFXMetadataSummarizer summarizer;
    return summarizer.run(argc, argv);
}


./VCFX_metadata_summarizer/VCFX_metadata_summarizer.h
#ifndef VCFX_METADATA_SUMMARIZER_H
#define VCFX_METADATA_SUMMARIZER_H

#include <iostream>
#include <string>
#include <vector>
#include <map>

// VCFX_metadata_summarizer: Header file for VCF metadata summarization tool
class VCFXMetadataSummarizer {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Processes VCF input and summarizes metadata
    void summarizeMetadata(std::istream& in);

    // Parses the header lines to extract metadata
    void parseHeader(const std::string& line, std::map<std::string, int>& metadata);

    // Prints the metadata summary
    void printSummary(const std::map<std::string, int>& metadata);
};

#endif // VCFX_METADATA_SUMMARIZER_H


./VCFX_missing_data_handler/VCFX_missing_data_handler.cpp
#include "VCFX_missing_data_handler.h"
#include <sstream>
#include <algorithm>
#include <stdexcept>

/**
 * @brief Displays the help message for the missing data handler tool.
 */
void printHelp() {
    std::cout << "VCFX_missing_data_handler\n"
              << "Usage: VCFX_missing_data_handler [OPTIONS] < input.vcf > output.vcf\n\n"
              << "Options:\n"
              << "  --fill-missing, -f           Impute missing genotypes with a default value (e.g., ./.).\n"
              << "  --default-genotype, -d GEN    Specify the default genotype for imputation (default: ./.).\n"
              << "  --help, -h                    Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Flags or imputes missing genotype data in a VCF file. By default, missing genotypes are flagged, "
              << "but can be imputed with a specified genotype using the --fill-missing option.\n\n"
              << "Examples:\n"
              << "  Flag missing data:\n"
              << "    ./VCFX_missing_data_handler < input.vcf > flagged_output.vcf\n\n"
              << "  Impute missing data with ./. :\n"
              << "    ./VCFX_missing_data_handler --fill-missing --default-genotype \"./.\" < input.vcf > imputed_output.vcf\n";
}

/**
 * @brief Splits a string by a given delimiter.
 *
 * @param str The input string to split.
 * @param delimiter The character delimiter.
 * @return A vector of split substrings.
 */
std::vector<std::string> splitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string temp;
    while (std::getline(ss, temp, delimiter)) {
        tokens.push_back(temp);
    }
    return tokens;
}

/**
 * @brief Parses command-line arguments.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @param args Reference to Arguments structure to populate.
 * @return true if parsing is successful, false otherwise.
 */
bool parseArguments(int argc, char* argv[], Arguments& args) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--fill-missing" || arg == "-f") {
            args.fill_missing = true;
        }
        else if ((arg == "--default-genotype" || arg == "-d") && i + 1 < argc) {
            args.default_genotype = argv[++i];
        }
        else if (arg == "--help" || arg == "-h") {
            printHelp();
            exit(0);
        }
        else {
            // Handle input files or other arguments if necessary
            args.input_files.push_back(arg);
        }
    }
    return true;
}

/**
 * @brief Processes the VCF file to handle missing genotype data.
 *
 * @param in Input stream (VCF file).
 * @param out Output stream (Modified VCF).
 * @param args Command-line arguments specifying behavior.
 * @return true if processing is successful, false otherwise.
 */
bool handleMissingData(std::istream& in, std::ostream& out, const Arguments& args) {
    std::string line;
    std::vector<std::string> header_fields;
    bool header_found = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            out << "\n";
            continue;
        }

        if (line[0] == '#') {
            out << line << "\n";
            if (line.find("#CHROM") == 0) {
                // Parse header to identify sample columns
                header_fields = splitString(line, '\t');
                if (header_fields.size() > 9) {
                    // Samples start from the 10th column (index 9)
                }
                header_found = true;
            }
            continue;
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        std::vector<std::string> fields = splitString(line, '\t');
        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping invalid VCF line with fewer than 9 fields.\n";
            continue;
        }

        // FORMAT column is the 9th field
        std::string format = fields[8];
        std::vector<std::string> format_fields = splitString(format, ':');
        int gt_index = -1;

        // Identify the index of the GT field within the FORMAT
        for (size_t i = 0; i < format_fields.size(); ++i) {
            if (format_fields[i] == "GT") {
                gt_index = static_cast<int>(i);
                break;
            }
        }

        if (gt_index == -1) {
            std::cerr << "Warning: GT field not found in FORMAT column. Skipping line.\n";
            continue;
        }

        // Process each sample
        for (size_t i = 9; i < fields.size(); ++i) {
            std::vector<std::string> sample_data = splitString(fields[i], ':');
            if (gt_index >= static_cast<int>(sample_data.size())) {
                // No GT data for this sample
                continue;
            }

            std::string genotype = sample_data[gt_index];
            bool is_missing = false;

            // Determine if genotype is missing
            if (genotype.empty() || genotype == "." || genotype == "./." || genotype == ".|.") {
                is_missing = true;
            }

            if (is_missing) {
                if (args.fill_missing) {
                    // Impute with default genotype
                    sample_data[gt_index] = args.default_genotype;
                } else {
                    // Flagging missing data can be customized here.
                    // For simplicity, we leave the genotype as is.
                    // Alternatively, annotations can be added to the INFO field.
                }

                // Reconstruct the sample data
                std::stringstream ss;
                for (size_t j = 0; j < sample_data.size(); ++j) {
                    ss << sample_data[j];
                    if (j != sample_data.size() - 1) {
                        ss << ":";
                    }
                }
                fields[i] = ss.str();
            }
        }

        // Reconstruct the VCF line
        std::stringstream ss_line;
        for (size_t i = 0; i < fields.size(); ++i) {
            ss_line << fields[i];
            if (i != fields.size() - 1) {
                ss_line << "\t";
            }
        }
        out << ss_line.str() << "\n";
    }

    return true;
}

/**
 * @brief Main function for the missing data handler tool.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return int Exit status.
 */
int main(int argc, char* argv[]) {
    Arguments args;
    parseArguments(argc, argv, args);

    if (args.fill_missing) {
        std::cerr << "Info: Missing genotypes will be imputed with genotype: " 
                  << args.default_genotype << "\n";
    }
    else {
        std::cerr << "Info: Missing genotypes will be flagged.\n";
    }

    bool success = handleMissingData(std::cin, std::cout, args);
    return success ? 0 : 1;
}


./VCFX_missing_data_handler/VCFX_missing_data_handler.h
#ifndef VCFX_MISSING_DATA_HANDLER_H
#define VCFX_MISSING_DATA_HANDLER_H

#include <iostream>
#include <string>
#include <vector>

/**
 * @brief Structure to hold command-line arguments for the missing data handler tool.
 */
struct Arguments {
    bool fill_missing = false;                 ///< Flag indicating whether to impute missing genotypes.
    std::string default_genotype = "./.";      ///< Default genotype to use for imputation.
    std::vector<std::string> input_files;      ///< List of input VCF files (if any).
};

/**
 * @brief Displays the help message for the missing data handler tool.
 */
void printHelp();

/**
 * @brief Parses command-line arguments.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @param args Reference to Arguments structure to populate.
 * @return true if parsing is successful, false otherwise.
 */
bool parseArguments(int argc, char* argv[], Arguments& args);

/**
 * @brief Splits a string by a given delimiter.
 *
 * @param str The input string to split.
 * @param delimiter The character delimiter.
 * @return A vector of split substrings.
 */
std::vector<std::string> splitString(const std::string& str, char delimiter);

/**
 * @brief Processes the VCF file to handle missing genotype data.
 *
 * @param in Input stream (VCF file).
 * @param out Output stream (Modified VCF).
 * @param args Command-line arguments specifying behavior.
 * @return true if processing is successful, false otherwise.
 */
bool handleMissingData(std::istream& in, std::ostream& out, const Arguments& args);

#endif // VCFX_MISSING_DATA_HANDLER_H


./VCFX_missing_detector/VCFX_missing_detector.cpp
#include "VCFX_missing_detector.h"
#include <getopt.h>
#include <sstream>

int VCFXMissingDetector::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {
        {"help",    no_argument,       0, 'h'},
        {0,         0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    // Perform missing genotype detection
    detectMissingGenotypes(std::cin, std::cout);

    return 0;
}

void VCFXMissingDetector::displayHelp() {
    std::cout << "VCFX_missing_detector: Detect variants with missing sample genotypes.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_missing_detector [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_missing_detector < input.vcf > flagged.vcf\n";
}

void VCFXMissingDetector::detectMissingGenotypes(std::istream& in, std::ostream& out) {
    std::string line;
    bool headerPassed = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            out << line << "\n";
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info, format;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        std::vector<std::string> genotypes;
        std::string genotype;
        while (ss >> genotype) {
            genotypes.push_back(genotype);
        }

        bool hasMissing = false;
        for (const auto& gt : genotypes) {
            if (gt.find("./.") != std::string::npos) {
                hasMissing = true;
                break;
            }
        }

        if (hasMissing) {
            // Append a flag to the INFO field
            if (info.back() != ';' && !info.empty()) {
                info += ";";
            }
            info += "MISSING_GENOTYPES=1";

            // Reconstruct the VCF line with updated INFO
            out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt
                << "\t" << qual << "\t" << filter << "\t" << info << "\t" << format;

            for (const auto& gt : genotypes) {
                out << "\t" << gt;
            }
            out << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXMissingDetector missingDetector;
    return missingDetector.run(argc, argv);
}

./VCFX_missing_detector/VCFX_missing_detector.h
#ifndef VCFX_MISSING_DETECTOR_H
#define VCFX_MISSING_DETECTOR_H

#include <iostream>
#include <string>
#include <vector>

// VCFXMissingDetector: Header file for Missing Sample Detection Tool
class VCFXMissingDetector {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Detects missing genotypes in VCF input
    void detectMissingGenotypes(std::istream& in, std::ostream& out);
};

#endif // VCFX_MISSING_DETECTOR_H


./VCFX_multiallelic_splitter/VCFX_multiallelic_splitter.cpp
#include "VCFX_multiallelic_splitter.h"
#include <sstream>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_multiallelic_splitter\n"
              << "Usage: VCFX_multiallelic_splitter [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h               Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Splits multi-allelic variants in a VCF file into multiple bi-allelic records.\n\n"
              << "Examples:\n"
              << "  ./VCFX_multiallelic_splitter < input.vcf > split_biallelic.vcf\n";
}

// Function to parse a VCF line into VCFVariant structure
bool parseVCFLine(const std::string& line, VCFVariant& variant) {
    if (line.empty() || line[0] == '#') return false;

    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> fields;

    // Split the line by tab
    while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
    }

    // VCF must have at least 8 fields before FORMAT
    if (fields.size() < 8) return false;

    variant.chrom = fields[0];
    try {
        variant.pos = std::stoi(fields[1]);
    } catch (...) {
        return false;
    }
    variant.id = fields[2];
    variant.ref = fields[3];
    // Split ALT alleles by comma
    std::stringstream alt_ss(fields[4]);
    std::string alt_allele;
    while (std::getline(alt_ss, alt_allele, ',')) {
        variant.alt.push_back(alt_allele);
    }
    variant.qual = fields[5];
    variant.filter = fields[6];
    variant.info = fields[7];

    // Collect FORMAT and SAMPLE fields if present
    if (fields.size() > 8) {
        variant.other_fields.assign(fields.begin() + 8, fields.end());
    }

    return true;
}

// Function to reconstruct INFO field
std::string reconstructInfo(const std::string& chrom, int pos, const VCFVariant& variant) {
    return variant.info; // For icity, retain the original INFO field
}

// Function to split multi-allelic variants into bi-allelic
bool splitMultiAllelicVariants(std::istream& in, std::ostream& out) {
    std::string line;
    bool header_found = false;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            // Pass through header lines unchanged
            out << line << "\n";
            if (line.find("#CHROM") == 0) {
                header_found = true;
            }
            continue;
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        VCFVariant variant;
        if (!parseVCFLine(line, variant)) {
            std::cerr << "Warning: Skipping malformed VCF line: " << line << "\n";
            continue;
        }

        // If only one ALT allele, write the record as is
        if (variant.alt.size() <= 1) {
            out << line << "\n";
            continue;
        }

        // More than one ALT allele, split into multiple records
        for (size_t i = 0; i < variant.alt.size(); ++i) {
            VCFVariant split_variant = variant;
            split_variant.alt = { variant.alt[i] };
            split_variant.info = reconstructInfo(variant.chrom, variant.pos, variant);

            // Reconstruct the VCF line
            std::stringstream split_ss;
            split_ss << split_variant.chrom << "\t"
                     << split_variant.pos << "\t"
                     << split_variant.id << "\t"
                     << split_variant.ref << "\t";

            // ALT field
            split_ss << split_variant.alt[0] << "\t"
                     << split_variant.qual << "\t"
                     << split_variant.filter << "\t"
                     << split_variant.info;

            // Append FORMAT and sample fields if present
            if (!split_variant.other_fields.empty()) {
                for (const auto& sample_field : split_variant.other_fields) {
                    split_ss << "\t" << sample_field;
                }
            }

            split_ss << "\n";
            out << split_ss.str();
        }
    }

    return true;
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

    // Split multi-allelic variants
    bool success = splitMultiAllelicVariants(std::cin, std::cout);
    return success ? 0 : 1;
}


./VCFX_multiallelic_splitter/VCFX_multiallelic_splitter.h
#ifndef VCFX_MULTIALLELIC_SPLITTER_H
#define VCFX_MULTIALLELIC_SPLITTER_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Structure to hold VCF variant information
struct VCFVariant {
    std::string chrom;
    int pos;
    std::string id;
    std::string ref;
    std::vector<std::string> alt;
    std::string qual;
    std::string filter;
    std::string info;
    std::vector<std::string> other_fields; // FORMAT and sample data
};

// Function to parse a VCF line into VCFVariant structure
bool parseVCFLine(const std::string& line, VCFVariant& variant);

// Function to split multi-allelic variants into bi-allelic
bool splitMultiAllelicVariants(std::istream& in, std::ostream& out);

#endif // VCFX_MULTIALLELIC_SPLITTER_H


