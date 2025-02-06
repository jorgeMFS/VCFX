./VCFX_af_subsetter/VCFX_af_subsetter.cpp
#include "VCFX_af_subsetter.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <cstdlib>

int VCFXAfSubsetter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    double minAF = 0.0;
    double maxAF = 1.0;

    static struct option long_options[] = {
        {"help",      no_argument,       0, 'h'},
        {"af-filter", required_argument, 0, 'a'},
        {0,           0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "ha:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'a': {
                std::string range(optarg);
                size_t dashPos = range.find('-');
                if (dashPos == std::string::npos) {
                    std::cerr << "Error: Invalid AF range format. Use <minAF>-<maxAF>.\n";
                    displayHelp();
                    return 1;
                }
                try {
                    minAF = std::stod(range.substr(0, dashPos));
                    maxAF = std::stod(range.substr(dashPos + 1));
                    if (minAF < 0.0 || maxAF > 1.0 || minAF > maxAF) {
                        throw std::invalid_argument("AF values out of range or minAF > maxAF.");
                    }
                } catch (const std::invalid_argument&) {
                    std::cerr << "Error: Invalid AF range values. Ensure they are numbers between 0.0 and 1.0 with minAF <= maxAF.\n";
                    displayHelp();
                    return 1;
                }
                break;
            }
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    // Perform AF-based subsetting from stdin to stdout
    subsetByAlleleFrequency(std::cin, std::cout, minAF, maxAF);

    return 0;
}

void VCFXAfSubsetter::displayHelp() {
    std::cout << "VCFX_af_subsetter: Subset variants based on alternate allele frequency (AF) ranges.\n\n"
              << "Usage:\n"
              << "  VCFX_af_subsetter [options] < input.vcf > output.vcf\n\n"
              << "Options:\n"
              << "  -h, --help                     Display this help message and exit\n"
              << "  -a, --af-filter <minAF>-<maxAF>  Specify the AF range for filtering (e.g., 0.01-0.05)\n\n"
              << "Example:\n"
              << "  VCFX_af_subsetter --af-filter 0.01-0.05 < input.vcf > subsetted.vcf\n";
}

bool VCFXAfSubsetter::parseAF(const std::string& infoField, std::vector<double>& afValues) {
    // Find "AF=" in the INFO string
    size_t pos = infoField.find("AF=");
    if (pos == std::string::npos) {
        return false;
    }

    // Extract substring up to the next semicolon or end of string
    size_t start = pos + 3; // skip 'AF='
    size_t end = infoField.find(';', start);
    std::string afStr = (end == std::string::npos)
                        ? infoField.substr(start)
                        : infoField.substr(start, end - start);

    // Split by comma (to handle multi-allelic AF like "0.2,0.8")
    std::stringstream ss(afStr);
    std::string token;
    while (std::getline(ss, token, ',')) {
        try {
            double afVal = std::stod(token);
            afValues.push_back(afVal);
        } catch (...) {
            return false; // parsing error
        }
    }
    return !afValues.empty();
}

void VCFXAfSubsetter::subsetByAlleleFrequency(std::istream& in, std::ostream& out, double minAF, double maxAF) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        // Print header lines (starting with '#') as-is
        if (line[0] == '#') {
            out << line << "\n";
            continue;
        }

        // Split the VCF line into fields
        std::vector<std::string> fields;
        {
            std::stringstream ss(line);
            std::string field;
            while (std::getline(ss, field, '\t')) {
                fields.push_back(field);
            }
        }

        if (fields.size() < 8) {
            std::cerr << "Warning: Skipping invalid VCF line (less than 8 fields): "
                      << line << "\n";
            continue;
        }

        // Parse AF values from the INFO field
        std::string info = fields[7];
        std::vector<double> afValues;
        if (!parseAF(info, afValues)) {
            std::cerr << "Warning: AF not found or invalid in INFO field. Skipping variant: "
                      << line << "\n";
            continue;
        }

        // If any AF is within [minAF, maxAF], keep this variant
        bool keepVariant = false;
        for (double af : afValues) {
            if (af >= minAF && af <= maxAF) {
                keepVariant = true;
                break;
            }
        }

        if (keepVariant) {
            out << line << "\n";
        }
    }
}

//
// Typical main():
//
int main(int argc, char* argv[]) {
    VCFXAfSubsetter afSubsetter;
    return afSubsetter.run(argc, argv);
}


./VCFX_af_subsetter/VCFX_af_subsetter.h
#ifndef VCFX_AF_SUBSETTER_H
#define VCFX_AF_SUBSETTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXAfSubsetter: Header file for Alternate Allele Frequency Subsetter Tool
class VCFXAfSubsetter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Subsets VCF input based on alternate allele frequency range
    void subsetByAlleleFrequency(std::istream& in, std::ostream& out, double minAF, double maxAF);

    // Parses the AF values from the INFO field (handles multi-allelic AF as comma-delimited)
    bool parseAF(const std::string& infoField, std::vector<double>& afValues);
};

#endif // VCFX_AF_SUBSETTER_H


./VCFX_alignment_checker/VCFX_alignment_checker.cpp
#include "VCFX_alignment_checker.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>

int VCFXAlignmentChecker::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    static struct option long_options[] = {
        {"help",                 no_argument,       0, 'h'},
        {"alignment-discrepancy", no_argument,      0, 'a'},
        {0,                      0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "ha", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'a':
                // Alignment-discrepancy mode (no extra actions needed here)
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    // We expect two positional arguments after the options:
    //   1) VCF file  2) Reference FASTA file
    if ((argc - optind) < 2) {
        std::cerr << "Error: Missing required arguments: <vcf_file> <reference.fasta>\n";
        displayHelp();
        return 1;
    }

    // Get file paths
    std::string vcfFile = argv[optind];
    std::string refFile = argv[optind + 1];

    // Open VCF input
    std::ifstream vcfStream(vcfFile);
    if (!vcfStream.is_open()) {
        std::cerr << "Error: Unable to open VCF file: " << vcfFile << "\n";
        return 1;
    }

    // Open reference genome file
    std::ifstream refStream(refFile);
    if (!refStream.is_open()) {
        std::cerr << "Error: Unable to open reference genome file: " << refFile << "\n";
        return 1;
    }

    // Load reference genome into memory
    if (!loadReferenceGenome(refStream)) {
        std::cerr << "Error: Failed to load reference genome.\n";
        return 1;
    }

    // Check discrepancies (results to stdout)
    checkDiscrepancies(vcfStream, std::cout);

    return 0;
}

void VCFXAlignmentChecker::displayHelp() {
    std::cout << "VCFX_alignment_checker: Identify discrepancies between VCF variants and a reference genome.\n\n"
              << "Usage:\n"
              << "  VCFX_alignment_checker --alignment-discrepancy <vcf_file> <reference.fasta>\n\n"
              << "Options:\n"
              << "  -h, --help                   Display this help message and exit\n"
              << "  -a, --alignment-discrepancy  Identify alignment discrepancies\n\n"
              << "Example:\n"
              << "  VCFX_alignment_checker --alignment-discrepancy input.vcf reference.fasta > discrepancies.txt\n";
}

bool VCFXAlignmentChecker::loadReferenceGenome(std::istream& in) {
    std::string line;
    std::string currentChrom;
    std::string seq;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            // If we already had a chromosome loaded, store its sequence
            if (!currentChrom.empty()) {
                referenceGenome[normalizeChromosome(currentChrom)] = seq;
            }
            // Start a new chromosome
            seq.clear();
            // Grab chromosome name (up to first space)
            size_t pos = line.find(' ');
            if (pos != std::string::npos) {
                currentChrom = line.substr(1, pos - 1);
            } else {
                currentChrom = line.substr(1);
            }
        } else {
            // Append this line to the sequence (uppercase)
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            seq += line;
        }
    }

    // Store the last chromosome read
    if (!currentChrom.empty()) {
        referenceGenome[normalizeChromosome(currentChrom)] = seq;
    }

    return true;
}

std::string VCFXAlignmentChecker::normalizeChromosome(const std::string& chrom) {
    // NOTE: This logic may cause mismatches if your reference is named "1" but your VCF says "chr1".
    // You may want to adjust this to match your actual naming conventions.
    std::string norm = chrom;
    if (norm.find("chr") != 0 && 
        !(norm == "X" || norm == "Y" || norm == "MT" ||
          std::all_of(norm.begin(), norm.end(), ::isdigit))) 
    {
        norm = "chr" + norm;
    }
    return norm;
}

std::string VCFXAlignmentChecker::getReferenceBases(const std::string& chrom, int pos, int length) {
    auto it = referenceGenome.find(normalizeChromosome(chrom));
    if (it == referenceGenome.end()) {
        return "";
    }

    const std::string& seq = it->second;
    // Convert VCF 1-based 'pos' to a 0-based index into the string
    size_t startIndex = static_cast<size_t>(pos - 1);
    if (pos < 1 || (startIndex + length) > seq.size()) {
        return "";
    }
    return seq.substr(startIndex, length);
}

void VCFXAlignmentChecker::checkDiscrepancies(std::istream& vcfIn, std::ostream& out) {
    // We do not need an extra reference stream parameter; the reference is already loaded in memory.

    std::string line;
    bool headerParsed = false;
    int chrIndex = -1, posIndex = -1, refIndex = -1, altIndex = -1;

    // Print a header for the discrepancies table
    out << "CHROM\tPOS\tID\tREF\tALT\tDiscrepancy_Type\tReference_Value\tVCF_Value\n";

    while (std::getline(vcfIn, line)) {
        if (line.empty()) {
            continue;
        }

        // Header lines
        if (line[0] == '#') {
            // If it's the #CHROM line, parse the column indices
            if (line.rfind("#CHROM", 0) == 0) {
                std::stringstream ss(line);
                std::string field;
                std::vector<std::string> headers;
                while (std::getline(ss, field, '\t')) {
                    headers.push_back(field);
                }
                for (size_t i = 0; i < headers.size(); ++i) {
                    if (headers[i] == "CHROM") chrIndex = static_cast<int>(i);
                    else if (headers[i] == "POS")   posIndex = static_cast<int>(i);
                    else if (headers[i] == "REF")   refIndex = static_cast<int>(i);
                    else if (headers[i] == "ALT")   altIndex = static_cast<int>(i);
                }
                if (chrIndex == -1 || posIndex == -1 || refIndex == -1 || altIndex == -1) {
                    std::cerr << "Error: VCF header does not contain required CHROM, POS, REF, ALT fields.\n";
                    return;
                }
                headerParsed = true;
            }
            continue;
        }

        // We should have a valid header by now
        if (!headerParsed) {
            std::cerr << "Error: VCF #CHROM header line not found before data lines.\n";
            return;
        }

        // Split line into fields
        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string field;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        // Basic sanity check
        if (fields.size() < static_cast<size_t>(altIndex + 1)) {
            std::cerr << "Warning: Skipping invalid VCF line (insufficient fields): " << line << "\n";
            continue;
        }

        std::string chrom = fields[chrIndex];
        int posVal = 0;
        try {
            posVal = std::stoi(fields[posIndex]);
        } catch (...) {
            std::cerr << "Warning: Invalid POS value. Skipping line: " << line << "\n";
            continue;
        }

        std::string ref = fields[refIndex];
        std::string alt = fields[altIndex];
        std::string id  = (fields.size() > 2) ? fields[2] : ".";

        // Handle multi-allelic ALTs
        std::vector<std::string> alts;
        {
            std::stringstream altSS(alt);
            while (std::getline(altSS, field, ',')) {
                alts.push_back(field);
            }
        }

        // Check each ALT allele
        for (const auto& allele : alts) {
            // If both REF and ALT are single bases, treat as SNP
            if (ref.size() == 1 && allele.size() == 1) {
                std::string ref_base = getReferenceBases(chrom, posVal, 1);
                if (ref_base.empty()) {
                    std::cerr << "Warning: Reference base not found for " << chrom << ":" << posVal << "\n";
                    continue;
                }
                // Compare REF in VCF vs reference genome
                if (ref != ref_base) {
                    out << chrom << "\t" << posVal << "\t" << id << "\t" << ref
                        << "\t" << allele << "\t" << "REF_MISMATCH"
                        << "\t" << ref_base << "\t" << ref << "\n";
                }
                // Compare ALT in VCF vs reference genome's same position
                // (Often for a standard SNP, the reference base is the only thing in the FASTA.)
                // This is somewhat conceptual: we're checking if the ALT base is the same as reference at that position.
                std::string alt_base = ref_base; // The reference at that position
                if (allele != alt_base) {
                    out << chrom << "\t" << posVal << "\t" << id << "\t" << ref
                        << "\t" << allele << "\t" << "ALT_MISMATCH"
                        << "\t" << alt_base << "\t" << allele << "\n";
                }
            } else {
                // Indel or complex
                size_t ref_len = ref.size();
                size_t alt_len = allele.size();
                // Compare as many bases as the shorter string has
                size_t len = std::min(ref_len, alt_len);

                std::string ref_seq = getReferenceBases(chrom, posVal, static_cast<int>(len));
                if (ref_seq.empty()) {
                    std::cerr << "Warning: Reference sequence not found for " << chrom << ":" << posVal << "\n";
                    continue;
                }

                std::string vcf_ref = ref.substr(0, len);
                std::string vcf_alt = allele.substr(0, len);

                if (vcf_ref != ref_seq) {
                    out << chrom << "\t" << posVal << "\t" << id << "\t" << ref
                        << "\t" << allele << "\t" << "REF_DISCREPANCY"
                        << "\t" << ref_seq << "\t" << vcf_ref << "\n";
                }
                if (vcf_alt != ref_seq) {
                    out << chrom << "\t" << posVal << "\t" << id << "\t" << ref
                        << "\t" << allele << "\t" << "ALT_DISCREPANCY"
                        << "\t" << ref_seq << "\t" << vcf_alt << "\n";
                }
            }
        }
    }
}

// Typical main(), linking to run()
int main(int argc, char* argv[]) {
    VCFXAlignmentChecker alignmentChecker;
    return alignmentChecker.run(argc, argv);
}


./VCFX_alignment_checker/VCFX_alignment_checker.h
#ifndef VCFX_ALIGNMENT_CHECKER_H
#define VCFX_ALIGNMENT_CHECKER_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// VCFXAlignmentChecker: Header file for Reference Alignment Discrepancy Finder Tool
class VCFXAlignmentChecker {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Loads the reference genome from a FASTA file
    bool loadReferenceGenome(std::istream& in);

    // Checks discrepancies between VCF variants and the in-memory reference genome
    void checkDiscrepancies(std::istream& vcfIn, std::ostream& out);

    // Retrieves the reference base(s) from the reference genome at a specific position
    std::string getReferenceBases(const std::string& chrom, int pos, int length = 1);

    // Stores the reference genome sequences, keyed by normalized chromosome name
    std::unordered_map<std::string, std::string> referenceGenome;

    // Helper function to convert chromosome names to a consistent format
    std::string normalizeChromosome(const std::string& chrom);
};

#endif // VCFX_ALIGNMENT_CHECKER_H


./VCFX_allele_balance_calc/VCFX_allele_balance_calc.cpp
// VCFX_allele_balance_calc.cpp

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <cstdlib>

// ---------------------------------------------------------------
// Header / Declarations
// ---------------------------------------------------------------
struct AlleleBalanceArguments {
    std::vector<std::string> samples;
};

// Forward declarations
void printHelp();
bool parseArguments(int argc, char* argv[], AlleleBalanceArguments& args);
bool calculateAlleleBalance(std::istream& in, std::ostream& out, const AlleleBalanceArguments& args);
double computeAlleleBalance(const std::string& genotype);

// ---------------------------------------------------------------
// Utility Functions
// ---------------------------------------------------------------
static std::vector<std::string> splitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string temp;
    while (std::getline(ss, temp, delimiter)) {
        tokens.push_back(temp);
    }
    return tokens;
}

// Print help message
void printHelp() {
    std::cout 
        << "VCFX_allele_balance_calc\n"
        << "Usage: VCFX_allele_balance_calc [OPTIONS] < input.vcf > allele_balance.tsv\n\n"
        << "Options:\n"
        << "  --samples, -s \"Sample1 Sample2\"   Specify the sample names to calculate allele balance for.\n"
        << "  --help, -h                        Display this help message and exit.\n\n"
        << "Description:\n"
        << "  Calculates the allele balance (ratio of reference to alternate alleles) for each sample.\n"
        << "  Allele balance is computed as (#RefAlleles / #AltAlleles), using the genotype field.\n"
        << "  This simple logic treats all non-zero alleles as 'alt' and 0 as 'ref',\n"
        << "  so multi-allelic sites are lumped into an overall alt count.\n\n"
        << "Examples:\n"
        << "  1) Calculate allele balance for SampleA and SampleB:\n"
        << "     ./VCFX_allele_balance_calc --samples \"SampleA SampleB\" < input.vcf > allele_balance.tsv\n\n"
        << "  2) Calculate allele balance for all samples:\n"
        << "     ./VCFX_allele_balance_calc < input.vcf > allele_balance_all.tsv\n\n";
}

// Parse command-line arguments
bool parseArguments(int argc, char* argv[], AlleleBalanceArguments& args) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--samples" || arg == "-s") && i + 1 < argc) {
            // Collect samples from next argument
            std::string samplesStr = argv[++i];
            args.samples = splitString(samplesStr, ' ');
            // Trim whitespace for each sample
            for (auto &sample : args.samples) {
                // left trim
                sample.erase(0, sample.find_first_not_of(" \t\n\r\f\v"));
                // right trim
                sample.erase(sample.find_last_not_of(" \t\n\r\f\v") + 1);
            }
        }
        else if (arg == "--help" || arg == "-h") {
            printHelp();
            std::exit(0);
        }
        else {
            std::cerr << "Warning: Unrecognized argument '" << arg << "'.\n";
        }
    }
    return true;
}

// Compute allele balance
// If genotype is "0/0", result= ref_count / alt_count => 2 / 0 => 0.0
// If genotype is "0/1", result= 1 / 1 => 1.0
// If genotype is "1/2", result= 0 / 2 => 0.0, etc.
double computeAlleleBalance(const std::string& genotype) {
    if (genotype.empty() || genotype == "." || genotype == "./." || genotype == ".|.") {
        return -1.0; // missing genotype
    }

    // Split by '/' or '|'. For a simple approach, we can replace '|' with '/' and split on '/':
    std::string g = genotype;
    for (auto &ch : g) {
        if (ch == '|') ch = '/';
    }
    auto alleles = splitString(g, '/');

    int ref_count = 0;
    int alt_count = 0;
    for (const auto& allele : alleles) {
        if (allele == "0") {
            ++ref_count;
        } else if (allele != "." && !allele.empty()) {
            ++alt_count;
        }
    }

    // If alt_count is 0 but we have some reference calls, ratio is 0.0
    if (alt_count == 0 && ref_count > 0) return 0.0;
    // If no meaningful data, return -1
    if (ref_count + alt_count == 0) return -1.0;

    return static_cast<double>(ref_count) / alt_count;
}

// Main function to read VCF, parse genotypes, and calculate allele balance
bool calculateAlleleBalance(std::istream& in, std::ostream& out, const AlleleBalanceArguments& args) {
    std::string line;
    bool headerFound = false;
    std::vector<std::string> headerFields;
    std::vector<int> sampleIndices; // columns (9-based) that we care about
    std::unordered_map<std::string, int> sampleMap;

    // Print exactly one TSV header
    out << "CHROM\tPOS\tID\tREF\tALT\tSample\tAllele_Balance\n";

    while (std::getline(in, line)) {
        // Skip empty lines
        if (line.empty()) {
            continue;
        }

        // Skip and do not print lines that start with '#'
        //   BUT parse #CHROM line to get sample columns
        if (line[0] == '#') {
            // If this is the #CHROM line, parse columns
            if (line.rfind("#CHROM", 0) == 0) {
                headerFields = splitString(line, '\t');
                // Build sample -> column index
                for (size_t i = 9; i < headerFields.size(); ++i) {
                    sampleMap[headerFields[i]] = static_cast<int>(i);
                }
                // If user specified samples, pick them out; otherwise take all
                if (!args.samples.empty()) {
                    for (const auto& sample : args.samples) {
                        auto it = sampleMap.find(sample);
                        if (it == sampleMap.end()) {
                            std::cerr << "Error: Sample '" << sample << "' not found in VCF header.\n";
                            return false;
                        }
                        sampleIndices.push_back(it->second);
                    }
                } else {
                    // all samples
                    for (size_t i = 9; i < headerFields.size(); ++i) {
                        sampleIndices.push_back(static_cast<int>(i));
                    }
                }
                headerFound = true;
            }
            continue; // do not output these lines
        }

        // We need a #CHROM header before data
        if (!headerFound) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        // Split the data line
        auto fields = splitString(line, '\t');
        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping invalid VCF line with fewer than 9 fields.\n";
            continue;
        }

        // Basic VCF columns
        const std::string &chrom = fields[0];
        const std::string &pos   = fields[1];
        const std::string &id    = fields[2];
        const std::string &ref   = fields[3];
        const std::string &alt   = fields[4];
        // The genotype info starts at index 9

        // Calculate AB for each requested sample
        for (auto idx : sampleIndices) {
            if (static_cast<size_t>(idx) >= fields.size()) {
                std::cerr << "Warning: sample index " << idx << " out of range.\n";
                continue;
            }
            const std::string &genotypeStr = fields[idx];
            // genotype is typically the first colon-delimited field
            // e.g. genotypeStr = "0/1:..." 
            // but we only need the actual genotype for computing ratio
            // split on ':' 
            auto subFields = splitString(genotypeStr, ':');
            if (subFields.empty()) {
                continue;
            }
            std::string genotype = subFields[0]; // e.g. "0/1"

            double ab = computeAlleleBalance(genotype);
            std::string abStr = (ab < 0.0) ? "NA" : std::to_string(ab);

            // The sample name is from the VCF header
            std::string sampleName = (idx < static_cast<int>(headerFields.size()))
                                     ? headerFields[idx]
                                     : ("Sample_" + std::to_string(idx));

            out << chrom << "\t"
                << pos   << "\t"
                << id    << "\t"
                << ref   << "\t"
                << alt   << "\t"
                << sampleName << "\t"
                << abStr << "\n";
        }
    }
    return true;
}

// ---------------------------------------------------------------
// main()
// ---------------------------------------------------------------
int main(int argc, char* argv[]) {
    AlleleBalanceArguments args;
    parseArguments(argc, argv, args);

    if (!args.samples.empty()) {
        std::cerr << "Info: Calculating allele balance for these samples: ";
        for (auto &s : args.samples) {
            std::cerr << s << " ";
        }
        std::cerr << "\n";
    } else {
        std::cerr << "Info: Calculating allele balance for ALL samples.\n";
    }

    bool success = calculateAlleleBalance(std::cin, std::cout, args);
    return success ? 0 : 1;
}


./VCFX_allele_balance_calc/VCFX_allele_balance_calc.h
#ifndef VCFX_ALLELE_BALANCE_CALC_H
#define VCFX_ALLELE_BALANCE_CALC_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Structure to hold command-line arguments
struct AlleleBalanceArguments {
    std::vector<std::string> samples;
};

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], AlleleBalanceArguments& args);

// Function to split a string by a delimiter
std::vector<std::string> splitString(const std::string& str, char delimiter);

// Function to process VCF and calculate allele balance
bool calculateAlleleBalance(std::istream& in, std::ostream& out, const AlleleBalanceArguments& args);

// Function to extract genotype from genotype string
std::string extractGenotype(const std::string& genotype_str);

// Function to calculate allele balance ratio
double computeAlleleBalance(const std::string& genotype);

#endif // VCFX_ALLELE_BALANCE_CALC_H


./VCFX_allele_balance_filter/VCFX_allele_balance_filter.cpp
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <getopt.h>

// ------------------------------------------------------
// Class Declaration
// ------------------------------------------------------
class VCFXAlleleBalanceFilter {
public:
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input based on allele balance threshold
    void filterByAlleleBalance(std::istream& in, std::ostream& out, double threshold);

    // Parses the genotype to calculate allele balance
    double calculateAlleleBalance(const std::string& genotype);
};

// ------------------------------------------------------
// Implementation
// ------------------------------------------------------

// Entry point for the tool
int VCFXAlleleBalanceFilter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    double threshold = -1.0; // invalid until set

    static struct option long_options[] = {
        {"help",                 no_argument,       0, 'h'},
        {"filter-allele-balance", required_argument, 0, 'f'},
        {0,                      0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hf:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'f':
                try {
                    threshold = std::stod(optarg);
                } catch (const std::invalid_argument&) {
                    std::cerr << "Error: Invalid threshold value.\n";
                    displayHelp();
                    return 1;
                }
                break;
            default:
                showHelp = true;
        }
    }

    // Validate threshold and possibly show help
    if (showHelp || threshold < 0.0 || threshold > 1.0) {
        displayHelp();
        return (showHelp ? 0 : 1);
    }

    // Perform allele balance filtering from stdin to stdout
    filterByAlleleBalance(std::cin, std::cout, threshold);

    return 0;
}

void VCFXAlleleBalanceFilter::displayHelp() {
    std::cout 
        << "VCFX_allele_balance_filter: Filter VCF variants based on allele balance ratios.\n\n"
        << "Usage:\n"
        << "  VCFX_allele_balance_filter --filter-allele-balance <THRESHOLD> [options]\n\n"
        << "Options:\n"
        << "  -h, --help                       Display this help message and exit\n"
        << "  -f, --filter-allele-balance VAL  Specify the allele balance threshold (0.0 - 1.0)\n\n"
        << "Example:\n"
        << "  VCFX_allele_balance_filter --filter-allele-balance 0.3 < input.vcf > filtered.vcf\n\n"
        << "Note:\n"
        << "  This filter lumps all non-'0' alleles (1,2,3,...) as ALT when calculating the ratio.\n"
        << "  If any sample's allele balance is < THRESHOLD, the entire variant line is skipped.\n";
}

// The core filter function
void VCFXAlleleBalanceFilter::filterByAlleleBalance(std::istream& in, std::ostream& out, double threshold) {
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        // Print header lines unchanged
        if (line[0] == '#') {
            out << line << "\n";
            continue;
        }

        // Parse the fixed VCF fields up to FORMAT
        // VCF line: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [samples...]
        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info, format;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        // The rest are genotype columns
        std::vector<std::string> genotypes;
        std::string genotype;
        while (ss >> genotype) {
            genotypes.push_back(genotype);
        }

        // "all or nothing" filter: if ANY genotype fails threshold => skip
        bool pass = true;
        for (const auto& gtField : genotypes) {
            double ab = calculateAlleleBalance(gtField);
            if (ab < threshold) {
                pass = false;
                break;
            }
        }

        if (pass) {
            out << line << "\n";
        }
    }
}

// Calculates allele balance as refCount / (refCount + altCount)
// *all* non-zero numeric alleles are counted as alt
double VCFXAlleleBalanceFilter::calculateAlleleBalance(const std::string& genotype) {
    // Genotype might be "0/1:...". We want the portion before the first ':'
    size_t colonPos = genotype.find(':');
    std::string gt = (colonPos == std::string::npos) ? genotype : genotype.substr(0, colonPos);

    // Replace '|' with '/' for simpler splitting
    for (auto &ch : gt) {
        if (ch == '|') ch = '/';
    }

    // Split on '/'
    std::stringstream ss(gt);
    std::string token;
    std::vector<std::string> alleles;
    while (std::getline(ss, token, '/')) {
        alleles.push_back(token);
    }

    int ref_count = 0;
    int alt_count = 0;

    // Count '0' as ref, any other number as alt
    for (auto & allele : alleles) {
        // If allele is empty or ".", skip it
        if (allele.empty() || allele == ".") {
            continue;
        }
        // If we can parse it as an integer:
        //   0 => ref, 1 or 2 or 3 => alt
        // If it's non-numeric, we skip it as missing.
        bool numeric = true;
        for (char c : allele) {
            if (!isdigit(c)) {
                numeric = false;
                break;
            }
        }
        if (!numeric) {
            // Non-numeric => skip?
            continue; 
        }
        // numeric => check if it's "0" or not
        if (allele == "0") {
            ref_count++;
        } else {
            alt_count++;
        }
    }

    int total = ref_count + alt_count;
    if (total == 0) {
        // No valid calls => ratio 0.0 or treat as "no data"
        return 0.0;
    }
    return static_cast<double>(ref_count) / total;
}

// ------------------------------------------------------
// main() linking to class
// ------------------------------------------------------
int main(int argc, char* argv[]) {
    VCFXAlleleBalanceFilter alleleBalanceFilter;
    return alleleBalanceFilter.run(argc, argv);
}


./VCFX_allele_balance_filter/VCFX_allele_balance_filter.h
#ifndef VCFX_ALLELE_BALANCE_FILTER_H
#define VCFX_ALLELE_BALANCE_FILTER_H

#include <iostream>
#include <string>
#include <vector>
#include <sstream>

// VCFXAlleleBalanceFilter: Header file for Allele Balance Filter Tool
class VCFXAlleleBalanceFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input based on allele balance threshold
    void filterByAlleleBalance(std::istream& in, std::ostream& out, double threshold);

    // Parses the genotype to calculate allele balance
    double calculateAlleleBalance(const std::string& genotype);
};

#endif // VCFX_ALLELE_BALANCE_FILTER_H


./VCFX_allele_counter/VCFX_allele_counter.cpp
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdlib>

// ---------------------------------------------------------------------
// Structures and Declarations
// ---------------------------------------------------------------------
struct AlleleCounterArguments {
    std::vector<std::string> samples;
};

static void printHelp();
static bool parseArguments(int argc, char* argv[], AlleleCounterArguments& args);
static std::vector<std::string> splitString(const std::string& str, char delimiter);
static bool countAlleles(std::istream& in, std::ostream& out, const AlleleCounterArguments& args);

// ---------------------------------------------------------------------
// printHelp
// ---------------------------------------------------------------------
static void printHelp() {
    std::cout
        << "VCFX_allele_counter\n"
        << "Usage: VCFX_allele_counter [OPTIONS] < input.vcf > allele_counts.tsv\n\n"
        << "Options:\n"
        << "  --samples, -s \"Sample1 Sample2\"   Specify the sample names to include.\n"
        << "                                     If not specified, all samples are processed.\n"
        << "  --help, -h                        Display this help message and exit.\n\n"
        << "Description:\n"
        << "  Reads a VCF from stdin and outputs a TSV file with the columns:\n"
        << "    CHROM  POS  ID  REF  ALT  Sample  Ref_Count  Alt_Count\n\n"
        << "  Each sample for each variant is listed. Alleles are determined by\n"
        << "  genotype strings (GT). This code treats any numeric allele that is\n"
        << "  not '0' as ALT, e.g. '1' or '2' or '3' => alt.\n\n"
        << "Examples:\n"
        << "  1) Count alleles for SampleA and SampleB:\n"
        << "     VCFX_allele_counter --samples \"SampleA SampleB\" < input.vcf > allele_counts.tsv\n\n"
        << "  2) Count alleles for all samples:\n"
        << "     VCFX_allele_counter < input.vcf > allele_counts_all.tsv\n";
}

// ---------------------------------------------------------------------
// parseArguments
// ---------------------------------------------------------------------
static bool parseArguments(int argc, char* argv[], AlleleCounterArguments& args) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--samples" || arg == "-s") && i + 1 < argc) {
            std::string samplesStr = argv[++i];
            args.samples = splitString(samplesStr, ' ');
            // Trim whitespace from each sample name
            for (auto& sample : args.samples) {
                // left trim
                sample.erase(0, sample.find_first_not_of(" \t\n\r\f\v"));
                // right trim
                sample.erase(sample.find_last_not_of(" \t\n\r\f\v") + 1);
            }
        }
        else if (arg == "--help" || arg == "-h") {
            printHelp();
            std::exit(0);
        }
        else {
            std::cerr << "Warning: Unrecognized argument '" << arg << "'.\n";
        }
    }
    return true;
}

// ---------------------------------------------------------------------
// splitString
// ---------------------------------------------------------------------
static std::vector<std::string> splitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string temp;
    while (std::getline(ss, temp, delimiter)) {
        tokens.push_back(temp);
    }
    return tokens;
}

// ---------------------------------------------------------------------
// countAlleles
// ---------------------------------------------------------------------
static bool countAlleles(std::istream& in, std::ostream& out, const AlleleCounterArguments& args) {
    std::string line;
    std::vector<std::string> headerFields;
    bool foundChromHeader = false;

    // We will store the columns (indexes) for the samples we care about
    std::vector<int> sampleIndices;

    // Print one TSV header (not the original VCF header lines)
    out << "CHROM\tPOS\tID\tREF\tALT\tSample\tRef_Count\tAlt_Count\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        // Skip lines that start with '#' (VCF headers), but parse #CHROM line
        if (line[0] == '#') {
            // If this is the #CHROM line, parse to find sample columns
            if (line.rfind("#CHROM", 0) == 0) {
                headerFields = splitString(line, '\t');
                if (headerFields.size() < 9) {
                    std::cerr << "Error: #CHROM line has fewer than 9 columns.\n";
                    return false;
                }
                // Build sample->index map
                std::unordered_map<std::string, int> sampleMap;
                for (size_t i = 9; i < headerFields.size(); ++i) {
                    sampleMap[headerFields[i]] = static_cast<int>(i);
                }

                // Decide which samples to process
                if (!args.samples.empty()) {
                    // Use only user-specified samples
                    for (const auto& s : args.samples) {
                        auto it = sampleMap.find(s);
                        if (it == sampleMap.end()) {
                            std::cerr << "Error: Sample '" << s << "' not found in VCF header.\n";
                            return false;
                        }
                        sampleIndices.push_back(it->second);
                    }
                } else {
                    // If no samples are specified, process all
                    for (size_t i = 9; i < headerFields.size(); ++i) {
                        sampleIndices.push_back((int)i);
                    }
                }
                foundChromHeader = true;
            }
            // We do NOT print the raw VCF headers to output
            continue;
        }

        // We expect a #CHROM line before data lines
        if (!foundChromHeader) {
            std::cerr << "Error: VCF #CHROM header not found before records.\n";
            return false;
        }

        // Split VCF data line
        auto fields = splitString(line, '\t');
        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping invalid VCF line with fewer than 9 fields:\n" 
                      << line << "\n";
            continue;
        }

        // Basic columns
        const std::string &chrom = fields[0];
        const std::string &pos   = fields[1];
        const std::string &id    = fields[2];
        const std::string &ref   = fields[3];
        const std::string &alt   = fields[4];

        // For each chosen sample, parse genotype
        for (int sIndex : sampleIndices) {
            if (sIndex >= (int)fields.size()) {
                std::cerr << "Warning: Sample index " << sIndex << " out of range for line:\n"
                          << line << "\n";
                continue;
            }
            const std::string &sampleField = fields[sIndex];

            // genotype is the portion before the first ':'
            size_t colonPos = sampleField.find(':');
            std::string gt = (colonPos == std::string::npos)
                             ? sampleField
                             : sampleField.substr(0, colonPos);

            // Replace '|' with '/' for a uniform split
            for (auto &c : gt) {
                if (c == '|') c = '/';
            }
            auto alleles = splitString(gt, '/');
            if (alleles.empty()) {
                // e.g. '.' or blank
                continue;
            }

            int refCount = 0;
            int altCount = 0;
            for (auto & allele : alleles) {
                // Skip empty or '.' allele
                if (allele.empty() || allele == ".") {
                    continue;
                }
                // If numeric & "0" => ref, else alt
                bool numeric = true;
                for (char c : allele) {
                    if (!isdigit(c)) { numeric = false; break; }
                }
                if (!numeric) {
                    // e.g. unknown character => skip
                    continue;
                }
                // If allele == "0" => ref, else alt
                if (allele == "0") {
                    refCount++;
                } else {
                    altCount++;
                }
            }

            // The sample's name
            std::string sampleName = (sIndex < (int)headerFields.size())
                                     ? headerFields[sIndex]
                                     : ("Sample_"+std::to_string(sIndex));

            // Output as TSV row
            out << chrom << "\t"
                << pos   << "\t"
                << id    << "\t"
                << ref   << "\t"
                << alt   << "\t"
                << sampleName << "\t"
                << refCount   << "\t"
                << altCount   << "\n";
        }
    }
    return true;
}

// ---------------------------------------------------------------------
// main()
// ---------------------------------------------------------------------
int main(int argc, char* argv[]) {
    AlleleCounterArguments args;
    parseArguments(argc, argv, args);

    if (!args.samples.empty()) {
        std::cerr << "Info: Counting alleles for these samples:\n";
        for (auto &s : args.samples) {
            std::cerr << "  " << s << "\n";
        }
    } else {
        std::cerr << "Info: Counting alleles for ALL samples.\n";
    }

    bool success = countAlleles(std::cin, std::cout, args);
    return success ? 0 : 1;
}


./VCFX_allele_counter/VCFX_allele_counter.h
#ifndef VCFX_ALLELE_COUNTER_H
#define VCFX_ALLELE_COUNTER_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Structure to hold command-line arguments
struct AlleleCounterArguments {
    std::vector<std::string> samples;
};

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], AlleleCounterArguments& args);

// Function to split a string by a delimiter
std::vector<std::string> splitString(const std::string& str, char delimiter);

// Function to process VCF and count alleles
bool countAlleles(std::istream& in, std::ostream& out, const AlleleCounterArguments& args);

// Function to extract genotype from genotype string
std::string extractGenotype(const std::string& genotype_str);

#endif // VCFX_ALLELE_COUNTER_H


./VCFX_allele_freq_calc/VCFX_allele_freq_calc.cpp
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cctype>

// ---------------------------------------------------------
// Print help message
// ---------------------------------------------------------
static void printHelp() {
    std::cout 
        << "VCFX_allele_freq_calc\n"
        << "Usage: VCFX_allele_freq_calc [OPTIONS]\n\n"
        << "Options:\n"
        << "  --help, -h   Display this help message and exit.\n\n"
        << "Description:\n"
        << "  Reads a VCF from stdin and outputs a TSV file:\n"
        << "    CHROM  POS  ID  REF  ALT  Allele_Frequency\n\n"
        << "  Allele frequency is computed as (#ALT alleles / total #alleles),\n"
        << "  counting any non-zero numeric allele (1,2,3,...) as ALT.\n\n"
        << "Example:\n"
        << "  ./VCFX_allele_freq_calc < input.vcf > allele_frequencies.tsv\n";
}

// ---------------------------------------------------------
// Utility: split a string by a delimiter
// ---------------------------------------------------------
static std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// ---------------------------------------------------------
// parseGenotype: counts how many are alt vs total among numeric alleles
// ---------------------------------------------------------
static void parseGenotype(const std::string& genotype, int &alt_count, int &total_count) {
    // e.g. genotype might be "0/1" or "2|3", etc.
    // First, unify the separators by replacing '|' with '/'
    std::string gt = genotype;
    for (char &c : gt) {
        if (c == '|') c = '/';
    }

    // Split on '/'
    auto alleles = split(gt, '/');

    for (const auto & allele : alleles) {
        if (allele.empty() || allele == ".") {
            // Skip missing
            continue;
        }
        // Check if it's numeric
        bool numeric = true;
        for (char c : allele) {
            if (!std::isdigit(c)) {
                numeric = false;
                break;
            }
        }
        if (!numeric) {
            // Not a numeric allele -> skip
            continue;
        }
        // If it's "0", it's reference; otherwise it's alt
        if (allele == "0") {
            // The reference allele still contributes to the total
            total_count += 1;
        } else {
            // Any non-zero numeric allele => alt
            alt_count += 1;
            total_count += 1;
        }
    }
}

// ---------------------------------------------------------
// Main calculation: read VCF from 'in', write results to 'out'
// ---------------------------------------------------------
static void calculateAlleleFrequency(std::istream& in, std::ostream& out) {
    std::string line;

    // We'll track how many samples are in the VCF, though we only need it
    // to confirm that columns 10..N contain sample data
    int sample_count = 0;
    bool foundChromHeader = false;

    // Print a single TSV header for our results
    out << "CHROM\tPOS\tID\tREF\tALT\tAllele_Frequency\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            // If this is the #CHROM line, figure out how many samples
            if (line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
                auto fields = split(line, '\t');
                if (fields.size() > 9) {
                    sample_count = static_cast<int>(fields.size() - 9);
                }
            }
            // Skip printing VCF headers in our output
            continue;
        }

        // We expect #CHROM line before actual data lines
        if (!foundChromHeader) {
            std::cerr << "Warning: Data line encountered before #CHROM header. Skipping line:\n"
                      << line << "\n";
            continue;
        }

        // Split the data line
        auto fields = split(line, '\t');
        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping invalid VCF line (fewer than 9 fields):\n"
                      << line << "\n";
            continue;
        }

        // Standard VCF columns:
        //  0:CHROM, 1:POS, 2:ID, 3:REF, 4:ALT, 5:QUAL, 6:FILTER, 7:INFO, 8:FORMAT, 9+:samples
        const std::string &chrom = fields[0];
        const std::string &pos   = fields[1];
        const std::string &id    = fields[2];
        const std::string &ref   = fields[3];
        const std::string &alt   = fields[4];
        const std::string &fmt   = fields[8];

        // Look for the "GT" field index in the FORMAT column
        auto formatFields = split(fmt, ':');
        int gtIndex = -1;
        for (size_t i = 0; i < formatFields.size(); ++i) {
            if (formatFields[i] == "GT") {
                gtIndex = static_cast<int>(i);
                break;
            }
        }
        if (gtIndex < 0) {
            // No GT field found => we cannot parse genotypes
            // skip this line
            continue;
        }

        int alt_count = 0;
        int total_count = 0;

        // For each sample column (starting at index 9)
        for (size_t i = 9; i < fields.size(); ++i) {
            // example sample string: "0/1:30,10:99:..." if GT is first
            auto sampleTokens = split(fields[i], ':');
            if (static_cast<size_t>(gtIndex) >= sampleTokens.size()) {
                // This sample has no GT data
                continue;
            }
            // parse its genotype
            parseGenotype(sampleTokens[gtIndex], alt_count, total_count);
        }

        // Compute allele frequency = alt_count / total_count
        double freq = 0.0;
        if (total_count > 0) {
            freq = static_cast<double>(alt_count) / static_cast<double>(total_count);
        }

        // Print the row with a fixed decimal precision
        out << chrom << "\t" << pos << "\t" << id << "\t"
            << ref   << "\t" << alt << "\t"
            << std::fixed << std::setprecision(4) << freq << "\n";
    }
}

// ---------------------------------------------------------
// main()
// ---------------------------------------------------------
int main(int argc, char* argv[]) {
    // Parse arguments for help
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // If no arguments (and no piped input), just show help
    // (This is optional; if the user pipes data, we can proceed.)
    if (argc == 1) {
        // Not strictly required, but user-friendly
        // Check if there's something on stdin, or just show help
        if (std::cin.peek() == EOF) {
            printHelp();
            return 1;
        }
    }

    // Perform allele frequency calculation
    calculateAlleleFrequency(std::cin, std::cout);
    return 0;
}


./VCFX_allele_freq_calc/VCFX_allele_freq_calc.h
#ifndef VCFX_ALLELE_FREQ_CALC_H
#define VCFX_ALLELE_FREQ_CALC_H

#include <iostream>
#include <string>

// Function to display help message
void printHelp();

// Function to perform allele frequency calculation on VCF records
void calculateAlleleFrequency(std::istream& in, std::ostream& out);

#endif // VCFX_ALLELE_FREQ_CALC_H


./VCFX_ancestry_assigner/VCFX_ancestry_assigner.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdlib>
#include <getopt.h>

// ---------------------------------------------------------
// A helper struct to store command-line options if needed
// ---------------------------------------------------------
struct AncestryOptions {
    bool showHelp = false;
    std::string freqFile;
};

// ---------------------------------------------------------
// Class: VCFXAncestryAssigner
// ---------------------------------------------------------
class VCFXAncestryAssigner {
public:
    // High-level entry point
    int run(int argc, char* argv[]);

private:
    // Command-line usage
    void displayHelp();

    // Parse freq file
    bool loadAncestralFrequencies(std::istream& in);

    // Parse one frequency line
    // Format (tab-separated): CHROM  POS  REF  ALT  POP1_FREQ  POP2_FREQ ...
    bool parseFrequencyLine(const std::string& line);

    // Assign ancestry from VCF
    void assignAncestry(std::istream& vcfIn, std::ostream& out);

private:
    // Populations in order
    std::vector<std::string> populations;

    // Frequencies by key: "chrom:pos:ref:alt" => (pop => freq)
    std::unordered_map<std::string, std::unordered_map<std::string, double>> variantFrequencies;
};

// ---------------------------------------------------------
// VCFXAncestryAssigner::run
// ---------------------------------------------------------
int VCFXAncestryAssigner::run(int argc, char* argv[]) {
    // 1. Parse arguments
    int opt;
    bool showHelp = false;
    std::string freqFile;

    static struct option long_options[] = {
        {"help",            no_argument,       0, 'h'},
        {"assign-ancestry", required_argument, 0, 'a'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "ha:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'a':
                freqFile = std::string(optarg);
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || freqFile.empty()) {
        displayHelp();
        return (showHelp ? 0 : 1);
    }

    // 2. Open frequency file
    std::ifstream freqStream(freqFile);
    if (!freqStream.is_open()) {
        std::cerr << "Error: Unable to open frequency file: " << freqFile << "\n";
        return 1;
    }

    // 3. Load frequencies
    if (!loadAncestralFrequencies(freqStream)) {
        std::cerr << "Error: Failed to load ancestral frequencies.\n";
        return 1;
    }

    // 4. Assign ancestry based on VCF (read from stdin, write to stdout)
    assignAncestry(std::cin, std::cout);
    return 0;
}

// ---------------------------------------------------------
// Show usage
// ---------------------------------------------------------
void VCFXAncestryAssigner::displayHelp() {
    std::cout << "VCFX_ancestry_assigner: Assign samples to ancestral populations based on variant frequencies.\n\n"
              << "Usage:\n"
              << "  VCFX_ancestry_assigner --assign-ancestry <freq_file> < input.vcf > ancestry.txt\n\n"
              << "Options:\n"
              << "  -h, --help                 Show this help message and exit\n"
              << "  -a, --assign-ancestry FILE Ancestral frequency file\n\n"
              << "Frequency File Format:\n"
              << "  The first line must be a header like:\n"
              << "    CHROM  POS  REF  ALT  POP1  POP2  ...\n"
              << "  Each subsequent line must have the same columns. For example:\n"
              << "    1   10000   A   C   0.10  0.20\n\n"
              << "Example:\n"
              << "  VCFX_ancestry_assigner --assign-ancestry ancestral_freq.tsv < input.vcf > ancestry_out.txt\n\n";
}

// ---------------------------------------------------------
// parseFrequencyLine
// Expects: CHROM  POS  REF  ALT  pop1Freq  pop2Freq ...
// ---------------------------------------------------------
bool VCFXAncestryAssigner::parseFrequencyLine(const std::string& line) {
    std::stringstream ss(line);
    std::vector<std::string> fields;
    std::string field;
    while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
    }

    // Must have at least CHROM, POS, REF, ALT, plus the populations
    if (fields.size() < 4 + populations.size()) {
        return false;
    }
    const std::string &chrom = fields[0];
    int pos = 0;
    try {
        pos = std::stoi(fields[1]);
    } catch (...) {
        return false;
    }
    const std::string &ref = fields[2];
    const std::string &alt = fields[3];

    // Build freq map for this variant
    std::unordered_map<std::string, double> freqMap;
    for (size_t i = 0; i < populations.size(); ++i) {
        double freqVal = 0.0;
        try {
            freqVal = std::stod(fields[4 + i]);
        } catch (...) {
            freqVal = 0.0; // default if missing
        }
        freqMap[populations[i]] = freqVal;
    }

    // Key = CHROM:POS:REF:ALT
    const std::string key = chrom + ":" + std::to_string(pos) + ":" + ref + ":" + alt;
    variantFrequencies[key] = freqMap;
    return true;
}

// ---------------------------------------------------------
// loadAncestralFrequencies
// First line is header with columns:
//   CHROM  POS  REF  ALT  pop1  pop2 ...
// ---------------------------------------------------------
bool VCFXAncestryAssigner::loadAncestralFrequencies(std::istream& in) {
    std::string line;
    if (!std::getline(in, line)) {
        std::cerr << "Error: Frequency file is empty.\n";
        return false;
    }

    // Parse the header
    {
        std::stringstream ss(line);
        std::vector<std::string> headers;
        std::string h;
        while (std::getline(ss, h, '\t')) {
            headers.push_back(h);
        }
        // We need at least 5 columns: CHROM, POS, REF, ALT, plus 1 pop
        if (headers.size() < 5) {
            std::cerr << "Error: Frequency header must have at least 5 columns.\n";
            return false;
        }
        // The populations start at column index 4
        for (size_t i = 4; i < headers.size(); ++i) {
            populations.push_back(headers[i]);
        }
    }

    // Parse each subsequent line
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (!parseFrequencyLine(line)) {
            std::cerr << "Warning: Skipping invalid frequency line:\n" << line << "\n";
        }
    }
    return true;
}

// ---------------------------------------------------------
// assignAncestry
// Reads VCF from vcfIn, writes "Sample <tab> AssignedPopulation" to out
// ---------------------------------------------------------
void VCFXAncestryAssigner::assignAncestry(std::istream& vcfIn, std::ostream& out) {
    std::string line;
    bool haveHeader = false;
    int chrIndex = -1, posIndex = -1, refIndex = -1, altIndex = -1;
    std::vector<std::string> sampleNames;

    // For each sample: sampleScores[sample][population] => cumulative ancestry score
    std::unordered_map<std::string, std::unordered_map<std::string, double>> sampleScores;

    while (std::getline(vcfIn, line)) {
        if (line.empty()) {
            continue;
        }

        // VCF headers
        if (line[0] == '#') {
            // If #CHROM line, parse for sample columns
            if (line.rfind("#CHROM", 0) == 0) {
                haveHeader = true;
                std::stringstream ss(line);
                std::vector<std::string> headers;
                std::string f;
                while (std::getline(ss, f, '\t')) {
                    headers.push_back(f);
                }
                // Identify columns
                for (size_t i = 0; i < headers.size(); ++i) {
                    if (headers[i] == "CHROM") chrIndex = (int)i;
                    else if (headers[i] == "POS") posIndex = (int)i;
                    else if (headers[i] == "REF") refIndex = (int)i;
                    else if (headers[i] == "ALT") altIndex = (int)i;
                }
                if (chrIndex < 0 || posIndex < 0 || refIndex < 0 || altIndex < 0) {
                    std::cerr << "Error: #CHROM header missing required columns CHROM POS REF ALT.\n";
                    return;
                }
                // Sample columns start at index 9
                for (size_t i = 9; i < headers.size(); ++i) {
                    sampleNames.push_back(headers[i]);
                    // Initialize all sample->pop scores to 0
                    for (auto &pop : populations) {
                        sampleScores[headers[i]][pop] = 0.0;
                    }
                }
            }
            continue; // skip printing these lines
        }

        // Must have #CHROM line before data
        if (!haveHeader) {
            std::cerr << "Error: VCF data encountered before #CHROM line.\n";
            return;
        }

        // Parse a data line
        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string f;
        while (std::getline(ss, f, '\t')) {
            fields.push_back(f);
        }
        if (fields.size() < 9) {
            // Not enough columns
            continue;
        }

        const std::string &chrom = fields[chrIndex];
        int posVal = 0;
        try {
            posVal = std::stoi(fields[posIndex]);
        } catch (...) {
            // invalid POS
            continue;
        }
        const std::string &ref = fields[refIndex];
        const std::string &altField = fields[altIndex];

        // Split ALT by comma for multi-allelic
        std::vector<std::string> alts;
        {
            std::stringstream altSS(altField);
            std::string a;
            while (std::getline(altSS, a, ',')) {
                alts.push_back(a);
            }
        }

        // Identify the GT field index from the FORMAT column
        const std::string &formatCol = fields[8];
        std::vector<std::string> formatParts;
        {
            std::stringstream fmts(formatCol);
            while (std::getline(fmts, f, ':')) {
                formatParts.push_back(f);
            }
        }
        int gtIndex = -1;
        for (size_t i = 0; i < formatParts.size(); ++i) {
            if (formatParts[i] == "GT") {
                gtIndex = (int)i;
                break;
            }
        }
        if (gtIndex < 0) {
            // No genotype data
            continue;
        }

        // For each sample
        for (size_t s = 0; s < sampleNames.size(); ++s) {
            // The sample field is at index (9 + s)
            size_t sampleFieldIndex = 9 + s;
            if (sampleFieldIndex >= fields.size()) {
                continue; // no data
            }
            const std::string &sampleStr = fields[sampleFieldIndex];
            // e.g. "0/1:..." => split by ':'
            std::vector<std::string> sampleParts;
            {
                std::stringstream sss(sampleStr);
                std::string x;
                while (std::getline(sss, x, ':')) {
                    sampleParts.push_back(x);
                }
            }
            if ((int)sampleParts.size() <= gtIndex) {
                // no GT for this sample
                continue;
            }
            // genotype string, e.g. "0/1"
            const std::string &genotype = sampleParts[gtIndex];

            // Replace '|' with '/' for consistency
            std::string gt = genotype;
            for (char &c : gt) {
                if (c == '|') c = '/';
            }

            // Split genotype on '/'
            // e.g. "1/2"
            std::vector<std::string> gtAlleles;
            {
                std::stringstream gts(gt);
                std::string tok;
                while (std::getline(gts, tok, '/')) {
                    gtAlleles.push_back(tok);
                }
            }
            if (gtAlleles.size() < 2) {
                // not a diploid or not parseable
                continue;
            }

            // For each allele in genotype:
            //   0 => REF, 1 => alts[0], 2 => alts[1], etc.
            // If an allele is numeric but >0, we find that alt, look up freq => add to ancestry
            for (size_t alleleIndex = 0; alleleIndex < gtAlleles.size(); ++alleleIndex) {
                const std::string &aStr = gtAlleles[alleleIndex];
                if (aStr.empty() || aStr == ".") {
                    continue; // missing
                }
                bool numeric = true;
                for (char c : aStr) {
                    if (!isdigit(c)) {
                        numeric = false;
                        break;
                    }
                }
                if (!numeric) {
                    continue;
                }
                int alleleInt = std::stoi(aStr);
                if (alleleInt <= 0) {
                    // 0 => ref
                    continue;
                }
                if ((size_t)alleleInt > alts.size()) {
                    // genotype says "3" but we only have 2 alts => skip
                    continue;
                }
                // The alt is alts[alleleInt - 1]
                const std::string &thisAlt = alts[alleleInt - 1];

                // Build key for freq
                std::string key = chrom + ":" + std::to_string(posVal) + ":" + ref + ":" + thisAlt;
                auto it = variantFrequencies.find(key);
                if (it == variantFrequencies.end()) {
                    // We have no freq data for this alt
                    continue;
                }
                // freqMap => pop => freq
                const auto &freqMap = it->second;

                // pick the population with the highest freq
                double maxFreq = -1.0;
                std::string bestPop;
                for (auto &popFreq : freqMap) {
                    if (popFreq.second > maxFreq) {
                        maxFreq = popFreq.second;
                        bestPop = popFreq.first;
                    }
                }
                if (bestPop.empty()) {
                    continue; // no freq data
                }
                // add score for bestPop
                // Each allele we see => + maxFreq
                // (2 alt copies if genotype=1/1 => we do it twice in loop)
                sampleScores[sampleNames[s]][bestPop] += maxFreq;
            }
        }
    }

    // After reading the VCF, pick the best population for each sample
    out << "Sample\tAssigned_Ancestry\n";
    for (auto &sample : sampleScores) {
        const std::string &sampleName = sample.first;
        double bestScore = -1.0;
        std::string bestPop = "NA";
        for (auto &p : sample.second) {
            if (p.second > bestScore) {
                bestScore = p.second;
                bestPop = p.first;
            }
        }
        out << sampleName << "\t" << bestPop << "\n";
    }
}

// ---------------------------------------------------------
// main() - just instantiate and run
// ---------------------------------------------------------
int main(int argc, char* argv[]) {
    VCFXAncestryAssigner assigner;
    return assigner.run(argc, argv);
}


./VCFX_ancestry_assigner/VCFX_ancestry_assigner.h
#ifndef VCFX_ANCESTRY_ASSIGNER_H
#define VCFX_ANCESTRY_ASSIGNER_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// VCFXAncestryAssigner: Header file for Ancestry Assignment Tool
class VCFXAncestryAssigner {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Parses a single line from the frequency file
    bool parseFrequencyLine(const std::string& line, std::string& chrom, int& pos, char& ref, char& alt, std::unordered_map<std::string, double>& freqMap, const std::vector<std::string>& populations);

    // Loads ancestral frequencies from the provided input stream
    bool loadAncestralFrequencies(std::istream& in);

    // Assigns ancestry to samples based on VCF input
    void assignAncestry(std::istream& vcfIn, std::ostream& out);

    // Stores populations
    std::vector<std::string> populations;

    // Stores variant frequencies: Key -> Population -> Frequency
    std::unordered_map<std::string, std::unordered_map<std::string, double>> variantFrequencies;
};

#endif // VCFX_ANCESTRY_ASSIGNER_H


./VCFX_ancestry_inferrer/VCFX_ancestry_inferrer.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdlib>
#include <getopt.h>

// ----------------------------------------------------
// A helper struct for storing freq data by population
// ----------------------------------------------------
struct PopFreqKey {
    // We store CHROM, POS, REF, ALT, POP as strings
    // so we can build a single key. Alternatively, we
    // can keep them separate or use a small struct.
    // We'll build a string key "chrom:pos:ref:alt:pop"
    // for direct indexing in an unordered_map.
    std::string key;
    double frequency;
};

// We’ll keep a direct structure: freqData["chr:pos:ref:alt:POP"] = frequency
typedef std::unordered_map<std::string, double> FrequencyMap;

// ----------------------------------------------------
// Class: VCFXAncestryInferrer
// ----------------------------------------------------
class VCFXAncestryInferrer {
public:
    int run(int argc, char* argv[]);

private:
    // Show usage
    void displayHelp();

    // Load population frequencies from a file that has lines:
    // CHROM  POS  REF  ALT  POPULATION  FREQUENCY
    bool loadPopulationFrequencies(const std::string& freqFilePath);

    // Infer ancestry from VCF
    void inferAncestry(std::istream& vcfInput, std::ostream& output);

private:
    // Frequencies keyed by "chr:pos:ref:alt:pop"
    FrequencyMap freqData;
};

// ----------------------------------------------------
// main() - create the inferrer and run
// ----------------------------------------------------
int main(int argc, char* argv[]) {
    VCFXAncestryInferrer inferrer;
    return inferrer.run(argc, argv);
}

// ----------------------------------------------------
// run() - parse arguments, load freq, run inference
// ----------------------------------------------------
int VCFXAncestryInferrer::run(int argc, char* argv[]) {
    int opt;
    bool showHelp = false;
    std::string freqFilePath;

    static struct option longOpts[] = {
        {"help",       no_argument,       0, 'h'},
        {"frequency",  required_argument, 0, 'f'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "hf:", longOpts, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'f':
                freqFilePath = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || freqFilePath.empty()) {
        displayHelp();
        // Return 0 if help was explicitly requested, otherwise 1
        return (showHelp ? 0 : 1);
    }

    // Load population frequencies
    if (!loadPopulationFrequencies(freqFilePath)) {
        std::cerr << "Error: Failed to load population frequencies from " << freqFilePath << "\n";
        return 1;
    }

    // Read VCF from stdin, write ancestry results to stdout
    inferAncestry(std::cin, std::cout);

    return 0;
}

// ----------------------------------------------------
// displayHelp
// ----------------------------------------------------
void VCFXAncestryInferrer::displayHelp() {
    std::cout 
        << "VCFX_ancestry_inferrer: Infer population ancestry based on allele frequencies.\n\n"
        << "Usage:\n"
        << "  VCFX_ancestry_inferrer --frequency <freq_file> [options]\n\n"
        << "Description:\n"
        << "  Reads a VCF from standard input and outputs a 2-column table:\n"
        << "    Sample  Inferred_Population\n\n"
        << "  The frequency file must have lines of the form:\n"
        << "    CHROM  POS  REF  ALT  POPULATION  FREQUENCY\n"
        << "  (tab-separated). For multi-allelic VCF sites, an ALT allele index 1\n"
        << "  corresponds to the first item in the comma-separated ALT list,\n"
        << "  index 2 => second ALT, etc.\n\n"
        << "Example:\n"
        << "  VCFX_ancestry_inferrer --frequency pop_frequencies.txt < input.vcf > ancestry_results.txt\n";
}

// ----------------------------------------------------
// loadPopulationFrequencies
//   freq file lines: CHROM, POS, REF, ALT, POP, FREQUENCY
// ----------------------------------------------------
bool VCFXAncestryInferrer::loadPopulationFrequencies(const std::string& freqFilePath) {
    std::ifstream freqFile(freqFilePath);
    if (!freqFile.is_open()) {
        std::cerr << "Error: Cannot open frequency file: " << freqFilePath << "\n";
        return false;
    }

    int lineNum = 0;
    std::string line;
    while (std::getline(freqFile, line)) {
        lineNum++;
        if (line.empty()) {
            continue;
        }

        // Parse line
        // e.g. "chr1   12345   A   G   EUR  0.45"
        std::stringstream ss(line);
        std::string chrom, pos, ref, alt, pop, freqStr;
        if (!(ss >> chrom >> pos >> ref >> alt >> pop >> freqStr)) {
            std::cerr << "Warning: Invalid line in frequency file (#" << lineNum << "): " << line << "\n";
            continue;
        }

        double freq = 0.0;
        try {
            freq = std::stod(freqStr);
        } catch (...) {
            std::cerr << "Warning: Invalid frequency value in line #" << lineNum << ": " << line << "\n";
            continue;
        }

        // Build a key: "chr:pos:ref:alt:pop"
        std::string key = chrom + ":" + pos + ":" + ref + ":" + alt + ":" + pop;
        freqData[key] = freq;
    }
    freqFile.close();

    if (freqData.empty()) {
        std::cerr << "Error: No valid population frequencies loaded.\n";
        return false;
    }
    return true;
}

// ----------------------------------------------------
// inferAncestry
//   1) parse the VCF header, get sample names
//   2) for each variant, parse REF + multi-ALT
//   3) for each sample's genotype, add to that sample's
//      population score = sum of freq for each ALT allele
//   4) after all lines, pick population with highest score
// ----------------------------------------------------
void VCFXAncestryInferrer::inferAncestry(std::istream& vcfInput, std::ostream& out) {
    std::string line;
    bool foundChromHeader = false;
    std::vector<std::string> headerFields;
    std::vector<std::string> sampleNames;

    // We'll accumulate ancestry scores: sample -> (pop -> score)
    // Then we pick the highest for each sample
    std::unordered_map<std::string, std::unordered_map<std::string, double>> sampleScores;

    while (std::getline(vcfInput, line)) {
        if (line.empty()) {
            continue;
        }

        // Check if line is a VCF header
        if (line[0] == '#') {
            // The #CHROM line is the last header line, we parse columns
            if (line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
                headerFields.clear();
                {
                    std::stringstream ss(line);
                    std::string tok;
                    while (std::getline(ss, tok, '\t')) {
                        headerFields.push_back(tok);
                    }
                }
                // sample columns start at index 9
                for (size_t c = 9; c < headerFields.size(); ++c) {
                    sampleNames.push_back(headerFields[c]);
                    // initialize each sample's population scores to 0
                    // We'll create the map on-the-fly later, so not strictly needed here
                }
            }
            continue;
        }

        // We must have #CHROM header before data
        if (!foundChromHeader) {
            std::cerr << "Error: Encountered VCF data before #CHROM header.\n";
            return;
        }

        // Parse data line
        // Format (minimally): CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [samples...]
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string f;
            while (std::getline(ss, f, '\t')) {
                fields.push_back(f);
            }
        }

        if (fields.size() < 10) {
            // not enough columns to have samples
            continue;
        }
        // Indices: 0=CHROM, 1=POS, 2=ID, 3=REF, 4=ALT, 5=QUAL, 6=FILTER, 7=INFO, 8=FORMAT, 9+ = samples
        const std::string &chrom = fields[0];
        const std::string &pos   = fields[1];
        // skip ID
        const std::string &ref   = fields[3];
        const std::string &altStr= fields[4];
        const std::string &format= fields[8];

        // Split ALT by comma for multi-allelic
        std::vector<std::string> altAlleles;
        {
            std::stringstream altSS(altStr);
            std::string altTok;
            while (std::getline(altSS, altTok, ',')) {
                altAlleles.push_back(altTok);
            }
        }

        // Find GT index in format
        // e.g. format might be "GT:DP:GQ:..."
        std::vector<std::string> formatParts;
        {
            std::stringstream fmts(format);
            std::string fmtTok;
            while (std::getline(fmts, fmtTok, ':')) {
                formatParts.push_back(fmtTok);
            }
        }
        int gtIndex = -1;
        for (size_t i = 0; i < formatParts.size(); ++i) {
            if (formatParts[i] == "GT") {
                gtIndex = (int)i;
                break;
            }
        }
        if (gtIndex < 0) {
            // no genotype in this line
            continue;
        }

        // For each sample
        for (size_t s = 0; s < sampleNames.size(); ++s) {
            // sample data is fields[9 + s]
            size_t sampleCol = 9 + s;
            if (sampleCol >= fields.size()) {
                continue; 
            }
            const std::string &sampleData = fields[sampleCol];

            // Split by ':'
            std::vector<std::string> sampleParts;
            {
                std::stringstream sampSS(sampleData);
                std::string part;
                while (std::getline(sampSS, part, ':')) {
                    sampleParts.push_back(part);
                }
            }
            if (gtIndex >= (int)sampleParts.size()) {
                // no GT
                continue;
            }

            std::string genotype = sampleParts[gtIndex]; // e.g. "0/1" or "2|1"
            // unify separators
            std::replace(genotype.begin(), genotype.end(), '|', '/');
            std::vector<std::string> alleleNums;
            {
                std::stringstream gtSS(genotype);
                std::string a;
                while (std::getline(gtSS, a, '/')) {
                    alleleNums.push_back(a);
                }
            }
            // We expect diploid => 2 allele calls, but handle if 1 or 2
            if (alleleNums.empty()) {
                continue;
            }
            // For each allele: 0 => REF, 1 => altAlleles[0], 2 => altAlleles[1], etc.
            for (auto &aStr : alleleNums) {
                if (aStr.empty() || aStr == ".") {
                    continue; // missing
                }
                bool numeric = true;
                for (char c : aStr) {
                    if (!isdigit(c)) {
                        numeric = false;
                        break;
                    }
                }
                if (!numeric) {
                    continue;
                }
                int aVal = std::stoi(aStr);
                if (aVal == 0) {
                    // ref => build a key with the REF as alt? 
                    // Typically we want alt freq, but let's say we do have lines for "chr:pos:ref:REF"? 
                    // Usually ancestry freq is about alt. If you want to incorporate ref, you'd store that
                    // in the freq file. We'll skip if we only track alt freqs. 
                    // For demonstration, let's skip 0 -> no alt allele => no ancestry added.
                    continue;
                }
                if (aVal > 0 && (size_t)aVal <= altAlleles.size()) {
                    // This alt
                    std::string actualAlt = altAlleles[aVal - 1];
                    // Now check each population freq. 
                    // We store freq in freqData keyed by "chr:pos:ref:alt:pop".
                    // We do not average them; we add the freq for each population? 
                    // Usually you'd pick the population with the highest freq. 
                    // Another approach is to add freq to that population's score. 
                    // Let's do "score += freq" for that population only. 
                    // That requires we look up each pop? 
                    // But we only have a single freq entry per pop in freqData. 
                    // Let's do a pass over freqData for all populations that have a key "chr:pos:ref:actualAlt:POP".
                    // Then we add that freq to sampleScores. This is a “sum all populations?” That doesn't make sense. 
                    // Usually you'd pick the single population with the highest freq. 
                    // Alternatively, you can add to them *all*, weighting by their freq. 
                    // But the code in question used "max freq" approach or "all freq"? 
                    // We'll do the "max freq" approach, consistent with #7 / #8 style.
                    
                    double bestFreq = -1.0;
                    std::string bestPop;
                    
                    // We must search each possible pop in freqData. But we store them individually as keys with pop at the end.
                    // Let's build a pattern "chrom:pos:ref:actualAlt:" and iterate over possible pops? 
                    // There's no direct iteration over freqData by prefix. So let's do a small trick:
                    
                    // We'll build the prefix:
                    std::string prefix = chrom + ":" + pos + ":" + ref + ":" + actualAlt + ":";
                    
                    // Now we can look up each population's freq by prefix + pop if we had a separate list of populations.
                    // But we do not store the list of populations here. So let's do an approach:
                    // We'll iterate over freqData and check if it starts with prefix, extracting the pop. This is not super efficient,
                    // but simple for demonstration. Or we can store a separate structure. 
                    // For a large dataset, we should store a multi-level structure or keep a set of populations. 
                    // We'll do the iteration approach for clarity.
                    
                    for (auto &kv : freqData) {
                        const std::string &fullKey = kv.first; // e.g. "chr1:12345:A:G:EUR"
                        if (fullKey.rfind(prefix, 0) == 0) {
                            // Extract the pop from the remainder
                            // prefix length is prefix.size(). The pop is what's after that
                            std::string pop = fullKey.substr(prefix.size());
                            double populationFreq = kv.second;
                            if (populationFreq > bestFreq) {
                                bestFreq = populationFreq;
                                bestPop = pop;
                            }
                        }
                    }
                    
                    if (!bestPop.empty() && bestFreq >= 0.0) {
                        // We add the bestFreq to the sample's pop
                        sampleScores[sampleNames[s]][bestPop] += bestFreq;
                    }
                } 
                // else if aVal > altAlleles.size() => skip 
            }
        }
    }

    // Done reading VCF. Now we pick the best population for each sample.
    out << "Sample\tInferred_Population\n";
    for (auto &sName : sampleNames) {
        // sampleScores[sName] => map<pop, score>
        auto it = sampleScores.find(sName);
        if (it == sampleScores.end()) {
            // no data => unknown
            out << sName << "\tUnknown\n";
            continue;
        }
        const auto &popMap = it->second;
        std::string bestPop = "Unknown";
        double bestScore = -1.0;
        for (auto &ps : popMap) {
            if (ps.second > bestScore) {
                bestScore = ps.second;
                bestPop = ps.first;
            }
        }
        out << sName << "\t" << bestPop << "\n";
    }
}


./VCFX_ancestry_inferrer/VCFX_ancestry_inferrer.h
#ifndef VCFX_ANCESTRY_INFERRER_H
#define VCFX_ANCESTRY_INFERRER_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// VCFXAncestryInferer: Header file for Ancestry Inference tool
class VCFXAncestryInferer {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Loads population allele frequencies from a file
    bool loadPopulationFrequencies(const std::string& freqFilePath);

    // Infers ancestry for each sample based on allele frequencies
    void inferAncestry(std::istream& vcfInput, std::ostream& ancestryOutput);

    // Structure to hold allele frequencies per population
    struct PopulationFrequencies {
        std::string population;
        std::unordered_map<std::string, double> variantFrequencies; // Key: "chr:pos:ref:alt", Value: frequency
    };

    std::vector<PopulationFrequencies> populations;
};

#endif // VCFX_ANCESTRY_INFERRER_H


./VCFX_annotation_extractor/VCFX_annotation_extractor.cpp
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <cstdlib>
#include <getopt.h>

// --------------------------------------------------------------
// A small struct to store command-line options
// --------------------------------------------------------------
struct AnnotationOptions {
    std::vector<std::string> annotations; // e.g. ["ANN", "Gene"]
};

// --------------------------------------------------------------
// Utility: split a string by a delimiter into a vector
// --------------------------------------------------------------
static std::vector<std::string> split(const std::string &str, char delim) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}

// --------------------------------------------------------------
// Utility: parse the INFO field into a map key->value
//   e.g. "ANN=xxx;Gene=YYY;DP=100" => {ANN:xxx, Gene:YYY, DP:100}
// --------------------------------------------------------------
static std::unordered_map<std::string, std::string> parseInfoToMap(const std::string &info) {
    std::unordered_map<std::string, std::string> infoMap;
    // split by ';'
    auto fields = split(info, ';');
    for (auto &f : fields) {
        if (f.empty()) {
            continue;
        }
        // e.g. f="ANN=..."
        auto eqPos = f.find('=');
        if (eqPos == std::string::npos) {
            // key without value? e.g. "SOMATIC"
            // You could store it as {SOMATIC: ""} if you want
            infoMap[f] = "";
            continue;
        }
        std::string key = f.substr(0, eqPos);
        std::string val = f.substr(eqPos + 1);
        infoMap[key] = val;
    }
    return infoMap;
}

// --------------------------------------------------------------
// Show usage/help
// --------------------------------------------------------------
static void printHelp() {
    std::cout 
        << "VCFX_annotation_extractor: Extract variant annotations from a VCF file.\n\n"
        << "Usage:\n"
        << "  VCFX_annotation_extractor --annotation-extract \"ANN,Gene\" < input.vcf > out.tsv\n\n"
        << "Options:\n"
        << "  -a, --annotation-extract   Comma-separated list of annotations to extract (e.g., ANN,Gene)\n"
        << "  -h, --help                 Display this help message and exit\n\n"
        << "Description:\n"
        << "  Reads a VCF from stdin and prints a tab-delimited output. For multi-ALT\n"
        << "  lines, each ALT allele is printed on its own line. If an annotation field (like\n"
        << "  'ANN=') has multiple comma-separated sub-entries, we attempt to align them with\n"
        << "  the ALT alleles in order.\n\n"
        << "Example:\n"
        << "  VCFX_annotation_extractor --annotation-extract \"ANN,Gene\" < input.vcf > out.tsv\n";
}

// --------------------------------------------------------------
// parseArguments: fill in AnnotationOptions
// --------------------------------------------------------------
static bool parseArguments(int argc, char* argv[], AnnotationOptions &opts) {
    bool showHelp = false;

    static struct option long_options[] = {
        {"annotation-extract", required_argument, 0, 'a'},
        {"help",               no_argument,       0, 'h'},
        {0,                    0,                 0,  0 }
    };

    while (true) {
        int optIdx = 0;
        int c = getopt_long(argc, argv, "a:h", long_options, &optIdx);
        if (c == -1) break;
        switch (c) {
            case 'a': {
                // comma-separated annotation names
                auto items = split(optarg, ',');
                for (auto &it : items) {
                    // trim spaces
                    while (!it.empty() && (it.front() == ' ' || it.front() == '\t')) {
                        it.erase(it.begin());
                    }
                    while (!it.empty() && (it.back() == ' ' || it.back() == '\t')) {
                        it.pop_back();
                    }
                    opts.annotations.push_back(it);
                }
            } break;
            case 'h':
            default:
                showHelp = true;
                break;
        }
    }

    if (showHelp) {
        printHelp();
        // Return false to indicate we should exit.
        return false;
    }
    // If no annotations, also show help
    if (opts.annotations.empty()) {
        printHelp();
        return false;
    }
    return true;
}

// --------------------------------------------------------------
// main extraction logic
//   1) Read VCF, pass headers unchanged
//   2) For data lines, parse ALT, parse INFO
//   3) For each ALT, line up any multi-comma annotation
//   4) Print to stdout: CHROM POS ID REF ALT [annotation1, annotation2, ...]
// --------------------------------------------------------------
static void processVCF(std::istream &in, const AnnotationOptions &opts) {
    std::string line;
    bool foundChromHeader = false;

    // Print a single header row for the TSV:
    // e.g. CHROM POS ID REF ALT ANN Gene ...
    // We'll do this only after we know how many annotation fields.
    // Or we can do it right away:
    std::cout << "CHROM\tPOS\tID\tREF\tALT";
    for (auto &annName : opts.annotations) {
        std::cout << "\t" << annName;
    }
    std::cout << "\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        // If header line
        if (line[0] == '#') {
            // If you want to preserve VCF headers, you can print them. 
            // But here we do NOT, since we produce a new TSV. 
            // If you prefer, uncomment below:
            // std::cout << line << "\n";
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
            }
            continue;
        }

        // We expect #CHROM line before data
        if (!foundChromHeader) {
            std::cerr << "Warning: Data encountered before #CHROM header: skipping\n";
            continue;
        }

        // Split the line by tabs
        // Minimal VCF => 8 fields: CHROM POS ID REF ALT QUAL FILTER INFO
        // plus optional FORMAT + samples
        auto fields = split(line, '\t');
        if (fields.size() < 8) {
            std::cerr << "Warning: Invalid VCF line (fewer than 8 fields): " << line << "\n";
            continue;
        }

        const std::string &chrom = fields[0];
        const std::string &pos   = fields[1];
        const std::string &id    = fields[2];
        const std::string &ref   = fields[3];
        const std::string &altStr= fields[4];
        // skip fields[5]=QUAL, fields[6]=FILTER
        const std::string &info  = fields[7];

        // Split ALT on comma for multiple alt alleles
        std::vector<std::string> alts = split(altStr, ',');

        // Parse INFO into a map: KEY -> value
        auto infoMap = parseInfoToMap(info);

        // For each annotation requested, retrieve the value or "NA"
        // BUT if the annotation might have multiple comma-separated sub-entries (like ANN=),
        // we handle that separately:
        // We'll do a 2-phase approach:
        //   1) gather strings for each annotation 
        //   2) if it's something like "ANN" that has multiple commas, we split them 
        //      and try to align them with ALT alleles.
        std::unordered_map<std::string, std::string> rawAnnValues;
        for (auto &annName : opts.annotations) {
            auto it = infoMap.find(annName);
            if (it == infoMap.end()) {
                rawAnnValues[annName] = "NA";
            } else {
                // entire string for that annotation
                rawAnnValues[annName] = it->second; 
            }
        }

        // Now let's produce lines, one per ALT. 
        // For each alt, we figure out the sub-annotation for any multi-value fields
        // like "ANN=val1,val2,val3" if we have 3 alts. 
        // We'll do that if rawAnnValues["ANN"] != "NA" 
        // We'll split by ',' and pick sub-annotation i for alt i 
        // If out-of-range, "NA."
        // We'll do the same for any annotation that looks comma separated. 
        // If your annotation is not meant to align with alt, you might skip this logic 
        // or check a list of known multi-value keys. 
        // We'll assume all requested keys might be multi-value.

        // For each annotation, we split on ',' 
        // We store them in a vector. Then for alt index i, we pick subAnn[i] if exists, else "NA."
        std::unordered_map<std::string, std::vector<std::string>> splittedAnnValues;
        for (auto &annName : opts.annotations) {
            std::string val = rawAnnValues[annName];
            if (val == "NA") {
                // no annotation => just store empty vector 
                splittedAnnValues[annName] = {};
                continue;
            }
            // split by ','
            auto subVals = split(val, ',');
            splittedAnnValues[annName] = subVals;
        }

        // Now produce lines
        for (size_t altIndex = 0; altIndex < alts.size(); ++altIndex) {
            // alt allele
            const std::string &thisAlt = alts[altIndex];

            // We'll prepare a line with columns: CHROM, POS, ID, REF, ALT, then each annotation
            // For each annotation, we see if splittedAnnValues[annName].size() > altIndex
            // If yes, output that sub-value, else "NA"
            std::ostringstream outLine;
            outLine << chrom << "\t" << pos << "\t" << id << "\t" 
                    << ref << "\t" << thisAlt;

            // Now each annotation
            for (auto &annName : opts.annotations) {
                const auto &subVals = splittedAnnValues[annName];
                std::string outVal = "NA";
                if (altIndex < subVals.size() && !subVals[altIndex].empty()) {
                    outVal = subVals[altIndex];
                }
                outLine << "\t" << outVal;
            }

            // Print
            std::cout << outLine.str() << "\n";
        }
    }
}

// --------------------------------------------------------------
// main()
// --------------------------------------------------------------
int main(int argc, char* argv[]) {
    AnnotationOptions opts;
    if (!parseArguments(argc, argv, opts)) {
        // parseArguments already printed help if needed
        return 0; 
    }
    // Process from stdin => produce to stdout
    processVCF(std::cin, opts);
    return 0;
}


./VCFX_annotation_extractor/VCFX_annotation_extractor.h
#ifndef VCFX_ANNOTATION_EXTRACTOR_H
#define VCFX_ANNOTATION_EXTRACTOR_H

#include <iostream>
#include <string>
#include <vector>

// VCFX_annotation_extractor: Header file for variant annotation extraction tool
class VCFXAnnotationExtractor {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Processes VCF input from the given stream
    void processVCF(std::istream& in, const std::vector<std::string>& annotations);

    // Parses annotations from the INFO field
    std::vector<std::string> parseINFO(const std::string& info);

    // Extracts specified annotations
    std::vector<std::string> extractAnnotations(const std::vector<std::string>& info_fields, const std::vector<std::string>& annotations);
};

#endif // VCFX_ANNOTATION_EXTRACTOR_H

./VCFX_compressor/VCFX_compressor.cpp
#include <iostream>
#include <string>
#include <zlib.h>
#include <cstring>
#include <getopt.h>

// ---------------------------------------------------------------------------
// Show help
// ---------------------------------------------------------------------------
static void printHelp() {
    std::cout 
        << "VCFX_compressor\n"
        << "Usage: VCFX_compressor [OPTIONS]\n\n"
        << "Options:\n"
        << "  --compress, -c         Compress the input VCF file (to stdout).\n"
        << "  --decompress, -d       Decompress the input VCF file (from stdin).\n"
        << "  --help, -h             Display this help message and exit.\n\n"
        << "Description:\n"
        << "  Compresses or decompresses data using zlib's raw DEFLATE (similar to gzip).\n"
        << "  Note that for .vcf.gz indexing via tabix, one typically needs BGZF blocks,\n"
        << "  which is not implemented here.\n\n"
        << "Examples:\n"
        << "  Compress:\n"
        << "    ./VCFX_compressor --compress < input.vcf > output.vcf.gz\n\n"
        << "  Decompress:\n"
        << "    ./VCFX_compressor --decompress < input.vcf.gz > output.vcf\n";
}

// ---------------------------------------------------------------------------
// compressDecompressVCF
//   compress = true  => read from 'in', produce gzip to 'out'
//   compress = false => read gzip from 'in', produce plain text to 'out'
// ---------------------------------------------------------------------------
static bool compressDecompressVCF(std::istream& in, std::ostream& out, bool compress) {
    constexpr int CHUNK = 16384;
    char inBuffer[CHUNK];
    char outBuffer[CHUNK];

    if (compress) {
        // Initialize for deflate
        z_stream strm;
        std::memset(&strm, 0, sizeof(strm));
        if (deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
            std::cerr << "Error: deflateInit2 failed.\n";
            return false;
        }

        // read -> deflate -> out
        int flush = Z_NO_FLUSH;
        do {
            in.read(inBuffer, CHUNK);
            std::streamsize bytesRead = in.gcount();
            flush = in.eof() ? Z_FINISH : Z_NO_FLUSH;

            strm.avail_in = static_cast<uInt>(bytesRead);
            strm.next_in  = reinterpret_cast<Bytef*>(inBuffer);

            // compress until input is used up
            do {
                strm.avail_out = CHUNK;
                strm.next_out  = reinterpret_cast<Bytef*>(outBuffer);

                int ret = deflate(&strm, flush);
                if (ret == Z_STREAM_ERROR) {
                    std::cerr << "Error: deflate failed.\n";
                    deflateEnd(&strm);
                    return false;
                }
                // # of bytes written
                size_t have = CHUNK - strm.avail_out;
                if (have > 0) {
                    out.write(outBuffer, have);
                    if (!out.good()) {
                        std::cerr << "Error: write to output stream failed.\n";
                        deflateEnd(&strm);
                        return false;
                    }
                }
            } while (strm.avail_out == 0);

        } while (flush != Z_FINISH);

        deflateEnd(&strm);
        return true;

    } else {
        // Initialize for inflate
        z_stream strm;
        std::memset(&strm, 0, sizeof(strm));
        // 15+32 to allow auto-detect of gzip/zlib
        if (inflateInit2(&strm, 15 + 32) != Z_OK) {
            std::cerr << "Error: inflateInit2 failed.\n";
            return false;
        }

        int ret = Z_OK;
        do {
            in.read(inBuffer, CHUNK);
            std::streamsize bytesRead = in.gcount();
            if (bytesRead == 0 && !in.eof()) {
                // Possibly an I/O error
                if (in.bad()) {
                    std::cerr << "Error: reading input stream.\n";
                    inflateEnd(&strm);
                    return false;
                }
            }
            strm.avail_in = static_cast<uInt>(bytesRead);
            strm.next_in  = reinterpret_cast<Bytef*>(inBuffer);

            // decompress until output buffer not needed
            do {
                strm.avail_out = CHUNK;
                strm.next_out  = reinterpret_cast<Bytef*>(outBuffer);

                ret = inflate(&strm, Z_NO_FLUSH);
                if (ret == Z_STREAM_ERROR || ret == Z_NEED_DICT ||
                    ret == Z_DATA_ERROR || ret == Z_MEM_ERROR)
                {
                    std::cerr << "Error: inflate failed with code " << ret << "\n";
                    inflateEnd(&strm);
                    return false;
                }
                size_t have = CHUNK - strm.avail_out;
                if (have > 0) {
                    out.write(outBuffer, have);
                    if (!out.good()) {
                        std::cerr << "Error: write to output stream failed.\n";
                        inflateEnd(&strm);
                        return false;
                    }
                }
            } while (strm.avail_out == 0);

        } while (ret != Z_STREAM_END && !in.eof());

        // Clean up
        inflateEnd(&strm);
        // If we didn't get Z_STREAM_END, the compressed data might be incomplete
        if (ret != Z_STREAM_END) {
            // Possibly truncated
            std::cerr << "Warning: stream ended prematurely or was truncated.\n";
        }
        return true;
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    bool compress = false;
    bool decompress = false;

    static struct option long_options[] = {
        {"compress",   no_argument, 0, 'c'},
        {"decompress", no_argument, 0, 'd'},
        {"help",       no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    while (true) {
        int optIndex = 0;
        int c = getopt_long(argc, argv, "cdh", long_options, &optIndex);
        if (c == -1) break;

        switch (c) {
            case 'c': compress = true;       break;
            case 'd': decompress = true;     break;
            case 'h': printHelp(); return 0;
            default:  printHelp(); return 1;
        }
    }

    if ((compress && decompress) || (!compress && !decompress)) {
        std::cerr << "Error: must specify exactly one of --compress or --decompress.\n";
        return 1;
    }

    if (!compressDecompressVCF(std::cin, std::cout, compress)) {
        std::cerr << "Error: Compression/Decompression failed.\n";
        return 1;
    }
    return 0;
}


./VCFX_compressor/VCFX_compressor.h
#ifndef VCFX_COMPRESSOR_H
#define VCFX_COMPRESSOR_H

#include <iostream>
#include <string>

// Function to display help message
void printHelp();

// Function to perform compression or decompression
bool compressDecompressVCF(std::istream& in, std::ostream& out, bool compress);

#endif // VCFX_COMPRESSOR_H


./VCFX_concordance_checker/VCFX_concordance_checker.cpp
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <cstdlib>
#include <getopt.h>

// ---------------------------------------------------------
// Data structure for command-line arguments
// ---------------------------------------------------------
struct ConcordanceArguments {
    std::string sample1;
    std::string sample2;
};

// ---------------------------------------------------------
// Print help
// ---------------------------------------------------------
static void printHelp() {
    std::cout 
        << "VCFX_concordance_checker\n"
        << "Usage: VCFX_concordance_checker [OPTIONS] < input.vcf > concordance_report.tsv\n\n"
        << "Options:\n"
        << "  --samples, -s \"Sample1 Sample2\"  Specify exactly two sample names to compare.\n"
        << "  --help, -h                        Display this help message and exit.\n\n"
        << "Description:\n"
        << "  Compares genotypes between two specified samples in a VCF file, including multi-allelic\n"
        << "  variants, and outputs per-variant concordance (Concordant or Discordant).\n\n"
        << "Example:\n"
        << "  ./VCFX_concordance_checker --samples \"SampleA SampleB\" < input.vcf > concordance_report.tsv\n";
}

// ---------------------------------------------------------
// Utility: split a string by delimiter
// ---------------------------------------------------------
static std::vector<std::string> splitString(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string temp;
    while (std::getline(ss, temp, delimiter)) {
        tokens.push_back(temp);
    }
    return tokens;
}

// ---------------------------------------------------------
// parseArguments: fill ConcordanceArguments
// ---------------------------------------------------------
static bool parseArguments(int argc, char* argv[], ConcordanceArguments& args) {
    bool showHelp = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--samples" || arg == "-s") && i + 1 < argc) {
            std::string samples_str = argv[++i];
            auto samples = splitString(samples_str, ' ');
            if (samples.size() != 2) {
                std::cerr << "Error: Please specify exactly two sample names (e.g., -s \"Sample1 Sample2\").\n";
                return false;
            }
            args.sample1 = samples[0];
            args.sample2 = samples[1];
        }
        else if (arg == "--help" || arg == "-h") {
            showHelp = true;
        }
        else {
            std::cerr << "Warning: Unrecognized argument '" << arg << "'.\n";
        }
    }

    if (showHelp) {
        printHelp();
        // Return false so main() can exit cleanly
        return false;
    }
    if (args.sample1.empty() || args.sample2.empty()) {
        std::cerr << "Error: Two sample names must be specified using --samples or -s.\n";
        std::cerr << "Use --help for usage information.\n";
        return false;
    }

    return true;
}

// ---------------------------------------------------------
// Convert a numeric allele index to a textual tag:
//   0 => "0"
//   1 => "1" (meaning altAlleles[0])
//   2 => "2" (meaning altAlleles[1]), etc.
// We'll store them as strings "0","1","2" for comparison
// ---------------------------------------------------------
static std::string normalizeDiploidGenotype(const std::string& genotypeField,
                                            const std::vector<std::string> &altAlleles)
{
    // genotypeField might be "0/1", "2|3", ".", ...
    // Step 1: replace '|' with '/'
    std::string gt = genotypeField;
    for (char &c : gt) {
        if (c == '|') c = '/';
    }
    // split by '/'
    auto alleles = splitString(gt, '/');
    if (alleles.size() < 2) {
        // not diploid or missing
        return "";
    }
    // parse each allele
    std::vector<int> numericAll;
    numericAll.reserve(2);
    for (auto &a : alleles) {
        if (a == "." || a.empty()) {
            // missing
            return "";
        }
        // must parse as an integer
        bool numeric = std::all_of(a.begin(), a.end(), ::isdigit);
        if (!numeric) {
            return "";
        }
        int val = std::stoi(a);
        // If val > altAlleles.size(), we can't interpret => skip
        // val==0 => reference
        // val==1 => altAlleles[0]
        // val==2 => altAlleles[1], etc.
        if (val < 0) {
            return "";
        }
        if (val > 0) {
            // if val=2 => we need altAlleles.size() >=2
            if ((size_t)val > altAlleles.size()) {
                // out of range
                return "";
            }
        }
        numericAll.push_back(val);
    }
    // We produce a string "x/y" sorted or not? Typically for genotype comparison, 
    // "1/0" is same as "0/1", but let's keep the original order. It's also common 
    // to store them in numeric order. We'll store them in sorted numeric order 
    // so that "0/1" == "1/0". Then we join with "/". 
    std::sort(numericAll.begin(), numericAll.end());
    std::ostringstream oss;
    oss << numericAll[0] << "/" << numericAll[1];
    return oss.str();
}

// ---------------------------------------------------------
// Calculate concordance
//   1) parse the #CHROM line, find sample1_index + sample2_index
//   2) for each variant, parse ALT, parse sample1 + sample2 genotype => normalized
//   3) if both genotypes are non-empty => compare => print row
//   4) track summary stats
// ---------------------------------------------------------
static bool calculateConcordance(std::istream &in, std::ostream &out, const ConcordanceArguments &args) {
    std::string line;
    bool foundChromHeader = false;
    std::vector<std::string> headerFields;
    int sample1_index = -1;
    int sample2_index = -1;

    // We'll track summary stats
    int totalVariants = 0;
    int concordant = 0;
    int discordant = 0;

    // Print a TSV header:
    // CHROM POS ID REF ALT(S) Sample1_GT Sample2_GT Concordance
    out << "CHROM\tPOS\tID\tREF\tALT\t" 
        << args.sample1 << "_GT\t" 
        << args.sample2 << "_GT\tConcordance\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            // If it's the #CHROM line, parse sample columns
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                headerFields = splitString(line, '\t');
                // find sample1, sample2
                std::unordered_map<std::string, int> sampleMap;
                for (size_t i = 9; i < headerFields.size(); ++i) {
                    sampleMap[headerFields[i]] = (int)i;
                }
                // check existence
                auto it1 = sampleMap.find(args.sample1);
                if (it1 == sampleMap.end()) {
                    std::cerr << "Error: Sample '" << args.sample1 << "' not found in VCF header.\n";
                    return false;
                }
                auto it2 = sampleMap.find(args.sample2);
                if (it2 == sampleMap.end()) {
                    std::cerr << "Error: Sample '" << args.sample2 << "' not found in VCF header.\n";
                    return false;
                }
                sample1_index = it1->second;
                sample2_index = it2->second;
                foundChromHeader = true;
            }
            continue;
        }

        if (!foundChromHeader) {
            std::cerr << "Error: VCF data encountered before #CHROM header.\n";
            return false;
        }

        // parse data line
        auto fields = splitString(line, '\t');
        if (fields.size() < 8) {
            // invalid
            continue;
        }
        // check we have sample columns for sample1, sample2
        if ((size_t)sample1_index >= fields.size() || 
            (size_t)sample2_index >= fields.size()) 
        {
            // line doesn't have enough columns
            continue;
        }

        const std::string &chrom = fields[0];
        const std::string &pos   = fields[1];
        const std::string &id    = fields[2];
        const std::string &ref   = fields[3];
        const std::string &alt   = fields[4];
        // sample columns:
        const std::string &sample1_field = fields[sample1_index];
        const std::string &sample2_field = fields[sample2_index];

        // Split the alt string (but for genotype normalization we only need the count)
        auto altAlleles = splitString(alt, ',');

        // Convert sample1 genotype to a normalized diploid string
        std::string s1_gt = normalizeDiploidGenotype(sample1_field, altAlleles);
        std::string s2_gt = normalizeDiploidGenotype(sample2_field, altAlleles);

        // If either is missing or un-parseable, skip
        if (s1_gt.empty() || s2_gt.empty()) {
            continue;
        }

        totalVariants++;
        bool same = (s1_gt == s2_gt);
        if (same) {
            concordant++;
        } else {
            discordant++;
        }
        std::string cc = same ? "Concordant" : "Discordant";

        // Print the row
        out << chrom << "\t" 
            << pos   << "\t" 
            << id    << "\t" 
            << ref   << "\t" 
            << alt   << "\t"
            << s1_gt << "\t"
            << s2_gt << "\t"
            << cc    << "\n";
    }

    // Print summary to stderr so we don't break the TSV
    std::cerr << "Total Variants Compared: " << totalVariants << "\n"
              << "Concordant Genotypes: " << concordant << "\n"
              << "Discordant Genotypes: " << discordant << "\n";

    return true;
}

// ---------------------------------------------------------
// main
// ---------------------------------------------------------
int main(int argc, char* argv[]) {
    ConcordanceArguments args;
    if (!parseArguments(argc, argv, args)) {
        // parseArguments prints error/help if needed
        return 1;
    }

    bool ok = calculateConcordance(std::cin, std::cout, args);
    return (ok ? 0 : 1);
}


./VCFX_concordance_checker/VCFX_concordance_checker.h
#ifndef VCFX_CONCORDANCE_CHECKER_H
#define VCFX_CONCORDANCE_CHECKER_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Structure to hold command-line arguments
struct ConcordanceArguments {
    std::string sample1;
    std::string sample2;
};

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], ConcordanceArguments& args);

// Function to split a string by a delimiter
std::vector<std::string> splitString(const std::string& str, char delimiter);

// Function to process VCF and calculate concordance
bool calculateConcordance(std::istream& in, std::ostream& out, const ConcordanceArguments& args);

// Function to extract genotype from genotype string
std::string extractGenotype(const std::string& genotype_str);

#endif // VCFX_CONCORDANCE_CHECKER_H


./VCFX_cross_sample_concordance/VCFX_cross_sample_concordance.cpp
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <getopt.h>

// --------------------------------------------------------------------------
// A helper function to print usage/help
// --------------------------------------------------------------------------
static void displayHelp() {
    std::cout
        << "VCFX_cross_sample_concordance: Check variant concordance across multiple samples.\n\n"
        << "Usage:\n"
        << "  VCFX_cross_sample_concordance [options] < input.vcf > concordance_results.txt\n\n"
        << "Options:\n"
        << "  -h, --help    Display this help message and exit\n\n"
        << "Description:\n"
        << "  Reads a multi-sample VCF from stdin, normalizes each sample's genotype\n"
        << "  (including multi-allelic variants), and determines if all samples that\n"
        << "  have a parseable genotype are in complete agreement. Prints one row per\n"
        << "  variant to stdout and a final summary to stderr.\n\n"
        << "Example:\n"
        << "  VCFX_cross_sample_concordance < input.vcf > results.tsv\n\n";
}

// --------------------------------------------------------------------------
// Utility: split a string by a delimiter
// --------------------------------------------------------------------------
static std::vector<std::string> split(const std::string &str, char delim) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}

// --------------------------------------------------------------------------
// Normalize a genotype string (e.g. "2|1", "0/1") to a canonical form
//   - Replace '|' with '/'
//   - Split by '/'
//   - Check for missing or invalid (out-of-range) => return empty
//   - Convert to int, store in vector
//   - Sort so "2/1" => "1/2"
//   - Return e.g. "1/2"
// --------------------------------------------------------------------------
static std::string normalizeGenotype(const std::string &sampleField,
                                     const std::vector<std::string> &altAlleles)
{
    // sampleField might be something like "0/1:.." or just "0/1"
    // We only care about the genotype portion (first colon-delimited field).
    auto parts = split(sampleField, ':');
    if (parts.empty()) {
        return "";
    }
    std::string gt = parts[0]; // e.g. "0/1", "2|3", etc.

    // unify separators
    for (char &c : gt) {
        if (c == '|') c = '/';
    }
    // split by '/'
    auto alleleStrs = split(gt, '/');
    if (alleleStrs.size() < 2) {
        // might be haploid or invalid
        return "";
    }

    std::vector<int> alleleInts;
    alleleInts.reserve(alleleStrs.size());
    for (auto &a : alleleStrs) {
        if (a == "." || a.empty()) {
            return ""; // missing
        }
        // parse
        for (char c : a) {
            if (!std::isdigit(c)) {
                return "";
            }
        }
        int val = std::stoi(a);
        // If val == 0 => REF
        // If val == 1 => altAlleles[0]
        // If val == 2 => altAlleles[1], etc.
        if (val < 0) {
            return "";
        }
        if (val > 0 && (size_t)val > altAlleles.size()) {
            // out of range for the alt array
            return "";
        }
        alleleInts.push_back(val);
    }

    // For genotype comparison, sort them so "2/1" => "1/2"
    std::sort(alleleInts.begin(), alleleInts.end());
    // Build a normalized string
    std::stringstream out;
    for (size_t i = 0; i < alleleInts.size(); ++i) {
        if (i > 0) out << "/";
        out << alleleInts[i];
    }
    return out.str();
}

// --------------------------------------------------------------------------
// The main function to process the VCF and produce cross-sample concordance
// --------------------------------------------------------------------------
static void calculateConcordance(std::istream &in, std::ostream &out) {
    std::string line;
    std::vector<std::string> sampleNames;
    bool gotChromHeader = false;

    // We'll print a header for our TSV
    // CHROM POS ID REF ALT Num_Samples Unique_Normalized_Genotypes Concordance_Status
    out << "CHROM\tPOS\tID\tREF\tALT\tNum_Samples\tUnique_Normalized_Genotypes\tConcordance_Status\n";

    // We track summary stats
    size_t totalVariants = 0;
    size_t allConcordantCount = 0;
    size_t discordantCount = 0;
    size_t skippedBecauseNoGenotypes = 0; // if all samples are missing/un-parseable

    // 1) Read until we find the #CHROM line
    //    then parse sample names from columns 9+
    while (!gotChromHeader && std::getline(in, line)) {
        if (line.rfind("#CHROM", 0) == 0) {
            // parse columns
            auto headers = split(line, '\t');
            // from index=9 onwards => sample names
            for (size_t i = 9; i < headers.size(); ++i) {
                sampleNames.push_back(headers[i]);
            }
            gotChromHeader = true;
            break; // proceed to reading data lines
        }
    }
    if (!gotChromHeader) {
        std::cerr << "Error: VCF header with #CHROM not found.\n";
        return;
    }

    // 2) Now read the data lines
    while (std::getline(in, line)) {
        // skip header lines or empties
        if (line.empty() || line[0] == '#') {
            continue;
        }
        auto fields = split(line, '\t');
        // minimal VCF => 8 columns + samples => need at least 10 if there's 1 sample
        if (fields.size() < 8) {
            // invalid line
            continue;
        }

        // parse standard columns
        const std::string &chrom = fields[0];
        const std::string &pos   = fields[1];
        const std::string &id    = fields[2];
        const std::string &ref   = fields[3];
        const std::string &alt   = fields[4];
        // we can skip qual/filter/info/format => we only need sample columns
        // sample columns start at index=9
        if (fields.size() < (9 + sampleNames.size())) {
            // line doesn't have enough columns for all samples
            // skip
            continue;
        }

        // parse alt for multi-allelic
        auto altAlleles = split(alt, ',');

        // For each sample, normalize genotype
        std::vector<std::string> normalizedGts;
        normalizedGts.reserve(sampleNames.size());
        size_t sampleCountHere = 0;

        for (size_t s = 0; s < sampleNames.size(); ++s) {
            auto &sampleField = fields[9 + s];
            std::string norm = normalizeGenotype(sampleField, altAlleles);
            if (!norm.empty()) {
                // we have a parseable genotype
                normalizedGts.push_back(norm);
                sampleCountHere++;
            }
        }

        // If no samples had a valid genotype, skip counting this variant
        if (sampleCountHere == 0) {
            skippedBecauseNoGenotypes++;
            continue;
        }

        // We have at least 1 genotype
        totalVariants++;

        // gather unique genotype strings
        std::unordered_map<std::string, int> freq;
        for (auto &g : normalizedGts) {
            freq[g]++;
        }
        size_t uniqueGtCount = freq.size();
        bool allConcordant = (uniqueGtCount == 1);

        if (allConcordant) {
            allConcordantCount++;
        } else {
            discordantCount++;
        }

        std::string status = allConcordant ? "Concordant" : "Discordant";
        // Print row
        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt
            << "\t" << sampleCountHere 
            << "\t" << uniqueGtCount
            << "\t" << status
            << "\n";
    }

    // Summaries to stderr so as not to pollute the TSV
    std::cerr << "Total Variants with >=1 parseable genotype: " << totalVariants << "\n";
    std::cerr << "   Concordant (all same genotype): " << allConcordantCount << "\n";
    std::cerr << "   Discordant (>=2 distinct genotypes): " << discordantCount << "\n";
    std::cerr << "Variants with no parseable genotypes (skipped): " << skippedBecauseNoGenotypes << "\n";
}

// --------------------------------------------------------------------------
// Command-line parsing + main
// --------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    bool showHelp = false;

    static struct option longOpts[] = {
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };

    while (true) {
        int optIndex=0;
        int c = getopt_long(argc, argv, "h", longOpts, &optIndex);
        if (c == -1) break;
        switch (c) {
            case 'h': showHelp = true; break;
            default:  showHelp = true; break;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    calculateConcordance(std::cin, std::cout);
    return 0;
}


./VCFX_cross_sample_concordance/VCFX_cross_sample_concordance.h
#ifndef VCFX_CROSS_SAMPLE_CONCORDANCE_H
#define VCFX_CROSS_SAMPLE_CONCORDANCE_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// VCFXCrossSampleConcordance: Header file for Cross-Sample Variant Concordance Tool
class VCFXCrossSampleConcordance {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Processes the VCF input and calculates concordance
    void calculateConcordance(std::istream& in, std::ostream& out);

    // Structure to hold variant information
    struct Variant {
        std::string chrom;
        std::string pos;
        std::string ref;
        std::string alt;
        std::vector<std::string> genotypes;
    };

    // Stores all variants
    std::vector<Variant> variants;
};

#endif // VCFX_CROSS_SAMPLE_CONCORDANCE_H

./VCFX_custom_annotator/VCFX_custom_annotator.cpp
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <getopt.h>

// ---------------------------------------------------------------------------
// Class: VCFXCustomAnnotator
// ---------------------------------------------------------------------------
class VCFXCustomAnnotator {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();
    bool loadAnnotations(const std::string& annotationFilePath,
                         std::unordered_map<std::string, std::string>& annotations);
    void addAnnotations(std::istream& in, std::ostream& out,
                        const std::unordered_map<std::string, std::string>& annotations);

    // Generate a unique key for a variant
    // Format: CHROM:POS:REF:ALT
    static std::string generateVariantKey(const std::string& chrom,
                                          const std::string& pos,
                                          const std::string& ref,
                                          const std::string& alt);
};

// ---------------------------------------------------------------------------
// Utility: split a string by a delimiter
// ---------------------------------------------------------------------------
static std::vector<std::string> split(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string tmp;
    while (std::getline(ss, tmp, delimiter)) {
        tokens.push_back(tmp);
    }
    return tokens;
}

// ---------------------------------------------------------------------------
// displayHelp
// ---------------------------------------------------------------------------
void VCFXCustomAnnotator::displayHelp() {
    std::cout
        << "VCFX_custom_annotator: Add custom annotations to the INFO field in a VCF file.\n\n"
        << "Usage:\n"
        << "  VCFX_custom_annotator --add-annotation <annotations.txt> [options]\n\n"
        << "Options:\n"
        << "  -h, --help                  Display this help message and exit\n"
        << "  -a, --add-annotation <file> Specify the annotation file\n\n"
        << "Description:\n"
        << "  Reads an annotation file with lines:\n"
        << "    CHROM  POS  REF  ALT  annotation...\n"
        << "  Then for each VCF variant, if it matches CHROM:POS:REF:ALT, inserts\n"
        << "  'CustomAnnotation=...' into the INFO field.\n"
        << "  Multi-allelic ALT fields are split on commas; we attempt to annotate\n"
        << "  each ALT separately. If no annotation is found for a given ALT, 'NA'\n"
        << "  is used for that allele's slot.\n\n"
        << "Example:\n"
        << "  VCFX_custom_annotator --add-annotation annotations.txt < input.vcf > annotated.vcf\n";
}

// ---------------------------------------------------------------------------
// generateVariantKey
// ---------------------------------------------------------------------------
std::string VCFXCustomAnnotator::generateVariantKey(const std::string& chrom,
                                                    const std::string& pos,
                                                    const std::string& ref,
                                                    const std::string& alt) {
    return chrom + ":" + pos + ":" + ref + ":" + alt;
}

// ---------------------------------------------------------------------------
// loadAnnotations
//   Reads a file with lines: CHROM POS REF ALT annotation...
//   Stores them in a map: "chrom:pos:ref:alt" -> annotation
// ---------------------------------------------------------------------------
bool VCFXCustomAnnotator::loadAnnotations(
    const std::string& annotationFilePath,
    std::unordered_map<std::string, std::string>& annotations)
{
    std::ifstream infile(annotationFilePath);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open annotation file " << annotationFilePath << "\n";
        return false;
    }

    std::string line;
    size_t line_num = 0;
    while (std::getline(infile, line)) {
        line_num++;
        if (line.empty() || line[0] == '#') {
            // skip comments or empty lines
            continue;
        }
        std::stringstream ss(line);
        std::string chrom, pos, ref, alt;
        if (!(ss >> chrom >> pos >> ref >> alt)) {
            std::cerr << "Warning: Skipping invalid annotation line " << line_num
                      << ": " << line << "\n";
            continue;
        }
        // The rest of the line (after the first four fields) is the annotation text
        std::string annotation;
        if (!std::getline(ss, annotation)) {
            // It's okay if there's no extra text, but we treat it as empty
            annotation = "";
        }
        // Trim leading whitespace from annotation
        if (!annotation.empty()) {
            size_t startPos = annotation.find_first_not_of(" \t");
            if (startPos == std::string::npos) {
                // it's all whitespace
                annotation.clear();
            } else {
                annotation.erase(0, startPos);
            }
        }
        // Build key
        std::string key = generateVariantKey(chrom, pos, ref, alt);
        annotations[key] = annotation;
    }
    return true;
}

// ---------------------------------------------------------------------------
// addAnnotations
//   Reads the VCF from 'in', writes to 'out' a new header line
//   and appends "CustomAnnotation=..." to each variant's INFO if found.
//   For multi-allelic ALT, we do one lookup per allele, merging them into
//   a single string "val1,val2" if there are multiple alt alleles.
//   If no annotation is found for an alt, we use "NA" for that slot.
// ---------------------------------------------------------------------------
void VCFXCustomAnnotator::addAnnotations(std::istream& in,
                                         std::ostream& out,
                                         const std::unordered_map<std::string, std::string>& annotations)
{
    bool infoHeaderInserted = false;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        // If it's a header line
        if (line[0] == '#') {
            // If it's #CHROM and we haven't inserted the new INFO header yet, do so
            if (!infoHeaderInserted && line.rfind("#CHROM", 0) == 0) {
                out << "##INFO=<ID=CustomAnnotation,Number=.,Type=String,Description=\"Custom annotations added by VCFX_custom_annotator (multi-allelic)\">\n";
                infoHeaderInserted = true;
            }
            out << line << "\n";
            continue;
        }

        // parse the 8 standard fields + possibly more
        // e.g. CHROM POS ID REF ALT QUAL FILTER INFO ...
        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        // read the remainder of the line as one string (could include FORMAT, sample columns, etc.)
        std::string rest;
        if (std::getline(ss, rest)) {
            // ' rest ' includes the leading space if any
            // We can just keep it as is
        }

        // If ALT is multi-allelic, e.g. "A,C,G", we do a separate lookup for each allele
        auto altAlleles = split(alt, ',');

        // For each alt allele, build the key, see if it's in the annotations
        // If not found, store "NA"
        std::vector<std::string> annVals;
        annVals.reserve(altAlleles.size());
        bool anyFound = false;
        for (auto &a : altAlleles) {
            std::string key = generateVariantKey(chrom, pos, ref, a);
            auto it = annotations.find(key);
            if (it == annotations.end()) {
                annVals.push_back("NA");
            } else {
                annVals.push_back(it->second.empty() ? "NA" : it->second);
                anyFound = true;
            }
        }

        // Build a single annotation field if needed:
        // e.g. "val1,val2" or "NA,val2" etc.
        // If we want to show the user which alt is which, we do them in order.
        std::string finalAnn;
        {
            // join annVals with commas
            std::ostringstream oss;
            for (size_t i = 0; i < annVals.size(); ++i) {
                if (i > 0) oss << ",";
                oss << annVals[i];
            }
            finalAnn = oss.str();
        }

        // Insert into INFO if we want it even if all are "NA"
        // Some users might prefer skipping if all alt are "NA".
        // We'll keep it consistent and always add it, to see that no annotation was found.
        // If the original info is ".", we replace it; else append ";"
        if (info == ".") {
            info = "CustomAnnotation=" + finalAnn;
        } else {
            info += ";CustomAnnotation=" + finalAnn;
        }

        // Reconstruct the VCF line
        // We print CHROM POS ID REF ALT QUAL FILTER INFO then the rest
        // e.g. leftover might be " FORMAT SAMPLE1 SAMPLE2..."
        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt
            << "\t" << qual << "\t" << filter << "\t" << info;

        // If there's anything left in 'rest', print it
        if (!rest.empty()) {
            out << rest; // includes leading space
        }
        out << "\n";
    }
}

// ---------------------------------------------------------------------------
// run() - parse arguments, load annotation map, apply to VCF
// ---------------------------------------------------------------------------
int VCFXCustomAnnotator::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string annotationFilePath;

    static struct option long_options[] = {
        {"help",           no_argument,       0, 'h'},
        {"add-annotation", required_argument, 0, 'a'},
        {0,0,0,0}
    };

    while ((opt = getopt_long(argc, argv, "ha:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'a':
                annotationFilePath = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || annotationFilePath.empty()) {
        displayHelp();
        return (showHelp ? 0 : 1);
    }

    // Load annotation map
    std::unordered_map<std::string, std::string> annotations;
    if (!loadAnnotations(annotationFilePath, annotations)) {
        std::cerr << "Error: Failed to load annotations from " << annotationFilePath << "\n";
        return 1;
    }

    // Annotate from stdin to stdout
    addAnnotations(std::cin, std::cout, annotations);
    return 0;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    VCFXCustomAnnotator annotator;
    return annotator.run(argc, argv);
}


./VCFX_custom_annotator/VCFX_custom_annotator.h
#ifndef VCFX_CUSTOM_ANNOTATOR_H
#define VCFX_CUSTOM_ANNOTATOR_H

#include <iostream>
#include <string>
#include <unordered_map>

// VCFXCustomAnnotator: Header file for Custom Annotation Addition Tool
class VCFXCustomAnnotator {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Loads annotations from a file into a map
    bool loadAnnotations(const std::string& annotationFilePath, std::unordered_map<std::string, std::string>& annotations);

    // Adds annotations to the VCF input
    void addAnnotations(std::istream& in, std::ostream& out, const std::unordered_map<std::string, std::string>& annotations);

    // Generates a unique key for a variant based on chromosome, position, ref, and alt
    std::string generateVariantKey(const std::string& chrom, const std::string& pos, const std::string& ref, const std::string& alt);
};

#endif // VCFX_CUSTOM_ANNOTATOR_H


./VCFX_diff_tool/VCFX_diff_tool.cpp
#include <iostream>
#include <string>
#include <fstream>
#include <unordered_set>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <getopt.h>

// ----------------------------------------------------------------------
// A helper function to split a string by a delimiter
// ----------------------------------------------------------------------
static std::vector<std::string> split(const std::string &str, char delim) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string tmp;
    while (std::getline(ss, tmp, delim)) {
        tokens.push_back(tmp);
    }
    return tokens;
}

// ----------------------------------------------------------------------
// Class: VCFXDiffTool
// ----------------------------------------------------------------------
class VCFXDiffTool {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();

    bool loadVariants(const std::string& filePath,
                      std::unordered_set<std::string>& variants);

    // We unify multi-allelic lines by splitting ALT on commas,
    // sorting them, and rejoining them. 
    // The final key is "chrom:pos:ref:sortedAltString"
    std::string generateVariantKey(const std::string& chrom,
                                   const std::string& pos,
                                   const std::string& ref,
                                   const std::string& altField);
};

// ----------------------------------------------------------------------
// displayHelp
// ----------------------------------------------------------------------
void VCFXDiffTool::displayHelp() {
    std::cout << "VCFX_diff_tool: Compare two VCF files and identify differences.\n\n"
              << "Usage:\n"
              << "  VCFX_diff_tool --file1 <file1.vcf> --file2 <file2.vcf>\n\n"
              << "Options:\n"
              << "  -h, --help                Display this help message and exit\n"
              << "  -a, --file1 <file1.vcf>   Specify the first VCF file\n"
              << "  -b, --file2 <file2.vcf>   Specify the second VCF file\n\n"
              << "Example:\n"
              << "  VCFX_diff_tool --file1 file1.vcf --file2 file2.vcf\n";
}

// ----------------------------------------------------------------------
// generateVariantKey
//   - altField might be e.g. "G" or "G,T" or "T,GA"
//   - We split on commas, sort them, rejoin => e.g. "G,T"
//   - Then produce "chrom:pos:ref:thatSortedAlt"
// ----------------------------------------------------------------------
std::string VCFXDiffTool::generateVariantKey(const std::string& chrom,
                                             const std::string& pos,
                                             const std::string& ref,
                                             const std::string& altField)
{
    auto alts = split(altField, ',');
    // sort them lexicographically
    std::sort(alts.begin(), alts.end());
    // rejoin
    std::ostringstream altSS;
    for (size_t i = 0; i < alts.size(); ++i) {
        if (i > 0) altSS << ",";
        altSS << alts[i];
    }
    std::string sortedAlt = altSS.str();

    return chrom + ":" + pos + ":" + ref + ":" + sortedAlt;
}

// ----------------------------------------------------------------------
// loadVariants
//   - Read each line, skip headers (#...)
//   - parse CHROM, POS, ID, REF, ALT
//   - generate key => add to set
// ----------------------------------------------------------------------
bool VCFXDiffTool::loadVariants(const std::string& filePath,
                                std::unordered_set<std::string>& variants)
{
    std::ifstream infile(filePath);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open file " << filePath << "\n";
        return false;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') {
            // skip headers and empty lines
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt;

        // minimal parse
        if (!(ss >> chrom >> pos >> id >> ref >> alt)) {
            std::cerr << "Warning: Skipping invalid VCF line:\n" << line << "\n";
            continue;
        }

        std::string key = generateVariantKey(chrom, pos, ref, alt);
        variants.insert(key);
    }
    return true;
}

// ----------------------------------------------------------------------
// run()
//   - parse args
//   - load sets from file1, file2
//   - compare => print differences
// ----------------------------------------------------------------------
int VCFXDiffTool::run(int argc, char* argv[]) {
    int opt;
    bool showHelp = false;
    std::string file1Path;
    std::string file2Path;

    static struct option long_options[] = {
        {"help",  no_argument,       0, 'h'},
        {"file1", required_argument, 0, 'a'},
        {"file2", required_argument, 0, 'b'},
        {0,0,0,0}
    };

    while ((opt = getopt_long(argc, argv, "ha:b:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'a':
                file1Path = optarg;
                break;
            case 'b':
                file2Path = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || file1Path.empty() || file2Path.empty()) {
        displayHelp();
        return (showHelp ? 0 : 1);
    }

    // Load variants
    std::unordered_set<std::string> variantsFile1;
    std::unordered_set<std::string> variantsFile2;

    if (!loadVariants(file1Path, variantsFile1)) {
        std::cerr << "Error: Failed to load variants from " << file1Path << "\n";
        return 1;
    }
    if (!loadVariants(file2Path, variantsFile2)) {
        std::cerr << "Error: Failed to load variants from " << file2Path << "\n";
        return 1;
    }

    // Identify differences
    std::vector<std::string> uniqueToFile1;
    std::vector<std::string> uniqueToFile2;

    for (const auto &v : variantsFile1) {
        if (variantsFile2.find(v) == variantsFile2.end()) {
            uniqueToFile1.push_back(v);
        }
    }
    for (const auto &v : variantsFile2) {
        if (variantsFile1.find(v) == variantsFile1.end()) {
            uniqueToFile2.push_back(v);
        }
    }

    // Print results
    std::cout << "Variants unique to " << file1Path << ":\n";
    for (auto &v : uniqueToFile1) {
        std::cout << v << "\n";
    }
    std::cout << "\nVariants unique to " << file2Path << ":\n";
    for (auto &v : uniqueToFile2) {
        std::cout << v << "\n";
    }

    return 0;
}

// ----------------------------------------------------------------------
// main
// ----------------------------------------------------------------------
int main(int argc, char* argv[]) {
    VCFXDiffTool diffTool;
    return diffTool.run(argc, argv);
}


./VCFX_diff_tool/VCFX_diff_tool.h
#ifndef VCFX_DIFF_TOOL_H
#define VCFX_DIFF_TOOL_H

#include <iostream>
#include <string>
#include <fstream>
#include <unordered_set>
#include <vector>

// VCFXDiffTool: Header file for VCF Diff Tool
class VCFXDiffTool {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Loads variants from a VCF file into a set
    bool loadVariants(const std::string& filePath, std::unordered_set<std::string>& variants);

    // Generates a unique key for a variant based on chromosome, position, ref, and alt
    std::string generateVariantKey(const std::string& chrom, const std::string& pos, const std::string& ref, const std::string& alt);
};

#endif // VCFX_DIFF_TOOL_H


./VCFX_distance_calculator/VCFX_distance_calculator.cpp
// VCFX_distance_calculator.cpp
#include "VCFX_distance_calculator.h"
#include <sstream>
#include <unordered_map>
#include <vector>
#include <limits>
#include <algorithm>
#include <cstdlib>

// --------------------------------------------------------------------------
// printHelp: Displays usage information.
// --------------------------------------------------------------------------
void printHelp() {
    std::cout << "VCFX_distance_calculator\n"
              << "Usage: VCFX_distance_calculator [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h               Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Calculates the distance between consecutive variants along each chromosome\n"
              << "  in a VCF file. Only the CHROM and POS columns are used.\n\n"
              << "Output (tab-delimited):\n"
              << "  CHROM   POS   PREV_POS   DISTANCE\n\n"
              << "Example:\n"
              << "  ./VCFX_distance_calculator < input.vcf > variant_distances.tsv\n";
}

// --------------------------------------------------------------------------
// parseVCFLine: Parses a VCF data line and extracts CHROM and POS.
// Returns false if the line is a header or cannot be parsed.
// --------------------------------------------------------------------------
bool parseVCFLine(const std::string& line, VCFVariant& variant) {
    // Skip header lines or empty lines.
    if (line.empty() || line[0] == '#')
        return false;

    std::stringstream ss(line);
    std::string chrom, pos_str;

    // We expect at least two tab-delimited columns: CHROM and POS.
    if (!std::getline(ss, chrom, '\t') ||
        !std::getline(ss, pos_str, '\t')) {
        return false;
    }

    try {
        variant.chrom = chrom;
        variant.pos = std::stoi(pos_str);
    } catch (...) {
        return false;
    }

    return true;
}

// --------------------------------------------------------------------------
// Structure to hold per-chromosome summary statistics.
// --------------------------------------------------------------------------
struct ChromStats {
    int count;             // Number of inter-variant distances computed
    long totalDistance;    // Sum of all distances
    int minDistance;       // Minimum distance seen
    int maxDistance;       // Maximum distance seen
    ChromStats() : count(0), totalDistance(0),
                   minDistance(std::numeric_limits<int>::max()),
                   maxDistance(0) {}
};

// --------------------------------------------------------------------------
// calculateDistances: Reads a VCF stream, calculates inter-variant distances,
// outputs a TSV line per variant, and writes summary statistics to stderr.
// --------------------------------------------------------------------------
bool calculateDistances(std::istream& in, std::ostream& out) {
    std::string line;
    bool headerFound = false;

    // Map to store the last variant position for each chromosome.
    std::unordered_map<std::string, int> lastPosMap;
    // Map to accumulate summary statistics per chromosome.
    std::unordered_map<std::string, ChromStats> chromStats;

    // Output TSV header.
    out << "CHROM\tPOS\tPREV_POS\tDISTANCE\n";

    while (std::getline(in, line)) {
        if (line.empty())
            continue;

        // Process header lines.
        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                headerFound = true;
            }
            continue;
        }

        if (!headerFound) {
            std::cerr << "Error: VCF header (#CHROM) not found before variant records.\n";
            return false;
        }

        VCFVariant variant;
        if (!parseVCFLine(line, variant)) {
            std::cerr << "Warning: Skipping malformed VCF line: " << line << "\n";
            continue;
        }

        // Check if we've seen a variant on this chromosome before.
        if (lastPosMap.find(variant.chrom) != lastPosMap.end()) {
            int prevPos = lastPosMap[variant.chrom];
            int distance = variant.pos - prevPos;
            out << variant.chrom << "\t" << variant.pos << "\t" << prevPos << "\t" << distance << "\n";

            // Update summary statistics.
            ChromStats &stats = chromStats[variant.chrom];
            stats.count++;
            stats.totalDistance += distance;
            stats.minDistance = std::min(stats.minDistance, distance);
            stats.maxDistance = std::max(stats.maxDistance, distance);
        } else {
            // No previous variant on this chromosome; output NA for distance.
            out << variant.chrom << "\t" << variant.pos << "\tNA\tNA\n";
        }

        // Update last position for this chromosome.
        lastPosMap[variant.chrom] = variant.pos;
    }

    // Output summary statistics to stderr.
    std::cerr << "\n=== Summary Statistics ===\n";
    for (const auto &entry : chromStats) {
        const std::string &chrom = entry.first;
        const ChromStats &stats = entry.second;
        double avgDistance = (stats.count > 0) ? static_cast<double>(stats.totalDistance) / stats.count : 0.0;
        std::cerr << "Chromosome: " << chrom << "\n"
                  << "  Variants compared: " << stats.count + 1 << "\n"
                  << "  Distances computed: " << stats.count << "\n"
                  << "  Total distance: " << stats.totalDistance << "\n"
                  << "  Min distance: " << stats.minDistance << "\n"
                  << "  Max distance: " << stats.maxDistance << "\n"
                  << "  Average distance: " << avgDistance << "\n\n";
    }

    return true;
}

// --------------------------------------------------------------------------
// main: Parses command-line arguments and calls calculateDistances.
// --------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // Check for help option.
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // Calculate distances from standard input to standard output.
    bool success = calculateDistances(std::cin, std::cout);
    return success ? 0 : 1;
}


./VCFX_distance_calculator/VCFX_distance_calculator.h
// VCFX_distance_calculator.h
#ifndef VCFX_DISTANCE_CALCULATOR_H
#define VCFX_DISTANCE_CALCULATOR_H

#include <iostream>
#include <string>

// Function to display help message
void printHelp();

// Structure to hold variant information (only CHROM and POS are used)
struct VCFVariant {
    std::string chrom;
    int pos;
};

// Function to parse a VCF line into a VCFVariant structure
bool parseVCFLine(const std::string& line, VCFVariant& variant);

// Function to calculate distances between consecutive variants
// Writes one output line per variant (with distance from the previous variant on that chromosome)
// and returns true on success.
bool calculateDistances(std::istream& in, std::ostream& out);

#endif // VCFX_DISTANCE_CALCULATOR_H


./VCFX_dosage_calculator/VCFX_dosage_calculator.cpp
#include "VCFX_dosage_calculator.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <cctype>

// ---------------------------------------------------------------------------
// displayHelp: Prints usage information to standard output.
// ---------------------------------------------------------------------------
void VCFXDosageCalculator::displayHelp() {
    std::cout << "VCFX_dosage_calculator: Calculate genotype dosage for each variant in a VCF file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_dosage_calculator [options] < input.vcf > dosage_output.txt\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Description:\n";
    std::cout << "  For each variant in the input VCF, the tool computes the dosage for each sample\n";
    std::cout << "  based on the genotype (GT) field. Dosage is defined as the number of alternate\n";
    std::cout << "  alleles (i.e. each allele > 0 counts as 1). Thus:\n";
    std::cout << "    0/0  => dosage 0\n";
    std::cout << "    0/1  => dosage 1\n";
    std::cout << "    1/1  => dosage 2\n";
    std::cout << "    1/2  => dosage 2  (each alternate, regardless of numeric value, counts as 1)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_dosage_calculator < input.vcf > dosage_output.txt\n";
}

// ---------------------------------------------------------------------------
// run: Parses command-line arguments and calls calculateDosage.
// ---------------------------------------------------------------------------
int VCFXDosageCalculator::run(int argc, char* argv[]) {
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {0,      0,           0,  0}
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
        return 0;
    }

    // Process VCF input from stdin and write dosage output to stdout.
    calculateDosage(std::cin, std::cout);

    return 0;
}

// ---------------------------------------------------------------------------
// calculateDosage: Processes the VCF line by line, computes dosage for each sample,
// and outputs a tab-delimited line per variant.
// ---------------------------------------------------------------------------
void VCFXDosageCalculator::calculateDosage(std::istream& in, std::ostream& out) {
    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerFields;

    // Print output header. The output columns are:
    // CHROM, POS, ID, REF, ALT, Dosages
    // where Dosages is a comma-separated list corresponding to each sample.
    out << "CHROM\tPOS\tID\tREF\tALT\tDosages\n";

    while (std::getline(in, line)) {
        if (line.empty())
            continue;

        // Process header lines
        if (line[0] == '#') {
            // Parse header line with "#CHROM" to get sample names.
            if (line.rfind("#CHROM", 0) == 0) {
                std::stringstream ss(line);
                std::string field;
                while (std::getline(ss, field, '\t')) {
                    headerFields.push_back(field);
                }
                headerParsed = true;
            }
            // Skip printing header lines in output.
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header (#CHROM) not found before variant records.\n";
            return;
        }

        // Split the variant line into fields (VCF standard: at least 10 fields expected)
        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string field;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        // We expect at least 10 columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, and at least one sample.
        if (fields.size() < 10) {
            std::cerr << "Warning: Skipping VCF line with fewer than 10 fields: " << line << "\n";
            continue;
        }

        // Extract variant-level fields.
        std::string chrom = fields[0];
        std::string pos = fields[1];
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        std::string format = fields[8];

        // For dosage calculation, we need to find the genotype (GT) field index.
        std::vector<std::string> formatFields = split(format, ':');
        int gtIndex = -1;
        for (size_t i = 0; i < formatFields.size(); ++i) {
            if (formatFields[i] == "GT") {
                gtIndex = static_cast<int>(i);
                break;
            }
        }
        if (gtIndex == -1) {
            std::cerr << "Warning: GT field not found in FORMAT column for variant " << chrom << ":" << pos << ".\n";
            out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\tNA\n";
            continue;
        }

        // Process each sample column (from index 9 onward)
        std::vector<std::string> dosageValues;
        for (size_t i = 9; i < fields.size(); ++i) {
            std::string sampleField = fields[i];
            // Split sample field by ':'; the first field is typically the genotype.
            std::vector<std::string> sampleParts = split(sampleField, ':');
            if (sampleParts.empty() || gtIndex >= static_cast<int>(sampleParts.size())) {
                dosageValues.push_back("NA");
                continue;
            }
            std::string genotype = sampleParts[gtIndex];
            // Normalize genotype separators: replace '|' with '/'
            std::replace(genotype.begin(), genotype.end(), '|', '/');

            if (genotype.empty() || genotype == ".") {
                dosageValues.push_back("NA");
                continue;
            }

            // Split genotype on '/'
            std::vector<std::string> alleles = split(genotype, '/');
            if (alleles.size() != 2) {
                // Non-diploid genotype or invalid format; mark as missing.
                dosageValues.push_back("NA");
                continue;
            }

            bool validGenotype = true;
            int dosage = 0;
            // For each allele, if it is "0", count 0; if it is any numeric >0, count 1.
            for (const std::string& a : alleles) {
                if (a == "." || a.empty()) {
                    validGenotype = false;
                    break;
                }
                // Ensure that the allele string is numeric.
                if (!std::all_of(a.begin(), a.end(), ::isdigit)) {
                    validGenotype = false;
                    break;
                }
                int alleleVal = std::stoi(a);
                // Count alternate alleles as 1 regardless of allele index.
                if (alleleVal > 0) {
                    dosage += 1;
                }
            }
            if (validGenotype) {
                dosageValues.push_back(std::to_string(dosage));
            } else {
                dosageValues.push_back("NA");
            }
        }

        // Join dosage values into a comma-separated string.
        std::ostringstream dosagesOSS;
        for (size_t i = 0; i < dosageValues.size(); ++i) {
            if (i > 0) {
                dosagesOSS << ",";
            }
            dosagesOSS << dosageValues[i];
        }

        // Output one line per variant with dosage information.
        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t"
            << dosagesOSS.str() << "\n";
    }
}

int main(int argc, char* argv[]) {
    VCFXDosageCalculator dosageCalculator;
    return dosageCalculator.run(argc, argv);
}


./VCFX_dosage_calculator/VCFX_dosage_calculator.h
#ifndef VCFX_DOSAGE_CALCULATOR_H
#define VCFX_DOSAGE_CALCULATOR_H

#include <iostream>
#include <string>
#include <vector>

// VCFX_dosage_calculator: Header file for Genotype Dosage Calculation tool
class VCFXDosageCalculator {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Calculates genotype dosage from VCF input and writes output
    void calculateDosage(std::istream& in, std::ostream& out);
};

#endif // VCFX_DOSAGE_CALCULATOR_H


./VCFX_duplicate_remover/VCFX_duplicate_remover.cpp
#include "VCFX_duplicate_remover.h"
#include <sstream>
#include <vector>
#include <algorithm>

// Utility function: Splits a string by a given delimiter.
static std::vector<std::string> splitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Function to display help message
void printHelp() {
    std::cout << "VCFX_duplicate_remover\n"
              << "Usage: VCFX_duplicate_remover [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h               Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Removes duplicate variants from a VCF file based on the combination of\n"
              << "  chromosome, position, REF, and ALT alleles. For multi-allelic records, the\n"
              << "  ALT field is normalized by sorting the comma-separated alleles so that the\n"
              << "  ordering does not affect duplicate detection.\n\n"
              << "Example:\n"
              << "  ./VCFX_duplicate_remover < input.vcf > unique_variants.vcf\n";
}

// Generates a unique key for a variant based on chrom, pos, ref, and alt.
// For multi-allelic ALT fields, the ALT alleles are split, sorted, and rejoined.
std::string generateNormalizedVariantKey(const std::string& chrom,
                                          const std::string& pos,
                                          const std::string& ref,
                                          const std::string& alt) {
    // Split ALT field on commas.
    std::vector<std::string> alts = splitString(alt, ',');
    // Sort alleles lexicographically.
    std::sort(alts.begin(), alts.end());
    // Rejoin sorted alleles into a single string.
    std::ostringstream oss;
    for (size_t i = 0; i < alts.size(); ++i) {
        if (i > 0) {
            oss << ",";
        }
        oss << alts[i];
    }
    std::string normalizedAlt = oss.str();
    return chrom + ":" + pos + ":" + ref + ":" + normalizedAlt;
}

// Helper function to generate a VariantKey from parsed values.
static VariantKey generateVariantKey(const std::string& chrom,
                                     const std::string& pos,
                                     const std::string& ref,
                                     const std::string& alt) {
    VariantKey key;
    key.chrom = chrom;
    try {
        key.pos = std::stoi(pos);
    } catch (...) {
        key.pos = 0;
    }
    key.ref = ref;
    key.alt = "";  // Will be set to normalized ALT.
    // Normalize ALT: sort multi-allelic entries.
    key.alt = generateNormalizedVariantKey(chrom, pos, ref, alt).substr(chrom.size() + pos.size() + ref.size() + 3); // skip prefix "chrom:pos:ref:"
    // Alternatively, simply:
    key.alt = generateNormalizedVariantKey(chrom, pos, ref, alt);
    // However, since generateNormalizedVariantKey already concatenates chrom:pos:ref:normalizedAlt,
    // we extract the normalizedAlt portion if needed. For simplicity, we can just store the full key.
    // For our VariantKey, we want: chrom, pos, ref, normalizedAlt.
    // We'll do that by re-parsing:
    std::vector<std::string> parts = splitString(generateNormalizedVariantKey(chrom, pos, ref, alt), ':');
    if (parts.size() >= 4) {
        key.chrom = parts[0];
        try {
            key.pos = std::stoi(parts[1]);
        } catch (...) {
            key.pos = 0;
        }
        key.ref = parts[2];
        key.alt = parts[3];
    }
    return key;
}

// Function to remove duplicate variants from a VCF file.
bool removeDuplicates(std::istream& in, std::ostream& out) {
    std::string line;
    // Print header lines as-is.
    // Use an unordered_set with our custom hash function to track seen variants.
    std::unordered_set<VariantKey, VariantKeyHash> seen_variants;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            // Print header lines unchanged.
            out << line << "\n";
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info;
        // Expect at least 8 columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info)) {
            std::cerr << "Warning: Skipping invalid VCF line.\n";
            continue;
        }

        // Generate a variant key that normalizes ALT ordering.
        VariantKey key = generateVariantKey(chrom, pos, ref, alt);
        if (seen_variants.find(key) == seen_variants.end()) {
            // Variant is unique; record and output the line.
            seen_variants.insert(key);
            out << line << "\n";
        } else {
            // Duplicate variant; skip.
        }
    }
    return true;
}

// ----------------------------------------------------------------------
// main: Parse command-line arguments and call removeDuplicates.
// ----------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // Simple argument parsing: if --help or -h is provided, print help.
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // Remove duplicates from standard input and write to standard output.
    bool success = removeDuplicates(std::cin, std::cout);
    return success ? 0 : 1;
}


./VCFX_duplicate_remover/VCFX_duplicate_remover.h
#ifndef VCFX_DUPLICATE_REMOVER_H
#define VCFX_DUPLICATE_REMOVER_H

#include <iostream>
#include <string>
#include <unordered_set>

// Function to display help message
void printHelp();

// Structure to represent a unique variant key
struct VariantKey {
    std::string chrom;
    int pos;
    std::string ref;
    std::string alt; // normalized: sorted, comma-separated alleles

    bool operator==(const VariantKey& other) const {
        return chrom == other.chrom &&
               pos == other.pos &&
               ref == other.ref &&
               alt == other.alt;
    }
};

// Custom hash function for VariantKey using a hash-combine approach
struct VariantKeyHash {
    std::size_t operator()(const VariantKey& k) const {
        std::size_t h1 = std::hash<std::string>()(k.chrom);
        std::size_t h2 = std::hash<int>()(k.pos);
        std::size_t h3 = std::hash<std::string>()(k.ref);
        std::size_t h4 = std::hash<std::string>()(k.alt);
        std::size_t seed = h1;
        seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h4 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

// Function to remove duplicate variants from a VCF file
bool removeDuplicates(std::istream& in, std::ostream& out);

#endif // VCFX_DUPLICATE_REMOVER_H


./VCFX_fasta_converter/VCFX_fasta_converter.cpp
#include "VCFX_fasta_converter.h"
#include <getopt.h>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <iomanip>
#include <cctype>
#include <map>

// A small map from two distinct bases to an IUPAC ambiguity code
// e.g. A + G => R, C + T => Y, etc.
static const std::map<std::string, char> IUPAC_ambiguities = {
    {"AG", 'R'}, {"GA", 'R'},
    {"CT", 'Y'}, {"TC", 'Y'},
    {"AC", 'M'}, {"CA", 'M'},
    {"GT", 'K'}, {"TG", 'K'},
    {"AT", 'W'}, {"TA", 'W'},
    {"CG", 'S'}, {"GC", 'S'}
};

// Utility to convert a single numeric allele index into the corresponding base
// returns '\0' on failure
static char alleleIndexToBase(int alleleIndex,
                              const std::string& ref,
                              const std::vector<std::string>& altAlleles)
{
    // 0 => ref, 1 => altAlleles[0], 2 => altAlleles[1], etc.
    if (alleleIndex == 0) {
        if (ref.size() == 1) {
            return std::toupper(ref[0]);
        } else {
            // multi-base or invalid for a single-locus representation
            return '\0';
        }
    } else {
        int altPos = alleleIndex - 1;
        if (altPos < 0 || (size_t)altPos >= altAlleles.size()) {
            return '\0'; // out of range
        }
        // altAlleles[altPos] must be a single base to be representable
        const std::string &a = altAlleles[altPos];
        if (a.size() == 1) {
            return std::toupper(a[0]);
        } else {
            // multi-base alt => can't represent as single base
            return '\0';
        }
    }
}

// If we have exactly two bases, see if there's a standard IUPAC code
// Otherwise returns 'N'
static char combineBasesIUPAC(char b1, char b2) {
    if (b1 == b2) {
        return b1;  // e.g. A + A => A
    }
    // build 2-char string in alphabetical order
    std::string pair;
    pair.push_back(std::min(b1, b2));
    pair.push_back(std::max(b1, b2));
    auto it = IUPAC_ambiguities.find(pair);
    if (it != IUPAC_ambiguities.end()) {
        return it->second;
    }
    return 'N'; // unknown combination
}

int VCFXFastaConverter::run(int argc, char* argv[]) {
    // Parse command-line arguments
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

    // Convert VCF input from stdin to FASTA output
    convertVCFtoFasta(std::cin, std::cout);
    return 0;
}

void VCFXFastaConverter::displayHelp() {
    std::cout << "VCFX_fasta_converter: Convert a variant-only VCF into simple per-sample FASTA.\n\n"
              << "Usage:\n"
              << "  VCFX_fasta_converter [options] < input.vcf > output.fasta\n\n"
              << "Description:\n"
              << "  Reads a VCF with diploid genotypes and writes a FASTA file. Each variant\n"
              << "  line becomes one position in the FASTA alignment. For multi-allelic sites,\n"
              << "  each sample's genotype is interpreted to produce a single IUPAC base\n"
              << "  (if heterozygous with different single-base alleles) or 'N' if ambiguous.\n\n"
              << "  Indels, multi-base alleles, or complicated genotypes default to 'N'.\n\n"
              << "Example:\n"
              << "  VCFX_fasta_converter < input.vcf > output.fasta\n\n";
}

void VCFXFastaConverter::convertVCFtoFasta(std::istream& in, std::ostream& out) {
    std::string line;
    std::vector<std::string> sampleNames;
    // Each sampleName -> sequence string
    std::unordered_map<std::string, std::string> sampleSequences;

    bool headerParsed = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            // Parse the #CHROM header to get sample columns
            if (line.rfind("#CHROM", 0) == 0) {
                std::stringstream ss(line);
                std::string field;
                // Skip the first 9 columns
                for (int i = 0; i < 9; ++i) {
                    if (!std::getline(ss, field, '\t')) {
                        break;
                    }
                }
                // Remaining fields are sample names
                while (std::getline(ss, field, '\t')) {
                    sampleNames.push_back(field);
                    sampleSequences[field] = "";
                }
                headerParsed = true;
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: #CHROM header not found before data lines.\n";
            return;
        }

        // Parse data line
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string fld;
            while (std::getline(ss, fld, '\t')) {
                fields.push_back(fld);
            }
        }
        // minimal check
        if (fields.size() < (9 + sampleNames.size())) {
            // not enough columns
            std::cerr << "Warning: Skipping malformed VCF line with insufficient columns.\n";
            continue;
        }

        // VCF standard columns
        //  0:CHROM, 1:POS, 2:ID, 3:REF, 4:ALT, 5:QUAL, 6:FILTER, 7:INFO, 8:FORMAT, 9+:samples
        const std::string &chrom = fields[0];
        // const std::string &pos = fields[1]; // not strictly needed for the FASTA
        // const std::string &id  = fields[2];
        const std::string &ref = fields[3];
        const std::string &altField = fields[4];
        const std::string &format   = fields[8];

        // Split alt on commas
        std::vector<std::string> altAlleles;
        {
            std::stringstream altSS(altField);
            std::string a;
            while (std::getline(altSS, a, ',')) {
                altAlleles.push_back(a);
            }
        }

        // Find GT index in format
        std::vector<std::string> formatFields;
        {
            std::stringstream fmts(format);
            std::string token;
            while (std::getline(fmts, token, ':')) {
                formatFields.push_back(token);
            }
        }
        int gtIndex = -1;
        for (size_t i = 0; i < formatFields.size(); ++i) {
            if (formatFields[i] == "GT") {
                gtIndex = static_cast<int>(i);
                break;
            }
        }
        bool hasGT = (gtIndex >= 0);

        // For each sample, figure out their genotype => one base or IUPAC or N
        for (size_t s = 0; s < sampleNames.size(); ++s) {
            const std::string &sampleName = sampleNames[s];
            // sample column is at fields[9 + s]
            if (9 + s >= fields.size()) {
                // missing sample column?
                sampleSequences[sampleName] += "N";
                continue;
            }
            const std::string &sampleData = fields[9 + s];
            if (!hasGT) {
                // no genotype => default to reference? or N?
                // We'll choose 'N' to avoid assumptions
                sampleSequences[sampleName] += "N";
                continue;
            }
            // parse sampleData by ':'
            std::vector<std::string> sampleParts;
            {
                std::stringstream sp(sampleData);
                std::string p;
                while (std::getline(sp, p, ':')) {
                    sampleParts.push_back(p);
                }
            }
            if (gtIndex >= (int)sampleParts.size()) {
                sampleSequences[sampleName] += "N";
                continue;
            }
            // genotype string
            std::string genotype = sampleParts[gtIndex];
            // unify separators
            for (char &c : genotype) {
                if (c == '|') c = '/';
            }
            if (genotype.empty() || genotype == ".") {
                sampleSequences[sampleName] += "N";
                continue;
            }
            // split genotype by '/'
            std::vector<std::string> alleles;
            {
                std::stringstream gtSS(genotype);
                std::string al;
                while (std::getline(gtSS, al, '/')) {
                    alleles.push_back(al);
                }
            }
            if (alleles.size() != 2) {
                // not diploid => 'N'
                sampleSequences[sampleName] += "N";
                continue;
            }

            // Convert each allele to a single base (char), or '\0' on fail
            char b1 = '\0';
            char b2 = '\0';
            {
                // if either is '.', skip
                if (alleles[0] == "." || alleles[1] == ".") {
                    sampleSequences[sampleName] += "N";
                    continue;
                }
                // parse numeric
                bool okA1 = true, okA2 = true;
                int a1 = 0, a2 = 0;
                // parse first allele
                {
                    for (char c : alleles[0]) {
                        if (!std::isdigit(c)) {okA1=false; break;}
                    }
                    if (okA1) a1 = std::stoi(alleles[0]);
                }
                // parse second allele
                {
                    for (char c : alleles[1]) {
                        if (!std::isdigit(c)) {okA2=false; break;}
                    }
                    if (okA2) a2 = std::stoi(alleles[1]);
                }
                if (!okA1 || !okA2) {
                    sampleSequences[sampleName] += "N";
                    continue;
                }
                b1 = alleleIndexToBase(a1, ref, altAlleles);
                b2 = alleleIndexToBase(a2, ref, altAlleles);
            }

            if (b1 == '\0' || b2 == '\0') {
                // means we couldn't interpret at least one allele
                sampleSequences[sampleName] += "N";
                continue;
            }

            // If same base => that base
            // If different => try IUPAC
            char finalBase = '\0';
            if (b1 == b2) {
                finalBase = b1;  // e.g. both 'A'
            } else {
                finalBase = combineBasesIUPAC(b1, b2); // might yield R, Y, etc., or 'N'
            }
            if (finalBase == '\0') finalBase = 'N';
            sampleSequences[sampleName] += finalBase;
        }
    }

    // Finally, output the sequences in FASTA format
    // e.g. >SampleName\n[sequence in 60-char lines]
    for (auto &kv : sampleSequences) {
        const std::string &sampleName = kv.first;
        const std::string &seq = kv.second;
        out << ">" << sampleName << "\n";
        // print in 60-char chunks
        for (size_t i = 0; i < seq.size(); i += 60) {
            out << seq.substr(i, 60) << "\n";
        }
    }
}


./VCFX_fasta_converter/VCFX_fasta_converter.h
#ifndef VCFX_FASTA_CONVERTER_H
#define VCFX_FASTA_CONVERTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXFastaConverter: Tool for converting a "variant-only" VCF into per-sample FASTA sequences
class VCFXFastaConverter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Converts VCF input to FASTA format
    void convertVCFtoFasta(std::istream& in, std::ostream& out);
};

#endif // VCFX_FASTA_CONVERTER_H


./VCFX_field_extractor/VCFX_field_extractor.cpp
#include "VCFX_field_extractor.h"
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <cctype>

// ------------------------------------------------------------------------
// printHelp
// ------------------------------------------------------------------------
void printHelp() {
    std::cout
        << "VCFX_field_extractor\n"
        << "Usage: VCFX_field_extractor --fields \"FIELD1,FIELD2,...\" [OPTIONS]\n\n"
        << "Description:\n"
        << "  Extracts specified fields from each VCF record. Fields can be:\n"
        << "    - Standard fields: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO\n"
        << "    - Subkeys in INFO (e.g. DP, AF, ANN). These are extracted from the INFO column.\n"
        << "    - Sample subfields: e.g. SampleName:GT or S2:DP, referencing the second sample's DP.\n"
        << "      You can use sample name as it appears in #CHROM line, or 'S' plus 1-based sample index.\n"
        << "If a requested field is not found or invalid, '.' is output.\n\n"
        << "Example:\n"
        << "  VCFX_field_extractor --fields \"CHROM,POS,ID,REF,ALT,DP,Sample1:GT\" < input.vcf > out.tsv\n\n"
        << "Options:\n"
        << "  --fields, -f   Comma-separated list of fields to extract\n"
        << "  --help, -h     Show this help message\n";
}

// ------------------------------------------------------------------------
// A utility function to parse "INFO" key-value pairs into a map.
// "INFO" might look like "DP=100;AF=0.5;ANN=some|stuff"
// ------------------------------------------------------------------------
static std::unordered_map<std::string, std::string> parseInfo(const std::string& infoField) {
    std::unordered_map<std::string, std::string> infoMap;
    if (infoField == "." || infoField.empty()) {
        return infoMap;
    }
    std::stringstream ss(infoField);
    std::string token;
    while (std::getline(ss, token, ';')) {
        if (token.empty()) {
            continue;
        }
        // We might have "KEY=VAL" or just "KEY" if it's a flag
        size_t eqPos = token.find('=');
        if (eqPos == std::string::npos) {
            // It's a flag
            infoMap[token] = "1";
        } else {
            std::string key = token.substr(0, eqPos);
            std::string val = token.substr(eqPos + 1);
            infoMap[key] = val;
        }
    }
    return infoMap;
}

// ------------------------------------------------------------------------
// parseLineExtract: Given a VCF data line (already split into columns),
// extracts the user-requested fields in 'fields' based on the logic:
//  - If field is one of CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO => direct from columns
//  - If it is an INFO subkey => parse info map
//  - If it is "SampleName:SUBFIELD" or "S<int>:SUBFIELD" => parse the genotype subfield
// ------------------------------------------------------------------------
static std::vector<std::string> parseLineExtract(
    const std::vector<std::string>& vcfCols,
    const std::vector<std::string>& fields,
    const std::unordered_map<std::string, int>& sampleNameToIndex)
{
    // The first 8 standard columns
    // 0:CHROM,1:POS,2:ID,3:REF,4:ALT,5:QUAL,6:FILTER,7:INFO,8:FORMAT,9+:samples
    std::vector<std::string> out;
    out.reserve(fields.size());

    // parse info subkeys
    std::unordered_map<std::string, std::string> infoMap;
    if (vcfCols.size() > 7) {
        infoMap = parseInfo(vcfCols[7]); // the INFO column
    }

    // parse the FORMAT column for sample subfields
    // e.g. if format= GT:DP:GQ, then subfield "DP" is index 1
    std::vector<std::string> formatTokens;
    if (vcfCols.size() > 8) {
        std::stringstream fmts(vcfCols[8]);
        std::string fmt;
        while (std::getline(fmts, fmt, ':')) {
            formatTokens.push_back(fmt);
        }
    }

    // For each requested field
    for (auto &fld : fields) {
        std::string value = "."; // default if not found

        // Check if it's a standard field
        if (fld == "CHROM") {
            if (vcfCols.size() > 0) value = vcfCols[0];
        } else if (fld == "POS") {
            if (vcfCols.size() > 1) value = vcfCols[1];
        } else if (fld == "ID") {
            if (vcfCols.size() > 2) value = vcfCols[2];
        } else if (fld == "REF") {
            if (vcfCols.size() > 3) value = vcfCols[3];
        } else if (fld == "ALT") {
            if (vcfCols.size() > 4) value = vcfCols[4];
        } else if (fld == "QUAL") {
            if (vcfCols.size() > 5) value = vcfCols[5];
        } else if (fld == "FILTER") {
            if (vcfCols.size() > 6) value = vcfCols[6];
        } else if (fld == "INFO") {
            if (vcfCols.size() > 7) value = vcfCols[7];
        } else {
            // Possibly an INFO subkey?
            if (infoMap.find(fld) != infoMap.end()) {
                value = infoMap[fld];
            } else {
                // Possibly a sample subfield: e.g. "SampleName:GT" or "S2:DP"
                // We'll parse something like "NAME:SUBFIELD"
                // or "S<index>:SUBFIELD"
                size_t colonPos = fld.find(':');
                if (colonPos != std::string::npos) {
                    std::string sampleNameOrID = fld.substr(0, colonPos);
                    std::string subfield = fld.substr(colonPos + 1);
                    // find sample index
                    int sampleColIndex = -1;
                    if (!sampleNameOrID.empty() && sampleNameOrID[0] == 'S'
                        && std::all_of(sampleNameOrID.begin()+1, sampleNameOrID.end(), ::isdigit))
                    {
                        // format S<int>
                        int idx = std::stoi(sampleNameOrID.substr(1));
                        // sample columns start at col=9 in the VCF, but idx is 1-based
                        sampleColIndex = 9 + (idx - 1);
                    } else {
                        // sample name
                        auto itS = sampleNameToIndex.find(sampleNameOrID);
                        if (itS != sampleNameToIndex.end()) {
                            sampleColIndex = itS->second;
                        }
                    }
                    // sampleColIndex is the VCF column with that sample
                    if (sampleColIndex >= 9 && (size_t)sampleColIndex < vcfCols.size()) {
                        // parse that sample field => split by ':'
                        std::vector<std::string> sampleTokens;
                        {
                            std::stringstream sss(vcfCols[sampleColIndex]);
                            std::string tkn;
                            while (std::getline(sss, tkn, ':')) {
                                sampleTokens.push_back(tkn);
                            }
                        }
                        // find subfield in the FORMAT
                        int subIx = -1;
                        for (int i=0; i<(int)formatTokens.size(); i++) {
                            if (formatTokens[i] == subfield) {
                                subIx = i;
                                break;
                            }
                        }
                        if (subIx >= 0 && subIx < (int)sampleTokens.size()) {
                            value = sampleTokens[subIx];
                        } else {
                            // subfield not found
                            value = ".";
                        }
                    } else {
                        // sample column not found
                        value = ".";
                    }
                } // end if colon found
            }
        } // end else standard field
        out.push_back(value);
    }

    return out;
}

// ------------------------------------------------------------------------
// extractFields: 
//   1) parse #CHROM line to identify sample names => build map sample->VCFcolumnIndex
//   2) for each data line, parse columns, parse info, parse sample subfields => output
// ------------------------------------------------------------------------
void extractFields(std::istream& in, std::ostream& out, const std::vector<std::string>& fields) {
    // Print the header row (the requested field names)
    for (size_t i = 0; i < fields.size(); i++) {
        out << fields[i];
        if (i + 1 < fields.size()) out << "\t";
    }
    out << "\n";

    std::string line;

    // We'll store sampleName -> columnIndex in this map
    std::unordered_map<std::string,int> sampleNameToIndex;
    bool foundChromHeader = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            // If it's the #CHROM line, parse out sample columns
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
                // parse columns
                std::stringstream ss(line);
                std::string col;
                int colIndex = 0;
                std::vector<std::string> hdrCols;
                while (std::getline(ss, col, '\t')) {
                    hdrCols.push_back(col);
                }
                // sample columns start at index=9
                for (int i = 9; i < (int)hdrCols.size(); i++) {
                    // sample name is hdrCols[i]
                    sampleNameToIndex[hdrCols[i]] = i;
                }
            }
            // skip printing header lines in the output
            continue;
        }

        // parse columns
        std::stringstream ss(line);
        std::vector<std::string> vcfCols;
        {
            std::string token;
            while (std::getline(ss, token, '\t')) {
                vcfCols.push_back(token);
            }
        }
        // parse the requested fields
        std::vector<std::string> extracted = parseLineExtract(vcfCols, fields, sampleNameToIndex);

        // print them as TSV
        for (size_t i=0; i<extracted.size(); i++) {
            out << extracted[i];
            if (i+1 < extracted.size()) out << "\t";
        }
        out << "\n";
    }
}

// ------------------------------------------------------------------------
// main: parse arguments, call extractFields
// ------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    std::vector<std::string> fields;
    bool showHelp = false;

    // Basic argument parsing
    for (int i=1; i<argc; i++) {
        std::string arg = argv[i];
        if (arg=="--help" || arg=="-h") {
            showHelp = true;
        } else if (arg.rfind("--fields",0)==0 || arg.rfind("-f",0)==0) {
            // parse next argument or after '='
            size_t eqPos = arg.find('=');
            if (eqPos != std::string::npos) {
                // e.g. --fields=CHROM,POS,ID
                std::string fieldsStr = arg.substr(eqPos+1);
                std::stringstream sss(fieldsStr);
                std::string f;
                while (std::getline(sss, f, ',')) {
                    fields.push_back(f);
                }
            } else {
                // next argument
                if (i+1<argc) {
                    i++;
                    std::string fieldsStr = argv[i];
                    std::stringstream sss(fieldsStr);
                    std::string f;
                    while (std::getline(sss, f, ',')) {
                        fields.push_back(f);
                    }
                }
            }
        }
    }

    if (showHelp) {
        printHelp();
        return 0;
    }
    if (fields.empty()) {
        std::cerr << "No fields specified. Use --fields or -f to specify.\n";
        std::cerr << "Use --help for usage.\n";
        return 1;
    }

    extractFields(std::cin, std::cout, fields);
    return 0;
}


./VCFX_field_extractor/VCFX_field_extractor.h
#ifndef VCFX_FIELD_EXTRACTOR_H
#define VCFX_FIELD_EXTRACTOR_H

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>

// Displays usage/help
void printHelp();

// Main extraction function that reads from 'in', writes to 'out', extracting
// user-specified fields (including standard fields, INFO subfields, and sample subfields).
void extractFields(std::istream& in, std::ostream& out, const std::vector<std::string>& fields);

#endif // VCFX_FIELD_EXTRACTOR_H


./VCFX_file_splitter/VCFX_file_splitter.cpp
#include "VCFX_file_splitter.h"
#include <getopt.h>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>

struct ChromFile {
    std::unique_ptr<std::ofstream> ofs;
    bool headerWritten;
};

int VCFXFileSplitter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string outputPrefix = "split";

    static struct option long_options[] = {
        {"help",   no_argument,       0, 'h'},
        {"prefix", required_argument, 0, 'p'},
        {0,        0,                 0,  0}
    };

    while ((opt = getopt_long(argc, argv, "hp:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'p':
                outputPrefix = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    splitVCFByChromosome(std::cin, outputPrefix);
    return 0;
}

void VCFXFileSplitter::displayHelp() {
    std::cout << "VCFX_file_splitter: Split a VCF file into multiple files based on chromosome.\n\n"
              << "Usage:\n"
              << "  VCFX_file_splitter [options] < input.vcf\n\n"
              << "Options:\n"
              << "  -h, --help            Display this help message and exit\n"
              << "  -p, --prefix <prefix> Output file prefix (default: 'split')\n\n"
              << "Example:\n"
              << "  VCFX_file_splitter --prefix \"chr\" < input.vcf\n";
}

// Splits the VCF by chromosome, writing the full header to each file.
void VCFXFileSplitter::splitVCFByChromosome(std::istream& in,
                                            const std::string& outputPrefix) {
    std::unordered_map<std::string, ChromFile> chromFiles;

    // We'll store all lines that begin with '#' (the header lines) until we hit
    // the first data line. We also handle the possibility of additional '#' lines
    // that appear after #CHROM, replicating them to all open chromosome files.
    std::vector<std::string> initialHeaderLines;
    bool foundFirstDataLine = false;

    // We'll read the entire file line by line, but we do two main phases:
    //   1) Collect all '#' lines (the main header). Once we see a data line,
    //      we treat it as phase 2, splitting lines by chromosome.
    //   2) If we encounter extra '#' lines in phase 2, we replicate them
    //      to all already-open files.

    std::string line;
    while (true) {
        std::streampos currentPos = in.tellg();
        if (!std::getline(in, line)) {
            // end of file
            break;
        }
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (!foundFirstDataLine) {
                // We haven't yet seen any data lines, so this is part of the "initial" header
                initialHeaderLines.push_back(line);
            } else {
                // This is an extra header line after data lines started, replicate to open files
                for (auto &kv : chromFiles) {
                    ChromFile &cf = kv.second;
                    if (cf.ofs && cf.ofs->is_open()) {
                        *(cf.ofs) << line << "\n";
                    }
                }
            }
        } else {
            // This is a data line
            if (!foundFirstDataLine) {
                // We just encountered the first data line
                // so from now on, we are in the "splitting" phase
                foundFirstDataLine = true;
            }
            // We parse the chromosome from the data line
            std::stringstream ss(line);
            std::string chrom;
            if (!std::getline(ss, chrom, '\t')) {
                std::cerr << "Warning: cannot parse CHROM from line: " << line << "\n";
                continue;
            }
            // Check or create file
            if (chromFiles.find(chrom) == chromFiles.end()) {
                // Create a new file
                std::string filename = outputPrefix + "_" + chrom + ".vcf";
                ChromFile cf;
                cf.ofs = std::make_unique<std::ofstream>(filename);
                if (!cf.ofs->is_open()) {
                    std::cerr << "Error: Unable to create file: " << filename << "\n";
                    continue;
                }
                cf.headerWritten = false;
                chromFiles[chrom] = std::move(cf);
            }
            ChromFile &cf = chromFiles[chrom];
            if (!cf.headerWritten) {
                // Write all initial header lines
                for (auto &hLine : initialHeaderLines) {
                    *(cf.ofs) << hLine << "\n";
                }
                cf.headerWritten = true;
            }
            // Write this data line
            *(cf.ofs) << line << "\n";
        }
    }

    // Close all files
    for (auto &kv : chromFiles) {
        ChromFile &cf = kv.second;
        if (cf.ofs && cf.ofs->is_open()) {
            cf.ofs->close();
        }
    }

    // If we never found any data line, it's possible the input was all header or empty
    // In that case, we do nothing or we could emit a warning
    if (!foundFirstDataLine) {
        // Possibly we can check if initialHeaderLines is non-empty => maybe we do something
        // But there's no data to split
        std::cerr << "Note: No variant data lines were found in the input.\n";
    }
}

int main(int argc, char* argv[]) {
    VCFXFileSplitter splitter;
    return splitter.run(argc, argv);
}


./VCFX_file_splitter/VCFX_file_splitter.h
#ifndef VCFX_FILE_SPLITTER_H
#define VCFX_FILE_SPLITTER_H

#include <iostream>
#include <string>

// VCFXFileSplitter: Splits a VCF file by chromosome into multiple smaller VCFs.
class VCFXFileSplitter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Splits the input VCF by chromosome using the specified prefix
    void splitVCFByChromosome(std::istream& in, const std::string& outputPrefix);
};

#endif // VCFX_FILE_SPLITTER_H


