#include "VCFX_alignment_checker.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>

// Implementation

int VCFXAlignmentChecker::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string vcfFile = "";
    std::string refFile = "";

    static struct option long_options[] = {
        {"help",                 no_argument,       0, 'h'},
        {"alignment-discrepancy", no_argument,       0, 'a'},
        {0,                      0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "ha", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'a':
                // Alignment discrepancy mode
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    // Ensure that reference genome is provided via stdin after VCF
    // Usage: VCFX_alignment_checker --alignment-discrepancy < input.vcf < reference.fasta

    // Open VCF input (first stdin)
    std::istream& vcfIn = std::cin;


    if (optind < argc) {
        refFile = std::string(argv[optind]);
    } else {
        // If no file provided, read from second stdin
        // This is not straightforward; inform the user to provide reference genome file
        std::cerr << "Error: Reference genome file must be provided as an argument.\n";
        displayHelp();
        return 1;
    }

    // Open reference genome file
    std::ifstream refStream(refFile);
    if (!refStream.is_open()) {
        std::cerr << "Error: Unable to open reference genome file: " << refFile << "\n";
        return 1;
    }

    // Load reference genome
    if (!loadReferenceGenome(refStream)) {
        std::cerr << "Error: Failed to load reference genome.\n";
        return 1;
    }

    // Check discrepancies
    checkDiscrepancies(vcfIn, refStream, std::cout);

    return 0;
}

void VCFXAlignmentChecker::displayHelp() {
    std::cout << "VCFX_alignment_checker: Identify discrepancies between VCF variants and a reference genome.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_alignment_checker --alignment-discrepancy <vcf_file> <reference.fasta>\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                        Display this help message and exit\n";
    std::cout << "  -a, --alignment-discrepancy        Identify alignment discrepancies\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_alignment_checker --alignment-discrepancy input.vcf reference.fasta > discrepancies.txt\n";
}

bool VCFXAlignmentChecker::loadReferenceGenome(std::istream& in) {
    std::string line;
    std::string currentChrom = "";
    std::string seq = "";

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!currentChrom.empty()) {
                // Store the previous chromosome sequence
                referenceGenome[normalizeChromosome(currentChrom)] = seq;
                seq.clear();
            }
            // Get chromosome name
            size_t pos = line.find(' ');
            if (pos != std::string::npos) {
                currentChrom = line.substr(1, pos - 1);
            } else {
                currentChrom = line.substr(1);
            }
        } else {
            // Append sequence lines
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            seq += line;
        }
    }

    // Store last chromosome
    if (!currentChrom.empty()) {
        referenceGenome[normalizeChromosome(currentChrom)] = seq;
    }

    return true;
}

std::string VCFXAlignmentChecker::normalizeChromosome(const std::string& chrom) {
    std::string norm = chrom;
    if (norm.find("chr") != 0 && 
        !(norm == "X" || norm == "Y" || norm == "MT" || 
          (std::all_of(norm.begin(), norm.end(), ::isdigit)))) {
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
    if (pos < 1 || static_cast<size_t>(pos -1 + length) > seq.length()) {
        return "";
    }

    return seq.substr(pos -1, length);
}

void VCFXAlignmentChecker::checkDiscrepancies(std::istream& vcfIn, std::istream& refIn, std::ostream& out) {
    std::string line;
    bool headerParsed = false;
    int chrIndex = -1, posIndex = -1, refIndex = -1, altIndex = -1;

    // Output header for discrepancies
    out << "CHROM\tPOS\tID\tREF\tALT\tDiscrepancy_Type\tReference_Value\tVCF_Value\n";

    while (std::getline(vcfIn, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.substr(0, 6) == "#CHROM") {
                // Parse header to get column indices
                std::stringstream ss(line);
                std::string field;
                std::vector<std::string> headers;
                while (std::getline(ss, field, '\t')) {
                    headers.push_back(field);
                }
                for (size_t i = 0; i < headers.size(); ++i) {
                    if (headers[i] == "CHROM") chrIndex = static_cast<int>(i);
                    else if (headers[i] == "POS") posIndex = static_cast<int>(i);
                    else if (headers[i] == "REF") refIndex = static_cast<int>(i);
                    else if (headers[i] == "ALT") altIndex = static_cast<int>(i);
                }
                if (chrIndex == -1 || posIndex == -1 || refIndex == -1 || altIndex == -1) {
                    std::cerr << "Error: VCF header does not contain required fields.\n";
                    return;
                }
                headerParsed = true;
            }
            // Continue to next line
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header line with #CHROM not found.\n";
            return;
        }

        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string field;

        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() < 5) {
            std::cerr << "Warning: Skipping invalid VCF line (less than 5 fields): " << line << "\n";
            continue;
        }

        std::string chrom = fields[chrIndex];
        int pos = 0;
        try {
            pos = std::stoi(fields[posIndex]);
        } catch (...) {
            std::cerr << "Warning: Invalid POS value. Skipping line: " << line << "\n";
            continue;
        }

        std::string ref = fields[refIndex];
        std::string alt = fields[altIndex];

        // Handle multi-allelic ALT fields
        std::vector<std::string> alts;
        std::stringstream altSS(alt);
        while (std::getline(altSS, field, ',')) {
            alts.push_back(field);
        }

        for (const auto& allele : alts) {
            if (ref.length() == 1 && allele.length() == 1) {
                // SNP
                std::string ref_base = getReferenceBases(chrom, pos, 1);
                if (ref_base.empty()) {
                    std::cerr << "Warning: Reference base not found for " << chrom << ":" << pos << "\n";
                    continue;
                }

                if (ref != ref_base) {
                    out << chrom << "\t" << pos << "\t" << fields[2] << "\t" << ref 
                        << "\t" << allele << "\t" << "REF_MISMATCH" 
                        << "\t" << ref_base << "\t" << ref << "\n";
                }

                std::string alt_base = getReferenceBases(chrom, pos, 1);
                if (alt_base.empty()) {
                    std::cerr << "Warning: Reference base not found for " << chrom << ":" << pos << "\n";
                    continue;
                }

                if (allele != alt_base) {
                    out << chrom << "\t" << pos << "\t" << fields[2] << "\t" << ref 
                        << "\t" << allele << "\t" << "ALT_MISMATCH" 
                        << "\t" << alt_base << "\t" << allele << "\n";
                }

            } else {
                // Indels or complex variants
                size_t ref_len = ref.length();
                size_t alt_len = allele.length();
                size_t len = std::min(ref_len, alt_len);

                std::string ref_seq = getReferenceBases(chrom, pos, len);
                if (ref_seq.empty()) {
                    std::cerr << "Warning: Reference sequence not found for " << chrom << ":" << pos << "\n";
                    continue;
                }

                std::string vcf_ref = ref.substr(0, len);
                std::string vcf_alt = allele.substr(0, len);

                if (vcf_ref != ref_seq) {
                    out << chrom << "\t" << pos << "\t" << fields[2] << "\t" << ref 
                        << "\t" << allele << "\t" << "REF_DISCREPANCY" 
                        << "\t" << ref_seq << "\t" << vcf_ref << "\n";
                }

                if (vcf_alt != ref_seq) {
                    out << chrom << "\t" << pos << "\t" << fields[2] << "\t" << ref 
                        << "\t" << allele << "\t" << "ALT_DISCREPANCY" 
                        << "\t" << ref_seq << "\t" << vcf_alt << "\n";
                }
            }
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXAlignmentChecker alignmentChecker;
    return alignmentChecker.run(argc, argv);
}