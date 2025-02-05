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
