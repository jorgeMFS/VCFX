#include "VCFX_ref_comparator.h"
#include <getopt.h>
#include <sstream>
#include <fstream>
#include <algorithm>

// Implementation of VCFXRefComparator
int VCFXRefComparator::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string referencePath;

    static struct option long_options[] = {
        {"help",        no_argument,       0, 'h'},
        {"reference",   required_argument, 0, 'r'},
        {0,             0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hr:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'r':
                referencePath = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || referencePath.empty()) {
        displayHelp();
        return 1;
    }

    // Load the reference genome
    if (!loadReference(referencePath)) {
        std::cerr << "Error: Failed to load reference genome from " << referencePath << "\n";
        return 1;
    }

    // Compare VCF variants against the reference genome
    compareWithReference(std::cin, std::cout);

    return 0;
}

void VCFXRefComparator::displayHelp() {
    std::cout << "VCFX_ref_comparator: Compare VCF variants against a reference genome.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_ref_comparator --reference <reference.fasta> [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                  Display this help message and exit\n";
    std::cout << "  -r, --reference <file>      Path to the reference genome FASTA file\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_ref_comparator --reference reference.fasta < input.vcf > comparison_output.vcf\n";
}

bool VCFXRefComparator::loadReference(const std::string& referencePath) {
    std::ifstream refFile(referencePath);
    if (!refFile.is_open()) {
        std::cerr << "Error: Cannot open reference genome file: " << referencePath << "\n";
        return false;
    }

    std::string line;
    std::string currentChrom;
    std::ostringstream sequenceStream;

    while (std::getline(refFile, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            // Save the previous chromosome sequence
            if (!currentChrom.empty()) {
                referenceGenome[currentChrom] = sequenceStream.str();
                sequenceStream.str("");
                sequenceStream.clear();
            }
            // Extract chromosome name
            std::stringstream ss(line.substr(1));
            ss >> currentChrom;
        } else {
            // Append sequence lines, removing any whitespace and converting to uppercase
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            sequenceStream << line;
        }
    }

    // Save the last chromosome sequence
    if (!currentChrom.empty()) {
        referenceGenome[currentChrom] = sequenceStream.str();
    }

    refFile.close();
    return true;
}

void VCFXRefComparator::compareWithReference(std::istream& vcfInput, std::ostream& vcfOutput) {
    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerFields;
    int infoIndex = -1;

    while (std::getline(vcfInput, line)) {
        if (line.empty()) {
            vcfOutput << "\n";
            continue;
        }

        if (line[0] == '#') {
            if (line.substr(0, 6) == "#CHROM") {
                // Parse header line
                std::stringstream ss(line);
                std::string field;
                headerFields.clear();
                while (std::getline(ss, field, '\t')) {
                    headerFields.push_back(field);
                }

                // Add new INFO field for comparison result
                std::cout << "##INFO=<ID=REF_COMPARISON,Number=1,Type=String,Description=\"Comparison of variant alleles against the reference genome\">\n";
                
                // Write modified header
                std::stringstream newHeader;
                for (size_t i = 0; i < headerFields.size(); ++i) {
                    newHeader << headerFields[i];
                    if (i != headerFields.size() - 1) {
                        newHeader << "\t";
                    }
                }
                newHeader << "\tREF_COMPARISON\n";
                vcfOutput << newHeader.str() << "\n";
                headerParsed = true;
            } else {
                // Write other header lines as-is
                vcfOutput << line << "\n";
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header line with #CHROM not found.\n";
            return;
        }

        // Parse VCF data line
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

        std::string chrom = fieldsVec[0];
        int pos = std::stoi(fieldsVec[1]);
        std::string id = fieldsVec[2];
        std::string ref = fieldsVec[3];
        std::string alt = fieldsVec[4];
        std::string qual = fieldsVec[5];
        std::string filter = fieldsVec[6];
        std::string info = fieldsVec[7];

        // Retrieve reference allele from the reference genome
        if (referenceGenome.find(chrom) == referenceGenome.end()) {
            std::cerr << "Warning: Chromosome " << chrom << " not found in reference genome.\n";
            info += ";REF_COMPARISON=UNKNOWN_CHROM";
            vcfOutput << line << "\t" << "UNKNOWN_CHROM" << "\n";
            continue;
        }

        const std::string& refSeq = referenceGenome[chrom];
        if (pos < 1 || pos > static_cast<int>(refSeq.size())) {
            std::cerr << "Warning: Position " << pos << " out of bounds for chromosome " << chrom << ".\n";
            info += ";REF_COMPARISON=INVALID_POS";
            vcfOutput << line << "\t" << "INVALID_POS" << "\n";
            continue;
        }

        // Extract reference allele from the reference genome
        std::string refFromGenome = refSeq.substr(pos - 1, ref.size());

        // Compare VCF ref allele with reference genome
        bool refMatch = (ref == refFromGenome);

        // Compare alternate alleles with reference genome
        std::vector<std::string> altAlleles;
        std::stringstream alt_ss(alt);
        std::string alt_allele;
        while (std::getline(alt_ss, alt_allele, ',')) {
            altAlleles.push_back(alt_allele);
        }

        std::vector<std::string> comparisonResults;
        for (const auto& altA : altAlleles) {
            if (altA == refFromGenome) {
                comparisonResults.push_back("REF_MATCH");
            } else {
                comparisonResults.push_back("NOVEL");
            }
        }

        // Prepare comparison result string
        std::string comparisonStr;
        for (size_t i = 0; i < comparisonResults.size(); ++i) {
            comparisonStr += comparisonResults[i];
            if (i != comparisonResults.size() - 1) {
                comparisonStr += ",";
            }
        }

        // Append to INFO field
        if (info.back() != ';') {
            info += ";";
        }
        info += "REF_COMPARISON=" + comparisonStr;

        // Reconstruct VCF line with updated INFO field
        std::stringstream newLine;
        for (size_t i = 0; i < 8; ++i) {
            newLine << fieldsVec[i] << "\t";
        }
        newLine << info << "\t";

        // Append any additional fields (e.g., FORMAT and samples)
        for (size_t i = 8; i < fieldsVec.size(); ++i) {
            newLine << fieldsVec[i];
            if (i != fieldsVec.size() - 1) {
                newLine << "\t";
            }
        }

        // Append the comparison result
        newLine << "\t" << comparisonStr;

        vcfOutput << newLine.str() << "\n";
    }
}

int main(int argc, char* argv[]) {
    VCFXRefComparator refComparator;
    return refComparator.run(argc, argv);
}