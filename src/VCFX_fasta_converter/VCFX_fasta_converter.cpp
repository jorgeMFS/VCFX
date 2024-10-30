#include "VCFX_fasta_converter.h"
#include <getopt.h>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <iomanip>

// Implementation of VCFXFastaConverter
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
    std::cout << "VCFX_fasta_converter: Convert VCF data into FASTA format.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_fasta_converter [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_fasta_converter < input.vcf > output.fasta\n";
}

void VCFXFastaConverter::convertVCFtoFasta(std::istream& in, std::ostream& out) {
    std::string line;
    std::vector<std::string> sampleNames;
    std::unordered_map<std::string, std::string> sampleSequences;
    bool headerParsed = false;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            // Parse header lines
            if (line.substr(0, 6) == "#CHROM") {
                std::stringstream ss(line);
                std::string field;
                // Skip the first 9 columns
                for (int i = 0; i < 9; ++i) {
                    std::getline(ss, field, '\t');
                }
                // The rest are sample names
                while (std::getline(ss, field, '\t')) {
                    sampleNames.push_back(field);
                    sampleSequences[field] = "";
                }
                headerParsed = true;
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "VCF header line with #CHROM not found.\n";
            return;
        }

        // Parse VCF data lines
        std::stringstream ss(line);
        std::string field;
        std::vector<std::string> fieldsVec;
        while (std::getline(ss, field, '\t')) {
            fieldsVec.push_back(field);
        }

        if (fieldsVec.size() < 10) {
            std::cerr << "Invalid VCF line with fewer than 10 fields.\n";
            continue;
        }

        std::string chrom = fieldsVec[0];
        std::string pos = fieldsVec[1];
        std::string id = fieldsVec[2];
        std::string ref = fieldsVec[3];
        std::string alt = fieldsVec[4];
        std::string format = fieldsVec[8];

        // Split FORMAT field to find the index of genotype (GT)
        std::vector<std::string> formatFields;
        std::stringstream fmt_ss(format);
        std::string fmt_field;
        while (std::getline(fmt_ss, fmt_field, ':')) {
            formatFields.push_back(fmt_field);
        }

        int gt_index = -1;
        for (size_t i = 0; i < formatFields.size(); ++i) {
            if (formatFields[i] == "GT") {
                gt_index = static_cast<int>(i);
                break;
            }
        }

        if (gt_index == -1) {
            std::cerr << "GT field not found in FORMAT column.\n";
            // Append reference allele if GT is not found
            for (const auto& sample : sampleNames) {
                sampleSequences[sample] += ref;
            }
            continue;
        }

        // Iterate over each sample and append nucleotide based on genotype
        for (size_t i = 0; i < sampleNames.size(); ++i) {
            std::string sample = fieldsVec[9 + i];
            std::vector<std::string> sampleFields;
            std::stringstream samp_ss(sample);
            std::string samp_field;
            while (std::getline(samp_ss, samp_field, ':')) {
                sampleFields.push_back(samp_field);
            }
            if (gt_index >= static_cast<int>(sampleFields.size())) {
                std::cerr << "GT index out of range in sample fields.\n";
                sampleSequences[sampleNames[i]] += "N"; // Unknown nucleotide
                continue;
            }
            std::string genotype = sampleFields[gt_index];
            // Replace '|' with '/' for consistency
            std::replace(genotype.begin(), genotype.end(), '|', '/');
            std::vector<std::string> alleles;
            std::stringstream genotype_ss(genotype);
            std::string allele;
            while (std::getline(genotype_ss, allele, '/')) {
                alleles.push_back(allele);
            }
            if (alleles.size() != 2) {
                std::cerr << "Non-diploid genotype encountered.\n";
                sampleSequences[sampleNames[i]] += "N"; // Unknown nucleotide
                continue;
            }

            // Determine nucleotide to append based on genotype
            // 0 -> reference allele, 1 -> first alternate allele, etc.
            std::string nucleotide = "N"; // Default unknown
            try {
                if (alleles[0] == "." || alleles[1] == ".") {
                    nucleotide = "N"; // Missing genotype
                } else {
                    int allele1 = std::stoi(alleles[0]);
                    int allele2 = std::stoi(alleles[1]);
                    if (allele1 == 0 && allele2 == 0) {
                        nucleotide = ref;
                    } else if ((allele1 == 0 && allele2 == 1) ||
                               (allele1 == 1 && allele2 == 0) ||
                               (allele1 == 1 && allele2 == 1)) {
                        // For homozygous or heterozygous alternate alleles, take first alternate allele
                        nucleotide = (alt.find(',') != std::string::npos) ? alt.substr(0, 1) : alt;
                    } else {
                        nucleotide = "N"; // For multi-allelic or unexpected
                    }
                }
            } catch (...) {
                nucleotide = "N"; // In case of conversion failure
            }
            sampleSequences[sampleNames[i]] += nucleotide;
        }
    }

    // Output FASTA sequences
    for (const auto& sample : sampleNames) {
        out << ">" << sample << "\n";
        std::string seq = sampleSequences[sample];
        // Print sequence in lines of 60 characters
        for (size_t i = 0; i < seq.length(); i += 60) {
            out << seq.substr(i, 60) << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXFastaConverter converter;
    return converter.run(argc, argv);
}
