#include "VCFX_dosage_calculator.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <iomanip>

// Implementation of VCFX_dosage_calculator
int VCFXDosageCalculator::run(int argc, char* argv[]) {
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

    // Calculate genotype dosage from stdin and output to stdout
    calculateDosage(std::cin, std::cout);

    return 0;
}

void VCFXDosageCalculator::displayHelp() {
    std::cout << "VCFX_dosage_calculator: Calculate genotype dosage for each variant.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_dosage_calculator [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_dosage_calculator < input.vcf > dosage_output.txt\n";
}

void VCFXDosageCalculator::calculateDosage(std::istream& in, std::ostream& out) {
    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerFields;

    // Print header for dosage output
    out << "CHROM\tPOS\tID\tREF\tALT\tDosages\n";

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            // Parse header lines
            if (line.substr(0, 6) == "#CHROM") {
                std::stringstream ss(line);
                std::string field;
                // Collect header fields
                while (std::getline(ss, field, '\t')) {
                    headerFields.push_back(field);
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
        std::string qual = fieldsVec[5];
        std::string filter = fieldsVec[6];
        std::string info = fieldsVec[7];
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
            // Append missing dosage if GT is not found
            out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\tNA\n";
            continue;
        }

        // Iterate over each sample to calculate dosage
        std::vector<int> dosages;
        bool valid = true;

        for (size_t i = 9; i < fieldsVec.size(); ++i) {
            std::string sample = fieldsVec[i];
            std::vector<std::string> sampleFields;
            std::stringstream samp_ss(sample);
            std::string samp_field;
            while (std::getline(samp_ss, samp_field, ':')) {
                sampleFields.push_back(samp_field);
            }

            if (gt_index >= static_cast<int>(sampleFields.size())) {
                std::cerr << "GT index out of range in sample fields.\n";
                dosages.push_back(-1); // Indicate missing dosage
                continue;
            }

            std::string genotype = sampleFields[gt_index];
            // Replace '|' with '/' for consistency
            std::replace(genotype.begin(), genotype.end(), '|', '/');

            if (genotype == "." || genotype.empty()) {
                dosages.push_back(-1); // Missing genotype
                continue;
            }

            std::vector<std::string> alleles;
            std::stringstream genotype_ss(genotype);
            std::string allele;
            while (std::getline(genotype_ss, allele, '/')) {
                alleles.push_back(allele);
            }

            if (alleles.size() != 2) {
                dosages.push_back(-1); // Non-diploid genotype
                continue;
            }

            int dosage = 0;
            bool genotypeValid = true;

            for (const auto& a : alleles) {
                if (a == ".") {
                    genotypeValid = false;
                    break;
                }
                try {
                    int alleleNum = std::stoi(a);
                    dosage += alleleNum;
                } catch (...) {
                    genotypeValid = false;
                    break;
                }
            }

            if (genotypeValid) {
                dosages.push_back(dosage);
            } else {
                dosages.push_back(-1); // Invalid genotype
            }
        }

        // Prepare dosage string
        std::ostringstream dosage_ss;
        for (size_t i = 0; i < dosages.size(); ++i) {
            if (dosages[i] == -1) {
                dosage_ss << "NA";
            } else {
                dosage_ss << dosages[i];
            }
            if (i != dosages.size() - 1) {
                dosage_ss << ",";
            }
        }

        // Output the dosage information
        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" 
            << dosage_ss.str() << "\n";
    }
}
