#include "VCFX_phase_checker.h"
#include <getopt.h>

// Implementation of VCFX_phase_checker
int VCFXPhaseChecker::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;

    while ((opt = getopt(argc, argv, "h")) != -1) {
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

    // Process VCF input from stdin
    processVCF(std::cin);

    return 0;
}

void VCFXPhaseChecker::displayHelp() {
    std::cout << "VCFX_phase_checker: Check if variants are phased in a VCF file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_phase_checker [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_phase_checker < input.vcf\n";
}

void VCFXPhaseChecker::processVCF(std::istream& in) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            std::cout << line << "\n"; // Preserve header lines
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

        if (fields.size() < 10) {
            std::cerr << "Invalid VCF line with fewer than 10 fields.\n";
            continue;
        }

        // GT is typically in FORMAT field followed by sample fields
        std::string format = fields[8];
        std::vector<std::string> format_fields;
        pos = 0;
        while ((pos = format.find(':')) != std::string::npos) {
            field = format.substr(0, pos);
            format_fields.push_back(field);
            format.erase(0, pos + 1);
        }
        format_fields.push_back(format);

        // Find the index of the GT field
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

        // Iterate over sample columns
        bool all_phased = true;
        for (size_t i = 9; i < fields.size(); ++i) {
            std::vector<std::string> sample_fields;
            std::string sample = fields[i];
            size_t s_pos = 0;
            while ((s_pos = sample.find(':')) != std::string::npos) {
                field = sample.substr(0, s_pos);
                sample_fields.push_back(field);
                sample.erase(0, s_pos + 1);
            }
            sample_fields.push_back(sample);

            if (gt_index >= static_cast<int>(sample_fields.size())) {
                std::cerr << "GT index out of range in sample fields.\n";
                continue;
            }

            std::string genotype = sample_fields[gt_index];
            if (!isPhased(genotype)) {
                all_phased = false;
                break;
            }
        }

        if (all_phased) {
            std::cout << line << "\n";
        } else {
            std::cerr << "Unphased genotype found at position " << fields[1] << "\n";
        }
    }
}

bool VCFXPhaseChecker::isPhased(const std::string& genotype) {
    return genotype.find('|') != std::string::npos;
}

int main(int argc, char* argv[]) {
    VCFXPhaseChecker checker;
    return checker.run(argc, argv);
}

