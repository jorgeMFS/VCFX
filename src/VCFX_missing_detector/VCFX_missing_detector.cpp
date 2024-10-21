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
