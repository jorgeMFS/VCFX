#include "VCFX_phred_filter.h"
#include <getopt.h>
#include <sstream>

// Implementation of VCFX_phred_filter
int VCFXPhredFilter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    double threshold = 30.0; // Default threshold
    bool showHelp = false;

    static struct option long_options[] = {
        {"phred-filter", required_argument, 0, 'p'},
        {"help",         no_argument,       0, 'h'},
        {0,              0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "p:h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'p':
                {
                    std::stringstream ss(optarg);
                    if (!(ss >> threshold)) {
                        std::cerr << "Invalid threshold value: " << optarg << "\n";
                        return 1;
                    }
                }
                break;
            case 'h':
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Process VCF input from stdin
    processVCF(std::cin, threshold);

    return 0;
}

void VCFXPhredFilter::displayHelp() {
    std::cout << "VCFX_phred_filter: Filter variants based on Phred quality score.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_phred_filter --phred-filter <threshold> [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -p, --phred-filter   Phred quality score threshold (e.g., 30)\n";
    std::cout << "  -h, --help           Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_phred_filter --phred-filter 30 < input.vcf > filtered.vcf\n";
}

void VCFXPhredFilter::processVCF(std::istream& in, double threshold) {
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

        if (fields.size() < 5) {
            std::cerr << "Invalid VCF line with fewer than 5 fields.\n";
            continue;
        }

        std::string qualStr = fields[5];
        double qual = parseQUAL(qualStr);

        if (qual >= threshold) {
            std::cout << line << "\n";
        }
    }
}

double VCFXPhredFilter::parseQUAL(const std::string& qualStr) {
    if (qualStr == ".") {
        return 0.0; // Treat missing QUAL as 0
    }
    try {
        return std::stod(qualStr);
    } catch (const std::invalid_argument&) {
        std::cerr << "Invalid QUAL value: " << qualStr << "\n";
        return 0.0;
    }
}

int main(int argc, char* argv[]) {
    VCFXPhredFilter filter;
    return filter.run(argc, argv);
}
