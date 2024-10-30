#include "VCFX_nonref_filter.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>

// Implementation of VCFXNonRefFilter
int VCFXNonRefFilter::run(int argc, char* argv[]) {
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

    // Filter VCF input from stdin and output to stdout
    filterNonRef(std::cin, std::cout);

    return 0;
}

void VCFXNonRefFilter::displayHelp() {
    std::cout << "VCFX_nonref_filter: Filter out variants where all genotypes are homozygous reference.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_nonref_filter [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_nonref_filter < input.vcf > filtered.vcf\n";
}

void VCFXNonRefFilter::filterNonRef(std::istream& in, std::ostream& out) {
    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerLines;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            out << line << "\n"; // Preserve header lines
            if (line.substr(0, 6) == "#CHROM") {
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
            continue;
        }

        bool allHomRef = true;

        // Iterate over each sample to check genotypes
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
                allHomRef = false;
                break;
            }

            std::string genotype = sampleFields[gt_index];
            // Replace '|' with '/' for consistency
            std::replace(genotype.begin(), genotype.end(), '|', '/');

            // Check if genotype is homozygous reference
            if (genotype == "0/0" || genotype == "0") {
                continue;
            } else {
                allHomRef = false;
                break;
            }
        }

        if (!allHomRef) {
            out << line << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXNonRefFilter filter;
    return filter.run(argc, argv);
}
