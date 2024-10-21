#include "VCFX_indel_normalizer.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>

// Implementation of VCFXIndelNormalizer
int VCFXIndelNormalizer::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {
        {"help",            no_argument,       0, 'h'},
        {0,                 0,                 0,  0 }
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

    // Perform indel normalization on stdin and output to stdout
    normalizeIndels(std::cin, std::cout);

    return 0;
}

void VCFXIndelNormalizer::displayHelp() {
    std::cout << "VCFX_indel_normalizer: Normalize indels to their left-most representation.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_indel_normalizer [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_indel_normalizer < input.vcf > normalized.vcf\n";
}

bool VCFXIndelNormalizer::isIndel(const std::string& ref, const std::string& alt) {
    return (ref.length() != alt.length());
}

bool VCFXIndelNormalizer::normalizeVariant(std::string& chrom, std::string& pos, std::string& ref, std::string& alt) {
    // Left-align the indel as per VCF specifications
    size_t start = 0;
    while (start < ref.length() && start < alt.length() && ref[ref.length() - 1 - start] == alt[alt.length() - 1 - start]) {
        start++;
    }

    if (start == 0) {
        return false; // No normalization needed
    }

    ref = ref.substr(0, ref.length() - start);
    alt = alt.substr(0, alt.length() - start);
    
    // Ensure there is at least one base before the indel
    if (ref.empty() || alt.empty()) {
        return false;
    }

    return true;
}

void VCFXIndelNormalizer::normalizeIndels(std::istream& in, std::ostream& out) {
    std::string line;
    bool headerPassed = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            out << line << "\n";
            if (line.substr(0, 6) == "#CHROM") {
                headerPassed = true;
            }
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info, format;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        if (isIndel(ref, alt)) {
            if (normalizeVariant(chrom, pos, ref, alt)) {
                // Reconstruct the VCF line with normalized indel
                std::string rest_of_line;
                getline(ss, rest_of_line);
                out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t"
                    << qual << "\t" << filter << "\t" << info << "\t" << format << rest_of_line << "\n";
                continue;
            }
        }

        // Output the original line if no normalization is needed
        out << line << "\n";
    }
}
