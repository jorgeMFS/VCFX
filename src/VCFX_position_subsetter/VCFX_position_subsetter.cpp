#include "VCFX_position_subsetter.h"
#include <sstream>
#include <vector>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_position_subsetter\n"
              << "Usage: VCFX_position_subsetter [OPTIONS]\n\n"
              << "Options:\n"
              << "  --region, -r \"CHR:START-END\"   Specify the genomic region to subset (e.g., \"chr1:10000-20000\").\n"
              << "  --help, -h                      Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Subsets VCF records based on the specified genomic region.\n\n"
              << "Examples:\n"
              << "  ./VCFX_position_subsetter --region \"chr1:10000-20000\" < input.vcf > subset.vcf\n";
}

// Helper structure to represent a genomic region
struct GenomicRegion {
    std::string chrom;
    int start;
    int end;
};

// Function to parse the region string
bool parseRegion(const std::string& region_str, GenomicRegion& region) {
    size_t colon = region_str.find(':');
    size_t dash = region_str.find('-');

    if (colon == std::string::npos || dash == std::string::npos || dash < colon) {
        std::cerr << "Error: Invalid region format. Expected format \"chrX:start-end\".\n";
        return false;
    }

    region.chrom = region_str.substr(0, colon);
    try {
        region.start = std::stoi(region_str.substr(colon + 1, dash - colon - 1));
        region.end = std::stoi(region_str.substr(dash + 1));
    } catch (...) {
        std::cerr << "Error: Unable to parse start or end positions.\n";
        return false;
    }

    if (region.start > region.end) {
        std::cerr << "Error: Start position is greater than end position.\n";
        return false;
    }

    return true;
}

// Function to subset VCF records based on genomic range
bool subsetVCFByPosition(std::istream& in, std::ostream& out, const std::string& region_str) {
    GenomicRegion region;
    if (!parseRegion(region_str, region)) {
        return false;
    }

    std::string line;
    bool header_found = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            out << line << "\n"; // Preserve header lines
            if (line.find("#CHROM") == 0) {
                header_found = true;
            }
            continue;
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        std::stringstream ss(line);
        std::string chrom, pos_str;
        // Extract CHROM and POS fields
        if (!(ss >> chrom >> pos_str)) {
            std::cerr << "Warning: Skipping invalid VCF line.\n";
            continue;
        }

        int pos = 0;
        try {
            pos = std::stoi(pos_str);
        } catch (...) {
            std::cerr << "Warning: Invalid POS value. Skipping line.\n";
            continue;
        }

        if (chrom == region.chrom && pos >= region.start && pos <= region.end) {
            out << line << "\n";
        }
    }

    return true;
}

int main(int argc, char* argv[]) {
    std::string region_str;

    // Argument parsing for help and region
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--region" || arg == "-r") && i + 1 < argc) {
            region_str = argv[++i];
        } else if (arg.find("--region=") == 0) {
            region_str = arg.substr(9);
        } else if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    if (region_str.empty()) {
        std::cerr << "Error: Genomic region not specified.\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    bool success = subsetVCFByPosition(std::cin, std::cout, region_str);
    return success ? 0 : 1;
}
