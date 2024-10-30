#include "VCFX_region_subsampler.h"
#include <getopt.h>
#include <sstream>
#include <fstream>
#include <algorithm>

// Implementation of VCFXRegionSubsampler
int VCFXRegionSubsampler::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string bedFilePath;

    static struct option long_options[] = {
        {"help",         no_argument,       0, 'h'},
        {"region-bed",   required_argument, 0, 'b'},
        {0,              0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hb:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'b':
                bedFilePath = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || bedFilePath.empty()) {
        displayHelp();
        return 1;
    }

    // Load regions from BED file
    std::unordered_map<std::string, std::vector<std::pair<int, int>>> regions;
    if (!loadRegions(bedFilePath, regions)) {
        std::cerr << "Error: Failed to load regions from " << bedFilePath << "\n";
        return 1;
    }

    // Subsample regions from VCF
    subsampleRegions(std::cin, std::cout, regions);

    return 0;
}

void VCFXRegionSubsampler::displayHelp() {
    std::cout << "VCFX_region_subsampler: Subsample variants from specific genomic regions defined in a BED file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_region_subsampler --region-bed <regions.bed> [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                Display this help message and exit\n";
    std::cout << "  -b, --region-bed <bed>    Specify the BED file with genomic regions\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_region_subsampler --region-bed regions.bed < input.vcf > subsampled.vcf\n";
}

bool VCFXRegionSubsampler::loadRegions(const std::string& bedFilePath, std::unordered_map<std::string, std::vector<std::pair<int, int>>>& regions) {
    std::ifstream infile(bedFilePath);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open BED file " << bedFilePath << "\n";
        return false;
    }

    std::string line;
    size_t line_num = 0;
    while (std::getline(infile, line)) {
        line_num++;
        if (line.empty() || line[0] == '#') {
            continue; // Skip comments and empty lines
        }

        std::stringstream ss(line);
        std::string chrom;
        int start, end;

        if (!(ss >> chrom >> start >> end)) {
            std::cerr << "Warning: Skipping invalid BED line " << line_num << ": " << line << "\n";
            continue;
        }

        // BED format uses 0-based start and 1-based end
        regions[chrom].emplace_back(std::make_pair(start + 1, end));
    }

    infile.close();
    return true;
}

bool VCFXRegionSubsampler::isVariantInRegions(const std::string& chrom, int pos, const std::unordered_map<std::string, std::vector<std::pair<int, int>>>& regions) {
    auto it = regions.find(chrom);
    if (it == regions.end()) {
        return false;
    }

    const auto& regionList = it->second;
    for (const auto& region : regionList) {
        if (pos >= region.first && pos <= region.second) {
            return true;
        }
    }

    return false;
}

void VCFXRegionSubsampler::subsampleRegions(std::istream& in, std::ostream& out, const std::unordered_map<std::string, std::vector<std::pair<int, int>>>& regions) {
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            out << line << "\n";
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, posStr, id, ref, alt, qual, filter, info, format;
        if (!(ss >> chrom >> posStr >> id >> ref >> alt >> qual >> filter >> info >> format)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        int pos;
        try {
            pos = std::stoi(posStr);
        } catch (const std::invalid_argument&) {
            std::cerr << "Warning: Invalid position \"" << posStr << "\" in line: " << line << "\n";
            continue;
        }

        if (isVariantInRegions(chrom, pos, regions)) {
            // Output the variant line
            out << line << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXRegionSubsampler regionSubsampler;
    return regionSubsampler.run(argc, argv);
}