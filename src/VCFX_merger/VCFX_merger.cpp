#include "VCFX_merger.h"
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <queue>
#include <algorithm>
#include <iostream>
#include <cstdlib>

// Implementation of VCFX_merger

int VCFXMerger::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::vector<std::string> inputFiles;

    static struct option long_options[] = {
        {"merge", required_argument, 0, 'm'},
        {"help",  no_argument,       0, 'h'},
        {0,       0,               0,   0 }
    };

    while ((opt = getopt_long(argc, argv, "m:h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'm': {
                // Split comma-separated file names
                std::string files = optarg;
                size_t pos = 0;
                while ((pos = files.find(',')) != std::string::npos) {
                    inputFiles.emplace_back(files.substr(0, pos));
                    files.erase(0, pos + 1);
                }
                inputFiles.emplace_back(files);
                break;
            }
            case 'h':
            default:
                showHelp = true;
        }
    }

    if (showHelp || inputFiles.empty()) {
        displayHelp();
        return 0;
    }

    // Merge VCF files and output to stdout
    mergeVCF(inputFiles, std::cout);
    return 0;
}

void VCFXMerger::displayHelp() {
    std::cout
        << "VCFX_merger: Merge multiple VCF files by variant position.\n\n"
        << "Usage:\n"
        << "  VCFX_merger --merge file1.vcf,file2.vcf,... [options]\n\n"
        << "Options:\n"
        << "  -m, --merge    Comma-separated list of VCF files to merge\n"
        << "  -h, --help     Display this help message and exit\n\n"
        << "Example:\n"
        << "  VCFX_merger --merge sample1.vcf,sample2.vcf > merged.vcf\n";
}

void VCFXMerger::mergeVCF(const std::vector<std::string>& inputFiles, std::ostream& out) {
    struct Variant {
        std::string chrom;
        long pos = 0;
        std::string line;
    };

    std::vector<Variant> variants;
    std::vector<std::string> headers;
    bool headersCaptured = false;

    for (const auto& file : inputFiles) {
        std::ifstream stream(file);
        if (!stream.is_open()) {
            std::cerr << "Failed to open file: " << file << "\n";
            continue;
        }

        std::string line;
        while (std::getline(stream, line)) {
            if (line.empty())
                continue;
            if (line[0] == '#') {
                if (!headersCaptured)
                    headers.push_back(line);
                continue;
            }

            std::istringstream ss(line);
            Variant v;
            std::getline(ss, v.chrom, '\t');
            std::string pos_str;
            std::getline(ss, pos_str, '\t');
            v.pos = std::strtol(pos_str.c_str(), nullptr, 10);
            v.line = line;
            variants.push_back(std::move(v));
        }

        if (!headersCaptured && !headers.empty())
            headersCaptured = true;
    }

    for (const auto& h : headers) {
        out << h << '\n';
    }

    std::sort(variants.begin(), variants.end(), [](const Variant& a, const Variant& b) {
        if (a.chrom == b.chrom) return a.pos < b.pos;
        return a.chrom < b.chrom;
    });

    for (const auto& v : variants) {
        out << v.line << '\n';
    }
}


int main(int argc, char* argv[]) {
    VCFXMerger merger;
    return merger.run(argc, argv);
}
