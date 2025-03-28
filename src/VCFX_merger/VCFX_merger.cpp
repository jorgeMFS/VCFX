#include "VCFX_merger.h"
#include <getopt.h>
#include <fstream>
#include <algorithm>
#include <map>
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
    std::vector<std::vector<std::string>> allVariants;
    std::vector<std::string> allHeaders;

    for (const auto& file : inputFiles) {
        std::vector<std::vector<std::string>> variants;
        std::vector<std::string> headerLines;
        parseVCF(file, variants, headerLines);

        // If no headers yet, copy the first file's headers
        if (allHeaders.empty()) {
            allHeaders = headerLines;
        }

        // Append all variants
        allVariants.insert(allVariants.end(), variants.begin(), variants.end());
    }

    // Sort all variants by chromosome and position
    std::sort(
        allVariants.begin(),
        allVariants.end(),
        [this](const std::vector<std::string>& a, const std::vector<std::string>& b) {
            if (a[0] == b[0]) {
                return std::stoi(a[1]) < std::stoi(b[1]);
            }
            return a[0] < b[0];
        }
    );

    // Output headers
    for (const auto& header : allHeaders) {
        out << header << "\n";
    }

    // Output merged variants
    for (const auto& variant : allVariants) {
        for (size_t i = 0; i < variant.size(); ++i) {
            out << variant[i];
            if (i < variant.size() - 1) {
                out << "\t";
            }
        }
        out << "\n";
    }
}

void VCFXMerger::parseVCF(const std::string& filename,
                          std::vector<std::vector<std::string>>& variants,
                          std::vector<std::string>& headerLines) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            headerLines.push_back(line);
            continue;
        }

        // Split by tab
        std::vector<std::string> fields;
        std::string field;
        size_t pos = 0;
        while ((pos = line.find('\t')) != std::string::npos) {
            field = line.substr(0, pos);
            fields.push_back(field);
            line.erase(0, pos + 1);
        }
        fields.push_back(line);

        variants.push_back(fields);
    }
    infile.close();
}

int main(int argc, char* argv[]) {
    VCFXMerger merger;
    return merger.run(argc, argv);
}
