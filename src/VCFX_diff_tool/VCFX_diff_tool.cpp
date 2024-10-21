#include "VCFX_diff_tool.h"
#include <getopt.h>
#include <sstream>

// Implementation of VCFXDiffTool
int VCFXDiffTool::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string file1Path;
    std::string file2Path;

    static struct option long_options[] = {
        {"help",    no_argument,       0, 'h'},
        {"file1",   required_argument, 0, 'a'},
        {"file2",   required_argument, 0, 'b'},
        {0,         0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "ha:b:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'a':
                file1Path = optarg;
                break;
            case 'b':
                file2Path = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || file1Path.empty() || file2Path.empty()) {
        displayHelp();
        return 1;
    }

    // Load variants from both files
    std::unordered_set<std::string> variantsFile1;
    std::unordered_set<std::string> variantsFile2;

    if (!loadVariants(file1Path, variantsFile1)) {
        std::cerr << "Error: Failed to load variants from " << file1Path << "\n";
        return 1;
    }

    if (!loadVariants(file2Path, variantsFile2)) {
        std::cerr << "Error: Failed to load variants from " << file2Path << "\n";
        return 1;
    }

    // Identify differences
    std::vector<std::string> uniqueToFile1;
    std::vector<std::string> uniqueToFile2;

    for (const auto& variant : variantsFile1) {
        if (variantsFile2.find(variant) == variantsFile2.end()) {
            uniqueToFile1.push_back(variant);
        }
    }

    for (const auto& variant : variantsFile2) {
        if (variantsFile1.find(variant) == variantsFile1.end()) {
            uniqueToFile2.push_back(variant);
        }
    }

    // Output differences
    std::cout << "Variants unique to " << file1Path << ":\n";
    for (const auto& variant : uniqueToFile1) {
        std::cout << variant << "\n";
    }

    std::cout << "\nVariants unique to " << file2Path << ":\n";
    for (const auto& variant : uniqueToFile2) {
        std::cout << variant << "\n";
    }

    return 0;
}

void VCFXDiffTool::displayHelp() {
    std::cout << "VCFX_diff_tool: Compare two VCF files and identify differences.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_diff_tool --file1 <file1.vcf> --file2 <file2.vcf>\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                Display this help message and exit\n";
    std::cout << "  -a, --file1 <file1.vcf>   Specify the first VCF file\n";
    std::cout << "  -b, --file2 <file2.vcf>   Specify the second VCF file\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_diff_tool --file1 file1.vcf --file2 file2.vcf\n";
}

bool VCFXDiffTool::loadVariants(const std::string& filePath, std::unordered_set<std::string>& variants) {
    std::ifstream infile(filePath);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open file " << filePath << "\n";
        return false;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip headers and empty lines
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt;

        if (!(ss >> chrom >> pos >> id >> ref >> alt)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        std::string key = generateVariantKey(chrom, pos, ref, alt);
        variants.insert(key);
    }

    infile.close();
    return true;
}

std::string VCFXDiffTool::generateVariantKey(const std::string& chrom, const std::string& pos, const std::string& ref, const std::string& alt) {
    return chrom + ":" + pos + ":" + ref + ":" + alt;
}
