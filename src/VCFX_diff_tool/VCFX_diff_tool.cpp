#include "vcfx_core.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

// ----------------------------------------------------------------------
// A helper function to split a string by a delimiter
// ----------------------------------------------------------------------
static std::vector<std::string> split(const std::string &str, char delim) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string tmp;
    while (std::getline(ss, tmp, delim)) {
        tokens.push_back(tmp);
    }
    return tokens;
}

// ----------------------------------------------------------------------
// Class: VCFXDiffTool
// ----------------------------------------------------------------------
class VCFXDiffTool {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();

    bool loadVariants(const std::string &filePath, std::unordered_set<std::string> &variants);

    // We unify multi-allelic lines by splitting ALT on commas,
    // sorting them, and rejoining them.
    // The final key is "chrom:pos:ref:sortedAltString"
    std::string generateVariantKey(const std::string &chrom, const std::string &pos, const std::string &ref,
                                   const std::string &altField);
};

// ----------------------------------------------------------------------
// displayHelp
// ----------------------------------------------------------------------
void VCFXDiffTool::displayHelp() {
    std::cout << "VCFX_diff_tool: Compare two VCF files and identify differences.\n\n"
              << "Usage:\n"
              << "  VCFX_diff_tool --file1 <file1.vcf> --file2 <file2.vcf>\n\n"
              << "Options:\n"
              << "  -h, --help                Display this help message and exit\n"
              << "  -a, --file1 <file1.vcf>   Specify the first VCF file\n"
              << "  -b, --file2 <file2.vcf>   Specify the second VCF file\n\n"
              << "Example:\n"
              << "  VCFX_diff_tool --file1 file1.vcf --file2 file2.vcf\n";
}

// ----------------------------------------------------------------------
// generateVariantKey
//   - altField might be e.g. "G" or "G,T" or "T,GA"
//   - We split on commas, sort them, rejoin => e.g. "G,T"
//   - Then produce "chrom:pos:ref:thatSortedAlt"
// ----------------------------------------------------------------------
std::string VCFXDiffTool::generateVariantKey(const std::string &chrom, const std::string &pos, const std::string &ref,
                                             const std::string &altField) {
    auto alts = split(altField, ',');
    // sort them lexicographically
    std::sort(alts.begin(), alts.end());
    // rejoin
    std::ostringstream altSS;
    for (size_t i = 0; i < alts.size(); ++i) {
        if (i > 0)
            altSS << ",";
        altSS << alts[i];
    }
    std::string sortedAlt = altSS.str();

    return chrom + ":" + pos + ":" + ref + ":" + sortedAlt;
}

// ----------------------------------------------------------------------
// loadVariants
//   - Read each line, skip headers (#...)
//   - parse CHROM, POS, ID, REF, ALT
//   - generate key => add to set
// ----------------------------------------------------------------------
bool VCFXDiffTool::loadVariants(const std::string &filePath, std::unordered_set<std::string> &variants) {
    std::ifstream infile(filePath);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open file " << filePath << "\n";
        return false;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') {
            // skip headers and empty lines
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt;

        // minimal parse
        if (!(ss >> chrom >> pos >> id >> ref >> alt)) {
            std::cerr << "Warning: Skipping invalid VCF line:\n" << line << "\n";
            continue;
        }

        std::string key = generateVariantKey(chrom, pos, ref, alt);
        variants.insert(key);
    }
    return true;
}

// ----------------------------------------------------------------------
// run()
//   - parse args
//   - load sets from file1, file2
//   - compare => print differences
// ----------------------------------------------------------------------
int VCFXDiffTool::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;
    std::string file1Path;
    std::string file2Path;

    static struct option long_options[] = {{"help", no_argument, 0, 'h'},
                                           {"file1", required_argument, 0, 'a'},
                                           {"file2", required_argument, 0, 'b'},
                                           {0, 0, 0, 0}};

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
        return (showHelp ? 0 : 1);
    }

    // Load variants
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

    for (const auto &v : variantsFile1) {
        if (variantsFile2.find(v) == variantsFile2.end()) {
            uniqueToFile1.push_back(v);
        }
    }
    for (const auto &v : variantsFile2) {
        if (variantsFile1.find(v) == variantsFile1.end()) {
            uniqueToFile2.push_back(v);
        }
    }

    // Print results
    std::cout << "Variants unique to " << file1Path << ":\n";
    for (auto &v : uniqueToFile1) {
        std::cout << v << "\n";
    }
    std::cout << "\nVariants unique to " << file2Path << ":\n";
    for (auto &v : uniqueToFile2) {
        std::cout << v << "\n";
    }

    return 0;
}

// ----------------------------------------------------------------------
// main
// ----------------------------------------------------------------------
static void show_help() {
    VCFXDiffTool obj;
    char arg0[] = "VCFX_diff_tool";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    if (vcfx::handle_common_flags(argc, argv, "VCFX_diff_tool", show_help))
        return 0;
    VCFXDiffTool diffTool;
    return diffTool.run(argc, argv);
}
