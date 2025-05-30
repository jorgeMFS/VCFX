#include "VCFX_missing_detector.h"
#include "vcfx_core.h"
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @brief Main logic runner
 */
int VCFXMissingDetector::run(int argc, char *argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {{"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};

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
        return 0;
    }

    // Perform missing genotype detection
    detectMissingGenotypes(std::cin, std::cout);
    return 0;
}

/**
 * @brief Prints usage
 */
void VCFXMissingDetector::displayHelp() {
    std::cout << "VCFX_missing_detector: Detect variants with missing sample genotypes.\n\n"
              << "Usage:\n"
              << "  VCFX_missing_detector [options] < input.vcf > flagged.vcf\n\n"
              << "Options:\n"
              << "  -h, --help    Display this help message and exit\n\n"
              << "Description:\n"
              << "  Reads a VCF from stdin, checks each sample's genotype for missing data,\n"
              << "  and if any sample has a missing genotype, appends 'MISSING_GENOTYPES=1'\n"
              << "  in the INFO field.\n\n"
              << "Example:\n"
              << "  VCFX_missing_detector < input.vcf > flagged_output.vcf\n";
}

/**
 * @brief Helper function to check if a genotype string has missing data
 *
 * We consider a genotype missing if:
 *   - The entire string is '.' or './.' or '.|.'
 *   - Either allele is '.' => e.g. '0/.', './1'
 */
static bool genotypeIsMissing(const std::string &gt_field) {
    if (gt_field.empty())
        return true;

    // Extract just the genotype part (before the first colon)
    std::string gt = gt_field;
    size_t colon_pos = gt_field.find(':');
    if (colon_pos != std::string::npos) {
        gt = gt_field.substr(0, colon_pos);
    }

    if (gt == "." || gt == "./." || gt == ".|.")
        return true;

    // also check partial => find slash or pipe
    std::string tmp(gt);
    for (char &c : tmp) {
        if (c == '|')
            c = '/';
    }
    auto delim = tmp.find('/');
    if (delim == std::string::npos) {
        return false; // not diploid => we do not handle
    }
    std::string a1 = tmp.substr(0, delim);
    std::string a2 = tmp.substr(delim + 1);
    if (a1 == "." || a2 == ".") {
        return true;
    }
    return false;
}

/**
 * @brief The main function to detect missing genotypes.
 *        If any sample genotype is missing, we append 'MISSING_GENOTYPES=1' to the INFO field.
 */
void VCFXMissingDetector::detectMissingGenotypes(std::istream &in, std::ostream &out) {
    std::string line;
    while (true) {
        if (!std::getline(in, line))
            break;
        if (line.empty()) {
            out << "\n";
            continue;
        }

        if (line[0] == '#') {
            // pass header lines unchanged
            out << line << "\n";
            continue;
        }

        // parse columns
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string token;
            while (std::getline(ss, token, '\t')) {
                fields.push_back(token);
            }
        }
        if (fields.size() < 9) {
            // not a valid data line => just pass
            out << line << "\n";
            continue;
        }

        // Check FORMAT field (index 8) to find the position of GT
        int gtIndex = -1;
        std::string format = fields[8];
        if (!format.empty()) {
            std::istringstream format_ss(format);
            std::string token;
            int index = 0;
            while (std::getline(format_ss, token, ':')) {
                if (token == "GT") {
                    gtIndex = index;
                    break;
                }
                index++;
            }
        }

        // gather sample columns from index=9 onward
        bool hasMissing = false;
        for (size_t s = 9; s < fields.size(); s++) {
            std::string sample_value = fields[s];

            // If GT is the only FORMAT field, use the entire sample value
            // Otherwise, extract just the GT part based on its index
            std::string gt;
            if (format == "GT") {
                gt = sample_value;
            } else if (gtIndex >= 0) {
                std::istringstream sample_ss(sample_value);
                std::string token;
                int index = 0;
                while (std::getline(sample_ss, token, ':')) {
                    if (index == gtIndex) {
                        gt = token;
                        break;
                    }
                    index++;
                }
            } else {
                // No GT field found, can't determine if missing
                continue;
            }

            if (genotypeIsMissing(gt)) {
                hasMissing = true;
                break;
            }
        }

        if (!hasMissing) {
            // no missing => output as is
            out << line << "\n";
            continue;
        }

        // We do have missing => we append MISSING_GENOTYPES=1 to the INFO field
        std::string &info = fields[7];
        // If info is empty => info=="."
        // some lines might do "info==."" or info== "..."
        // if info=="." => we replace with "MISSING_GENOTYPES=1"
        if (info == "." || info.empty()) {
            info = "MISSING_GENOTYPES=1";
        } else {
            // ensure we have a semicolon if not already
            if (!info.empty() && info.back() != ';') {
                info.push_back(';');
            }
            info += "MISSING_GENOTYPES=1";
        }

        // rejoin
        // we produce the entire line
        // fields[0..8], then sample columns
        std::stringstream outLine;
        for (size_t i = 0; i < fields.size(); i++) {
            if (i > 0)
                outLine << "\t";
            outLine << fields[i];
        }
        out << outLine.str() << "\n";
    }
}

static void show_help() {
    VCFXMissingDetector obj;
    char arg0[] = "VCFX_missing_detector";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    if (vcfx::handle_common_flags(argc, argv, "VCFX_missing_detector", show_help))
        return 0;
    VCFXMissingDetector missingDetector;
    return missingDetector.run(argc, argv);
}
