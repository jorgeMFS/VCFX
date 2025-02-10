#include "VCFX_phase_checker.h"
#include <getopt.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>

int VCFXPhaseChecker::run(int argc, char* argv[]) {
    // If no arguments are provided, display help.
    if (argc == 1) {
        displayHelp();
        return 0;
    }
    
    bool showHelp = false;
    int opt;
    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    while ((opt = getopt(argc, argv, "h")) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            default:
                showHelp = true;
                break;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }
    processVCF(std::cin, std::cout);
    return 0;
}

void VCFXPhaseChecker::displayHelp() {
    std::cout <<
"VCFX_phase_checker: Output only VCF variant lines in which every sample genotype is fully phased.\n\n"
"Usage:\n"
"  ./VCFX_phase_checker [options] < input.vcf > phased_output.vcf\n\n"
"Options:\n"
"  -h, --help   Display this help message and exit\n\n"
"Description:\n"
"  The tool reads a VCF from standard input and checks the GT field (genotype) for each sample.\n"
"  A genotype is considered fully phased if it uses the '|' separator (e.g., \"0|1\") and contains\n"
"  no missing alleles. If every sample in a variant line is fully phased, the line is printed to\n"
"  standard output; otherwise, it is skipped and a warning is written to standard error.\n\n"
"Examples:\n"
"  To extract only fully phased variants:\n"
"    ./VCFX_phase_checker < input.vcf > phased_output.vcf\n";
}

bool VCFXPhaseChecker::isFullyPhased(const std::string &gt) const {
    if(gt.empty() || gt == "." || gt == "./." || gt == ".|.") return false;
    // A fully phased genotype should use '|' exclusively.
    if(gt.find('/') != std::string::npos) return false;
    std::vector<std::string> alleles;
    std::istringstream ss(gt);
    std::string token;
    while(std::getline(ss, token, '|')) {
        alleles.push_back(token);
    }
    if(alleles.empty()) return false;
    for(const auto &al : alleles) {
        if(al.empty() || al == ".") return false;
    }
    return true;
}

void VCFXPhaseChecker::processVCF(std::istream &in, std::ostream &out) {
    bool headerFound = false;
    std::string line;
    while (std::getline(in, line)) {
        if(line.empty()) {
            out << line << "\n";
            continue;
        }
        if(line[0] == '#') {
            out << line << "\n";
            if(line.rfind("#CHROM", 0) == 0) {
                headerFound = true;
            }
            continue;
        }
        if(!headerFound) {
            std::cerr << "Warning: Data line encountered before #CHROM header; skipping line.\n";
            continue;
        }
        std::vector<std::string> fields;
        std::istringstream ss(line);
        std::string f;
        while(std::getline(ss, f, '\t')) {
            fields.push_back(f);
        }
        if(fields.size() < 10) {
            std::cerr << "Warning: Invalid VCF line with fewer than 10 columns; skipping line.\n";
            continue;
        }
        std::vector<std::string> fmt;
        {
            std::istringstream fs(fields[8]);
            std::string token;
            while(std::getline(fs, token, ':')) {
                fmt.push_back(token);
            }
        }
        int gtIndex = -1;
        for (size_t i = 0; i < fmt.size(); ++i) {
            if (fmt[i] == "GT") {
                gtIndex = static_cast<int>(i);
                break;
            }
        }
        if (gtIndex == -1) {
            std::cerr << "Warning: GT field not found; skipping line.\n";
            continue;
        }
        bool allPhased = true;
        for (size_t i = 9; i < fields.size(); ++i) {
            std::vector<std::string> sampleFields;
            std::istringstream ssSample(fields[i]);
            std::string sub;
            while(std::getline(ssSample, sub, ':')) {
                sampleFields.push_back(sub);
            }
            if (gtIndex >= static_cast<int>(sampleFields.size())) {
                std::cerr << "Warning: GT index out of range in sample; skipping line.\n";
                allPhased = false;
                break;
            }
            if (!isFullyPhased(sampleFields[gtIndex])) {
                allPhased = false;
                break;
            }
        }
        if (allPhased) {
            out << line << "\n";
        } else {
            std::cerr << "Unphased genotype found at CHROM=" << fields[0] 
                      << ", POS=" << fields[1] << "; line skipped.\n";
        }
    }
}

int main(int argc, char* argv[]){
    VCFXPhaseChecker checker;
    return checker.run(argc, argv);
}
