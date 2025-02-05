#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdlib>
#include <getopt.h>

// ---------------------------------------------------------
// A helper struct to store command-line options if needed
// ---------------------------------------------------------
struct AncestryOptions {
    bool showHelp = false;
    std::string freqFile;
};

// ---------------------------------------------------------
// Class: VCFXAncestryAssigner
// ---------------------------------------------------------
class VCFXAncestryAssigner {
public:
    // High-level entry point
    int run(int argc, char* argv[]);

private:
    // Command-line usage
    void displayHelp();

    // Parse freq file
    bool loadAncestralFrequencies(std::istream& in);

    // Parse one frequency line
    // Format (tab-separated): CHROM  POS  REF  ALT  POP1_FREQ  POP2_FREQ ...
    bool parseFrequencyLine(const std::string& line);

    // Assign ancestry from VCF
    void assignAncestry(std::istream& vcfIn, std::ostream& out);

private:
    // Populations in order
    std::vector<std::string> populations;

    // Frequencies by key: "chrom:pos:ref:alt" => (pop => freq)
    std::unordered_map<std::string, std::unordered_map<std::string, double>> variantFrequencies;
};

// ---------------------------------------------------------
// VCFXAncestryAssigner::run
// ---------------------------------------------------------
int VCFXAncestryAssigner::run(int argc, char* argv[]) {
    // 1. Parse arguments
    int opt;
    bool showHelp = false;
    std::string freqFile;

    static struct option long_options[] = {
        {"help",            no_argument,       0, 'h'},
        {"assign-ancestry", required_argument, 0, 'a'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "ha:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'a':
                freqFile = std::string(optarg);
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || freqFile.empty()) {
        displayHelp();
        return (showHelp ? 0 : 1);
    }

    // 2. Open frequency file
    std::ifstream freqStream(freqFile);
    if (!freqStream.is_open()) {
        std::cerr << "Error: Unable to open frequency file: " << freqFile << "\n";
        return 1;
    }

    // 3. Load frequencies
    if (!loadAncestralFrequencies(freqStream)) {
        std::cerr << "Error: Failed to load ancestral frequencies.\n";
        return 1;
    }

    // 4. Assign ancestry based on VCF (read from stdin, write to stdout)
    assignAncestry(std::cin, std::cout);
    return 0;
}

// ---------------------------------------------------------
// Show usage
// ---------------------------------------------------------
void VCFXAncestryAssigner::displayHelp() {
    std::cout << "VCFX_ancestry_assigner: Assign samples to ancestral populations based on variant frequencies.\n\n"
              << "Usage:\n"
              << "  VCFX_ancestry_assigner --assign-ancestry <freq_file> < input.vcf > ancestry.txt\n\n"
              << "Options:\n"
              << "  -h, --help                 Show this help message and exit\n"
              << "  -a, --assign-ancestry FILE Ancestral frequency file\n\n"
              << "Frequency File Format:\n"
              << "  The first line must be a header like:\n"
              << "    CHROM  POS  REF  ALT  POP1  POP2  ...\n"
              << "  Each subsequent line must have the same columns. For example:\n"
              << "    1   10000   A   C   0.10  0.20\n\n"
              << "Example:\n"
              << "  VCFX_ancestry_assigner --assign-ancestry ancestral_freq.tsv < input.vcf > ancestry_out.txt\n\n";
}

// ---------------------------------------------------------
// parseFrequencyLine
// Expects: CHROM  POS  REF  ALT  pop1Freq  pop2Freq ...
// ---------------------------------------------------------
bool VCFXAncestryAssigner::parseFrequencyLine(const std::string& line) {
    std::stringstream ss(line);
    std::vector<std::string> fields;
    std::string field;
    while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
    }

    // Must have at least CHROM, POS, REF, ALT, plus the populations
    if (fields.size() < 4 + populations.size()) {
        return false;
    }
    const std::string &chrom = fields[0];
    int pos = 0;
    try {
        pos = std::stoi(fields[1]);
    } catch (...) {
        return false;
    }
    const std::string &ref = fields[2];
    const std::string &alt = fields[3];

    // Build freq map for this variant
    std::unordered_map<std::string, double> freqMap;
    for (size_t i = 0; i < populations.size(); ++i) {
        double freqVal = 0.0;
        try {
            freqVal = std::stod(fields[4 + i]);
        } catch (...) {
            freqVal = 0.0; // default if missing
        }
        freqMap[populations[i]] = freqVal;
    }

    // Key = CHROM:POS:REF:ALT
    const std::string key = chrom + ":" + std::to_string(pos) + ":" + ref + ":" + alt;
    variantFrequencies[key] = freqMap;
    return true;
}

// ---------------------------------------------------------
// loadAncestralFrequencies
// First line is header with columns:
//   CHROM  POS  REF  ALT  pop1  pop2 ...
// ---------------------------------------------------------
bool VCFXAncestryAssigner::loadAncestralFrequencies(std::istream& in) {
    std::string line;
    if (!std::getline(in, line)) {
        std::cerr << "Error: Frequency file is empty.\n";
        return false;
    }

    // Parse the header
    {
        std::stringstream ss(line);
        std::vector<std::string> headers;
        std::string h;
        while (std::getline(ss, h, '\t')) {
            headers.push_back(h);
        }
        // We need at least 5 columns: CHROM, POS, REF, ALT, plus 1 pop
        if (headers.size() < 5) {
            std::cerr << "Error: Frequency header must have at least 5 columns.\n";
            return false;
        }
        // The populations start at column index 4
        for (size_t i = 4; i < headers.size(); ++i) {
            populations.push_back(headers[i]);
        }
    }

    // Parse each subsequent line
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (!parseFrequencyLine(line)) {
            std::cerr << "Warning: Skipping invalid frequency line:\n" << line << "\n";
        }
    }
    return true;
}

// ---------------------------------------------------------
// assignAncestry
// Reads VCF from vcfIn, writes "Sample <tab> AssignedPopulation" to out
// ---------------------------------------------------------
void VCFXAncestryAssigner::assignAncestry(std::istream& vcfIn, std::ostream& out) {
    std::string line;
    bool haveHeader = false;
    int chrIndex = -1, posIndex = -1, refIndex = -1, altIndex = -1;
    std::vector<std::string> sampleNames;

    // For each sample: sampleScores[sample][population] => cumulative ancestry score
    std::unordered_map<std::string, std::unordered_map<std::string, double>> sampleScores;

    while (std::getline(vcfIn, line)) {
        if (line.empty()) {
            continue;
        }

        // VCF headers
        if (line[0] == '#') {
            // If #CHROM line, parse for sample columns
            if (line.rfind("#CHROM", 0) == 0) {
                haveHeader = true;
                std::stringstream ss(line);
                std::vector<std::string> headers;
                std::string f;
                while (std::getline(ss, f, '\t')) {
                    headers.push_back(f);
                }
                // Identify columns
                for (size_t i = 0; i < headers.size(); ++i) {
                    if (headers[i] == "CHROM") chrIndex = (int)i;
                    else if (headers[i] == "POS") posIndex = (int)i;
                    else if (headers[i] == "REF") refIndex = (int)i;
                    else if (headers[i] == "ALT") altIndex = (int)i;
                }
                if (chrIndex < 0 || posIndex < 0 || refIndex < 0 || altIndex < 0) {
                    std::cerr << "Error: #CHROM header missing required columns CHROM POS REF ALT.\n";
                    return;
                }
                // Sample columns start at index 9
                for (size_t i = 9; i < headers.size(); ++i) {
                    sampleNames.push_back(headers[i]);
                    // Initialize all sample->pop scores to 0
                    for (auto &pop : populations) {
                        sampleScores[headers[i]][pop] = 0.0;
                    }
                }
            }
            continue; // skip printing these lines
        }

        // Must have #CHROM line before data
        if (!haveHeader) {
            std::cerr << "Error: VCF data encountered before #CHROM line.\n";
            return;
        }

        // Parse a data line
        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string f;
        while (std::getline(ss, f, '\t')) {
            fields.push_back(f);
        }
        if (fields.size() < 9) {
            // Not enough columns
            continue;
        }

        const std::string &chrom = fields[chrIndex];
        int posVal = 0;
        try {
            posVal = std::stoi(fields[posIndex]);
        } catch (...) {
            // invalid POS
            continue;
        }
        const std::string &ref = fields[refIndex];
        const std::string &altField = fields[altIndex];

        // Split ALT by comma for multi-allelic
        std::vector<std::string> alts;
        {
            std::stringstream altSS(altField);
            std::string a;
            while (std::getline(altSS, a, ',')) {
                alts.push_back(a);
            }
        }

        // Identify the GT field index from the FORMAT column
        const std::string &formatCol = fields[8];
        std::vector<std::string> formatParts;
        {
            std::stringstream fmts(formatCol);
            while (std::getline(fmts, f, ':')) {
                formatParts.push_back(f);
            }
        }
        int gtIndex = -1;
        for (size_t i = 0; i < formatParts.size(); ++i) {
            if (formatParts[i] == "GT") {
                gtIndex = (int)i;
                break;
            }
        }
        if (gtIndex < 0) {
            // No genotype data
            continue;
        }

        // For each sample
        for (size_t s = 0; s < sampleNames.size(); ++s) {
            // The sample field is at index (9 + s)
            size_t sampleFieldIndex = 9 + s;
            if (sampleFieldIndex >= fields.size()) {
                continue; // no data
            }
            const std::string &sampleStr = fields[sampleFieldIndex];
            // e.g. "0/1:..." => split by ':'
            std::vector<std::string> sampleParts;
            {
                std::stringstream sss(sampleStr);
                std::string x;
                while (std::getline(sss, x, ':')) {
                    sampleParts.push_back(x);
                }
            }
            if ((int)sampleParts.size() <= gtIndex) {
                // no GT for this sample
                continue;
            }
            // genotype string, e.g. "0/1"
            const std::string &genotype = sampleParts[gtIndex];

            // Replace '|' with '/' for consistency
            std::string gt = genotype;
            for (char &c : gt) {
                if (c == '|') c = '/';
            }

            // Split genotype on '/'
            // e.g. "1/2"
            std::vector<std::string> gtAlleles;
            {
                std::stringstream gts(gt);
                std::string tok;
                while (std::getline(gts, tok, '/')) {
                    gtAlleles.push_back(tok);
                }
            }
            if (gtAlleles.size() < 2) {
                // not a diploid or not parseable
                continue;
            }

            // For each allele in genotype:
            //   0 => REF, 1 => alts[0], 2 => alts[1], etc.
            // If an allele is numeric but >0, we find that alt, look up freq => add to ancestry
            for (size_t alleleIndex = 0; alleleIndex < gtAlleles.size(); ++alleleIndex) {
                const std::string &aStr = gtAlleles[alleleIndex];
                if (aStr.empty() || aStr == ".") {
                    continue; // missing
                }
                bool numeric = true;
                for (char c : aStr) {
                    if (!isdigit(c)) {
                        numeric = false;
                        break;
                    }
                }
                if (!numeric) {
                    continue;
                }
                int alleleInt = std::stoi(aStr);
                if (alleleInt <= 0) {
                    // 0 => ref
                    continue;
                }
                if ((size_t)alleleInt > alts.size()) {
                    // genotype says "3" but we only have 2 alts => skip
                    continue;
                }
                // The alt is alts[alleleInt - 1]
                const std::string &thisAlt = alts[alleleInt - 1];

                // Build key for freq
                std::string key = chrom + ":" + std::to_string(posVal) + ":" + ref + ":" + thisAlt;
                auto it = variantFrequencies.find(key);
                if (it == variantFrequencies.end()) {
                    // We have no freq data for this alt
                    continue;
                }
                // freqMap => pop => freq
                const auto &freqMap = it->second;

                // pick the population with the highest freq
                double maxFreq = -1.0;
                std::string bestPop;
                for (auto &popFreq : freqMap) {
                    if (popFreq.second > maxFreq) {
                        maxFreq = popFreq.second;
                        bestPop = popFreq.first;
                    }
                }
                if (bestPop.empty()) {
                    continue; // no freq data
                }
                // add score for bestPop
                // Each allele we see => + maxFreq
                // (2 alt copies if genotype=1/1 => we do it twice in loop)
                sampleScores[sampleNames[s]][bestPop] += maxFreq;
            }
        }
    }

    // After reading the VCF, pick the best population for each sample
    out << "Sample\tAssigned_Ancestry\n";
    for (auto &sample : sampleScores) {
        const std::string &sampleName = sample.first;
        double bestScore = -1.0;
        std::string bestPop = "NA";
        for (auto &p : sample.second) {
            if (p.second > bestScore) {
                bestScore = p.second;
                bestPop = p.first;
            }
        }
        out << sampleName << "\t" << bestPop << "\n";
    }
}

// ---------------------------------------------------------
// main() - just instantiate and run
// ---------------------------------------------------------
int main(int argc, char* argv[]) {
    VCFXAncestryAssigner assigner;
    return assigner.run(argc, argv);
}
