#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdlib>
#include <getopt.h>

// ----------------------------------------------------
// A helper struct for storing freq data by population
// ----------------------------------------------------
struct PopFreqKey {
    // We store CHROM, POS, REF, ALT, POP as strings
    // so we can build a single key. Alternatively, we
    // can keep them separate or use a small struct.
    // We'll build a string key "chrom:pos:ref:alt:pop"
    // for direct indexing in an unordered_map.
    std::string key;
    double frequency;
};

// We’ll keep a direct structure: freqData["chr:pos:ref:alt:POP"] = frequency
typedef std::unordered_map<std::string, double> FrequencyMap;

// ----------------------------------------------------
// Class: VCFXAncestryInferrer
// ----------------------------------------------------
class VCFXAncestryInferrer {
public:
    int run(int argc, char* argv[]);

private:
    // Show usage
    void displayHelp();

    // Load population frequencies from a file that has lines:
    // CHROM  POS  REF  ALT  POPULATION  FREQUENCY
    bool loadPopulationFrequencies(const std::string& freqFilePath);

    // Infer ancestry from VCF
    void inferAncestry(std::istream& vcfInput, std::ostream& output);

private:
    // Frequencies keyed by "chr:pos:ref:alt:pop"
    FrequencyMap freqData;
};

// ----------------------------------------------------
// main() - create the inferrer and run
// ----------------------------------------------------
int main(int argc, char* argv[]) {
    VCFXAncestryInferrer inferrer;
    return inferrer.run(argc, argv);
}

// ----------------------------------------------------
// run() - parse arguments, load freq, run inference
// ----------------------------------------------------
int VCFXAncestryInferrer::run(int argc, char* argv[]) {
    int opt;
    bool showHelp = false;
    std::string freqFilePath;

    static struct option longOpts[] = {
        {"help",       no_argument,       0, 'h'},
        {"frequency",  required_argument, 0, 'f'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "hf:", longOpts, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'f':
                freqFilePath = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || freqFilePath.empty()) {
        displayHelp();
        // Return 0 if help was explicitly requested, otherwise 1
        return (showHelp ? 0 : 1);
    }

    // Load population frequencies
    if (!loadPopulationFrequencies(freqFilePath)) {
        std::cerr << "Error: Failed to load population frequencies from " << freqFilePath << "\n";
        return 1;
    }

    // Read VCF from stdin, write ancestry results to stdout
    inferAncestry(std::cin, std::cout);

    return 0;
}

// ----------------------------------------------------
// displayHelp
// ----------------------------------------------------
void VCFXAncestryInferrer::displayHelp() {
    std::cout 
        << "VCFX_ancestry_inferrer: Infer population ancestry based on allele frequencies.\n\n"
        << "Usage:\n"
        << "  VCFX_ancestry_inferrer --frequency <freq_file> [options]\n\n"
        << "Description:\n"
        << "  Reads a VCF from standard input and outputs a 2-column table:\n"
        << "    Sample  Inferred_Population\n\n"
        << "  The frequency file must have lines of the form:\n"
        << "    CHROM  POS  REF  ALT  POPULATION  FREQUENCY\n"
        << "  (tab-separated). For multi-allelic VCF sites, an ALT allele index 1\n"
        << "  corresponds to the first item in the comma-separated ALT list,\n"
        << "  index 2 => second ALT, etc.\n\n"
        << "Example:\n"
        << "  VCFX_ancestry_inferrer --frequency pop_frequencies.txt < input.vcf > ancestry_results.txt\n";
}

// ----------------------------------------------------
// loadPopulationFrequencies
//   freq file lines: CHROM, POS, REF, ALT, POP, FREQUENCY
// ----------------------------------------------------
bool VCFXAncestryInferrer::loadPopulationFrequencies(const std::string& freqFilePath) {
    std::ifstream freqFile(freqFilePath);
    if (!freqFile.is_open()) {
        std::cerr << "Error: Cannot open frequency file: " << freqFilePath << "\n";
        return false;
    }

    int lineNum = 0;
    std::string line;
    while (std::getline(freqFile, line)) {
        lineNum++;
        if (line.empty()) {
            continue;
        }

        // Parse line
        // e.g. "chr1   12345   A   G   EUR  0.45"
        std::stringstream ss(line);
        std::string chrom, pos, ref, alt, pop, freqStr;
        if (!(ss >> chrom >> pos >> ref >> alt >> pop >> freqStr)) {
            std::cerr << "Warning: Invalid line in frequency file (#" << lineNum << "): " << line << "\n";
            continue;
        }

        double freq = 0.0;
        try {
            freq = std::stod(freqStr);
        } catch (...) {
            std::cerr << "Warning: Invalid frequency value in line #" << lineNum << ": " << line << "\n";
            continue;
        }

        // Build a key: "chr:pos:ref:alt:pop"
        std::string key = chrom + ":" + pos + ":" + ref + ":" + alt + ":" + pop;
        freqData[key] = freq;
    }
    freqFile.close();

    if (freqData.empty()) {
        std::cerr << "Error: No valid population frequencies loaded.\n";
        return false;
    }
    return true;
}

// ----------------------------------------------------
// inferAncestry
//   1) parse the VCF header, get sample names
//   2) for each variant, parse REF + multi-ALT
//   3) for each sample's genotype, add to that sample's
//      population score = sum of freq for each ALT allele
//   4) after all lines, pick population with highest score
// ----------------------------------------------------
void VCFXAncestryInferrer::inferAncestry(std::istream& vcfInput, std::ostream& out) {
    std::string line;
    bool foundChromHeader = false;
    std::vector<std::string> headerFields;
    std::vector<std::string> sampleNames;

    // We'll accumulate ancestry scores: sample -> (pop -> score)
    // Then we pick the highest for each sample
    std::unordered_map<std::string, std::unordered_map<std::string, double>> sampleScores;

    while (std::getline(vcfInput, line)) {
        if (line.empty()) {
            continue;
        }

        // Check if line is a VCF header
        if (line[0] == '#') {
            // The #CHROM line is the last header line, we parse columns
            if (line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
                headerFields.clear();
                {
                    std::stringstream ss(line);
                    std::string tok;
                    while (std::getline(ss, tok, '\t')) {
                        headerFields.push_back(tok);
                    }
                }
                // sample columns start at index 9
                for (size_t c = 9; c < headerFields.size(); ++c) {
                    sampleNames.push_back(headerFields[c]);
                    // initialize each sample's population scores to 0
                    // We'll create the map on-the-fly later, so not strictly needed here
                }
            }
            continue;
        }

        // We must have #CHROM header before data
        if (!foundChromHeader) {
            std::cerr << "Error: Encountered VCF data before #CHROM header.\n";
            return;
        }

        // Parse data line
        // Format (minimally): CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [samples...]
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string f;
            while (std::getline(ss, f, '\t')) {
                fields.push_back(f);
            }
        }

        if (fields.size() < 10) {
            // not enough columns to have samples
            continue;
        }
        // Indices: 0=CHROM, 1=POS, 2=ID, 3=REF, 4=ALT, 5=QUAL, 6=FILTER, 7=INFO, 8=FORMAT, 9+ = samples
        const std::string &chrom = fields[0];
        const std::string &pos   = fields[1];
        // skip ID
        const std::string &ref   = fields[3];
        const std::string &altStr= fields[4];
        const std::string &format= fields[8];

        // Split ALT by comma for multi-allelic
        std::vector<std::string> altAlleles;
        {
            std::stringstream altSS(altStr);
            std::string altTok;
            while (std::getline(altSS, altTok, ',')) {
                altAlleles.push_back(altTok);
            }
        }

        // Find GT index in format
        // e.g. format might be "GT:DP:GQ:..."
        std::vector<std::string> formatParts;
        {
            std::stringstream fmts(format);
            std::string fmtTok;
            while (std::getline(fmts, fmtTok, ':')) {
                formatParts.push_back(fmtTok);
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
            // no genotype in this line
            continue;
        }

        // For each sample
        for (size_t s = 0; s < sampleNames.size(); ++s) {
            // sample data is fields[9 + s]
            size_t sampleCol = 9 + s;
            if (sampleCol >= fields.size()) {
                continue; 
            }
            const std::string &sampleData = fields[sampleCol];

            // Split by ':'
            std::vector<std::string> sampleParts;
            {
                std::stringstream sampSS(sampleData);
                std::string part;
                while (std::getline(sampSS, part, ':')) {
                    sampleParts.push_back(part);
                }
            }
            if (gtIndex >= (int)sampleParts.size()) {
                // no GT
                continue;
            }

            std::string genotype = sampleParts[gtIndex]; // e.g. "0/1" or "2|1"
            // unify separators
            std::replace(genotype.begin(), genotype.end(), '|', '/');
            std::vector<std::string> alleleNums;
            {
                std::stringstream gtSS(genotype);
                std::string a;
                while (std::getline(gtSS, a, '/')) {
                    alleleNums.push_back(a);
                }
            }
            // We expect diploid => 2 allele calls, but handle if 1 or 2
            if (alleleNums.empty()) {
                continue;
            }
            // For each allele: 0 => REF, 1 => altAlleles[0], 2 => altAlleles[1], etc.
            for (auto &aStr : alleleNums) {
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
                int aVal = std::stoi(aStr);
                if (aVal == 0) {
                    // ref => build a key with the REF as alt? 
                    // Typically we want alt freq, but let's say we do have lines for "chr:pos:ref:REF"? 
                    // Usually ancestry freq is about alt. If you want to incorporate ref, you'd store that
                    // in the freq file. We'll skip if we only track alt freqs. 
                    // For demonstration, let's skip 0 -> no alt allele => no ancestry added.
                    continue;
                }
                if (aVal > 0 && (size_t)aVal <= altAlleles.size()) {
                    // This alt
                    std::string actualAlt = altAlleles[aVal - 1];
                    // Now check each population freq. 
                    // We store freq in freqData keyed by "chr:pos:ref:alt:pop".
                    // We do not average them; we add the freq for each population? 
                    // Usually you'd pick the population with the highest freq. 
                    // Another approach is to add freq to that population's score. 
                    // Let's do "score += freq" for that population only. 
                    // That requires we look up each pop? 
                    // But we only have a single freq entry per pop in freqData. 
                    // Let's do a pass over freqData for all populations that have a key "chr:pos:ref:actualAlt:POP".
                    // Then we add that freq to sampleScores. This is a “sum all populations?” That doesn't make sense. 
                    // Usually you'd pick the single population with the highest freq. 
                    // Alternatively, you can add to them *all*, weighting by their freq. 
                    // But the code in question used "max freq" approach or "all freq"? 
                    // We'll do the "max freq" approach, consistent with #7 / #8 style.
                    
                    double bestFreq = -1.0;
                    std::string bestPop;
                    
                    // We must search each possible pop in freqData. But we store them individually as keys with pop at the end.
                    // Let's build a pattern "chrom:pos:ref:actualAlt:" and iterate over possible pops? 
                    // There's no direct iteration over freqData by prefix. So let's do a small trick:
                    
                    // We'll build the prefix:
                    std::string prefix = chrom + ":" + pos + ":" + ref + ":" + actualAlt + ":";
                    
                    // Now we can look up each population's freq by prefix + pop if we had a separate list of populations.
                    // But we do not store the list of populations here. So let's do an approach:
                    // We'll iterate over freqData and check if it starts with prefix, extracting the pop. This is not super efficient,
                    // but simple for demonstration. Or we can store a separate structure. 
                    // For a large dataset, we should store a multi-level structure or keep a set of populations. 
                    // We'll do the iteration approach for clarity.
                    
                    for (auto &kv : freqData) {
                        const std::string &fullKey = kv.first; // e.g. "chr1:12345:A:G:EUR"
                        if (fullKey.rfind(prefix, 0) == 0) {
                            // Extract the pop from the remainder
                            // prefix length is prefix.size(). The pop is what's after that
                            std::string pop = fullKey.substr(prefix.size());
                            double populationFreq = kv.second;
                            if (populationFreq > bestFreq) {
                                bestFreq = populationFreq;
                                bestPop = pop;
                            }
                        }
                    }
                    
                    if (!bestPop.empty() && bestFreq >= 0.0) {
                        // We add the bestFreq to the sample's pop
                        sampleScores[sampleNames[s]][bestPop] += bestFreq;
                    }
                } 
                // else if aVal > altAlleles.size() => skip 
            }
        }
    }

    // Done reading VCF. Now we pick the best population for each sample.
    out << "Sample\tInferred_Population\n";
    for (auto &sName : sampleNames) {
        // sampleScores[sName] => map<pop, score>
        auto it = sampleScores.find(sName);
        if (it == sampleScores.end()) {
            // no data => unknown
            out << sName << "\tUnknown\n";
            continue;
        }
        const auto &popMap = it->second;
        std::string bestPop = "Unknown";
        double bestScore = -1.0;
        for (auto &ps : popMap) {
            if (ps.second > bestScore) {
                bestScore = ps.second;
                bestPop = ps.first;
            }
        }
        out << sName << "\t" << bestPop << "\n";
    }
}
