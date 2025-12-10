#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

// SIMD support for fast newline scanning
#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#define VCFX_HAS_SSE2 1
#if defined(__AVX2__)
#define VCFX_HAS_AVX2 1
#endif
#elif defined(__aarch64__)
#include <arm_neon.h>
#define VCFX_HAS_NEON 1
#endif

// RAII wrapper for memory-mapped files
struct MappedFile {
    const char* data = nullptr;
    size_t size = 0;
    int fd = -1;

    MappedFile() = default;
    ~MappedFile() { close(); }

    MappedFile(const MappedFile&) = delete;
    MappedFile& operator=(const MappedFile&) = delete;

    bool open(const char* path) {
        fd = ::open(path, O_RDONLY);
        if (fd < 0) return false;

        struct stat st;
        if (fstat(fd, &st) < 0) {
            ::close(fd);
            fd = -1;
            return false;
        }
        size = static_cast<size_t>(st.st_size);
        if (size == 0) {
            ::close(fd);
            fd = -1;
            return true; // Empty file is valid
        }

        data = static_cast<const char*>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED) {
            data = nullptr;
            ::close(fd);
            fd = -1;
            return false;
        }

        // Advise kernel for sequential access
        madvise(const_cast<char*>(data), size, MADV_SEQUENTIAL | MADV_WILLNEED);
        return true;
    }

    void close() {
        if (data && size > 0) {
            munmap(const_cast<char*>(data), size);
        }
        if (fd >= 0) {
            ::close(fd);
        }
        data = nullptr;
        size = 0;
        fd = -1;
    }
};

// SIMD-optimized newline finder
static inline const char* findNewline(const char* start, const char* end) {
#if defined(VCFX_HAS_AVX2)
    const __m256i nl = _mm256_set1_epi8('\n');
    while (start + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(start));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, nl);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask != 0) {
            return start + __builtin_ctz(mask);
        }
        start += 32;
    }
#elif defined(VCFX_HAS_SSE2)
    const __m128i nl = _mm_set1_epi8('\n');
    while (start + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i*>(start));
        __m128i cmp = _mm_cmpeq_epi8(chunk, nl);
        int mask = _mm_movemask_epi8(cmp);
        if (mask != 0) {
            return start + __builtin_ctz(mask);
        }
        start += 16;
    }
#elif defined(VCFX_HAS_NEON)
    const uint8x16_t nl = vdupq_n_u8('\n');
    while (start + 16 <= end) {
        uint8x16_t chunk = vld1q_u8(reinterpret_cast<const uint8_t*>(start));
        uint8x16_t cmp = vceqq_u8(chunk, nl);
        uint64_t mask0 = vgetq_lane_u64(vreinterpretq_u64_u8(cmp), 0);
        uint64_t mask1 = vgetq_lane_u64(vreinterpretq_u64_u8(cmp), 1);
        if (mask0) {
            return start + (__builtin_ctzll(mask0) >> 3);
        }
        if (mask1) {
            return start + 8 + (__builtin_ctzll(mask1) >> 3);
        }
        start += 16;
    }
#endif
    // Scalar fallback
    const char* nl_pos = static_cast<const char*>(std::memchr(start, '\n', end - start));
    return nl_pos ? nl_pos : end;
}

// Extract nth tab-delimited field from a line (0-indexed)
static inline std::string_view getNthField(std::string_view line, size_t n) {
    size_t start = 0;
    size_t fieldIdx = 0;

    for (size_t i = 0; i <= line.size(); ++i) {
        if (i == line.size() || line[i] == '\t') {
            if (fieldIdx == n) {
                return line.substr(start, i - start);
            }
            fieldIdx++;
            start = i + 1;
        }
    }
    return std::string_view{};
}

// ----------------------------------------------------
// A helper struct for storing freq data by population
// ----------------------------------------------------
struct PopFreqKey {
    std::string key;
    double frequency;
};

// Map from variant key to map of pop->frequency for efficient lookup
typedef std::unordered_map<std::string, double> FrequencyMap;
typedef std::unordered_map<std::string, std::unordered_map<std::string, double>> VariantPopFreqMap;

// ----------------------------------------------------
// Class: VCFXAncestryInferrer
// ----------------------------------------------------
class VCFXAncestryInferrer {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();
    bool loadPopulationFrequencies(const std::string &freqFilePath);
    bool inferAncestry(std::istream &vcfInput, std::ostream &output);
    bool inferAncestryMmap(const char* filepath, std::ostream &output);

  private:
    FrequencyMap freqData;
    VariantPopFreqMap variantPopFreqs;
    std::set<std::string> populations;
    std::string inputFile;
    bool quiet = false;
};

// ----------------------------------------------------
// main()
// ----------------------------------------------------
static void show_help() {
    VCFXAncestryInferrer obj;
    char arg0[] = "VCFX_ancestry_inferrer";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_ancestry_inferrer", show_help))
        return 0;
    VCFXAncestryInferrer inferrer;
    return inferrer.run(argc, argv);
}

// ----------------------------------------------------
// run()
// ----------------------------------------------------
int VCFXAncestryInferrer::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;
    std::string freqFilePath;

    static struct option longOpts[] = {
        {"help", no_argument, 0, 'h'},
        {"frequency", required_argument, 0, 'f'},
        {"input", required_argument, 0, 'i'},
        {"quiet", no_argument, 0, 'q'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "hf:i:q", longOpts, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'f':
            freqFilePath = optarg;
            break;
        case 'i':
            inputFile = optarg;
            break;
        case 'q':
            quiet = true;
            break;
        default:
            showHelp = true;
        }
    }

    // Check for positional argument (input file)
    if (inputFile.empty() && optind < argc) {
        inputFile = argv[optind];
    }

    if (showHelp || freqFilePath.empty()) {
        displayHelp();
        return (showHelp ? 0 : 1);
    }

    // Load population frequencies
    if (!loadPopulationFrequencies(freqFilePath)) {
        std::cerr << "Error: Failed to load population frequencies from " << freqFilePath << "\n";
        return 1;
    }

    // Use mmap if input file specified
    if (!inputFile.empty()) {
        if (!inferAncestryMmap(inputFile.c_str(), std::cout)) {
            return 1;
        }
    } else {
        // Read VCF from stdin
        if (!inferAncestry(std::cin, std::cout)) {
            return 1;
        }
    }

    return 0;
}

// ----------------------------------------------------
// displayHelp
// ----------------------------------------------------
void VCFXAncestryInferrer::displayHelp() {
    std::cout << "VCFX_ancestry_inferrer: Infer population ancestry based on allele frequencies.\n\n"
              << "Usage:\n"
              << "  VCFX_ancestry_inferrer --frequency <freq_file> -i input.vcf > ancestry.txt\n"
              << "  VCFX_ancestry_inferrer --frequency <freq_file> < input.vcf > ancestry.txt\n\n"
              << "Options:\n"
              << "  -f, --frequency FILE   Population frequency file (required)\n"
              << "  -i, --input FILE       Input VCF file (uses mmap for better performance)\n"
              << "  -q, --quiet            Suppress warning messages\n"
              << "  -h, --help             Display this help message\n\n"
              << "Description:\n"
              << "  Reads a VCF and outputs a 2-column table:\n"
              << "    Sample  Inferred_Population\n\n"
              << "  The frequency file must have lines of the form:\n"
              << "    CHROM  POS  REF  ALT  POPULATION  FREQUENCY\n"
              << "  (tab-separated). For multi-allelic VCF sites, an ALT allele index 1\n"
              << "  corresponds to the first item in the comma-separated ALT list,\n"
              << "  index 2 => second ALT, etc.\n\n"
              << "Example:\n"
              << "  VCFX_ancestry_inferrer -f pop_frequencies.txt -i input.vcf > ancestry.txt\n";
}

// ----------------------------------------------------
// loadPopulationFrequencies
// ----------------------------------------------------
bool VCFXAncestryInferrer::loadPopulationFrequencies(const std::string &freqFilePath) {
    std::ifstream freqFile(freqFilePath);
    if (!freqFile.is_open()) {
        std::cerr << "Error: Cannot open frequency file: " << freqFilePath << "\n";
        return false;
    }

    int lineNum = 0;
    std::string line;
    while (std::getline(freqFile, line)) {
        lineNum++;
        if (line.empty()) continue;

        std::vector<std::string> fields;
        vcfx::split_tabs(line, fields);
        if (fields.size() < 6) {
            if (!quiet) {
                std::cerr << "Warning: Invalid line in frequency file (#" << lineNum << "): " << line << "\n";
            }
            continue;
        }

        const std::string &chrom = fields[0];
        const std::string &pos = fields[1];
        const std::string &ref = fields[2];
        const std::string &alt = fields[3];
        const std::string &pop = fields[4];
        const std::string &freqStr = fields[5];

        double freq = 0.0;
        try {
            freq = std::stod(freqStr);
        } catch (...) {
            if (!quiet) {
                std::cerr << "Warning: Invalid frequency value in line #" << lineNum << ": " << line << "\n";
            }
            continue;
        }

        std::string key = chrom + ":" + pos + ":" + ref + ":" + alt + ":" + pop;
        freqData[key] = freq;

        std::string variantKey = chrom + ":" + pos + ":" + ref + ":" + alt;
        variantPopFreqs[variantKey][pop] = freq;

        populations.insert(pop);
    }
    freqFile.close();

    if (freqData.empty()) {
        std::cerr << "Error: No valid population frequencies loaded.\n";
        return false;
    }
    return true;
}

// ----------------------------------------------------
// inferAncestryMmap - Memory-mapped version
// ----------------------------------------------------
bool VCFXAncestryInferrer::inferAncestryMmap(const char* filepath, std::ostream &out) {
    MappedFile mf;
    if (!mf.open(filepath)) {
        std::cerr << "Error: Cannot open file: " << filepath << "\n";
        return false;
    }

    if (mf.size == 0) {
        std::cerr << "Error: Empty VCF file\n";
        return false;
    }

    const char* pos = mf.data;
    const char* end = mf.data + mf.size;
    bool foundChromHeader = false;
    std::vector<std::string> sampleNames;

    // Sample -> (pop -> score)
    std::unordered_map<std::string, std::unordered_map<std::string, double>> sampleScores;

    while (pos < end) {
        const char* lineEnd = findNewline(pos, end);
        std::string_view line(pos, lineEnd - pos);

        if (line.empty()) {
            pos = lineEnd + 1;
            continue;
        }

        if (line[0] == '#') {
            // Check for #CHROM header
            if (line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
                foundChromHeader = true;

                // Parse sample columns (columns 9+)
                size_t colIdx = 0;
                size_t fieldStart = 0;
                for (size_t i = 0; i <= line.size(); ++i) {
                    if (i == line.size() || line[i] == '\t') {
                        if (colIdx >= 9) {
                            std::string sampleName(line.substr(fieldStart, i - fieldStart));
                            sampleNames.push_back(sampleName);
                            for (const auto &pop : populations) {
                                sampleScores[sampleName][pop] = 0.0;
                            }
                        }
                        colIdx++;
                        fieldStart = i + 1;
                    }
                }
            }
            pos = lineEnd + 1;
            continue;
        }

        if (!foundChromHeader) {
            std::cerr << "Error: Encountered VCF data before #CHROM header.\n";
            return false;
        }

        // Parse data line
        std::string_view chrom = getNthField(line, 0);
        std::string_view posField = getNthField(line, 1);
        std::string_view ref = getNthField(line, 3);
        std::string_view altStr = getNthField(line, 4);
        std::string_view format = getNthField(line, 8);

        if (chrom.empty() || posField.empty() || altStr.empty()) {
            pos = lineEnd + 1;
            continue;
        }

        // Split ALT by comma for multi-allelic
        std::vector<std::string_view> altAlleles;
        size_t altStart = 0;
        for (size_t i = 0; i <= altStr.size(); ++i) {
            if (i == altStr.size() || altStr[i] == ',') {
                altAlleles.push_back(altStr.substr(altStart, i - altStart));
                altStart = i + 1;
            }
        }

        if (altAlleles.empty()) {
            pos = lineEnd + 1;
            continue;
        }

        // Find GT index in format
        int gtIndex = -1;
        size_t fmtIdx = 0;
        size_t fmtStart = 0;
        for (size_t i = 0; i <= format.size(); ++i) {
            if (i == format.size() || format[i] == ':') {
                if (format.substr(fmtStart, i - fmtStart) == "GT") {
                    gtIndex = static_cast<int>(fmtIdx);
                    break;
                }
                fmtIdx++;
                fmtStart = i + 1;
            }
        }

        if (gtIndex < 0) {
            pos = lineEnd + 1;
            continue;
        }

        // Parse samples (columns 9+)
        size_t colIdx = 0;
        size_t fieldStart = 0;
        size_t sampleIdx = 0;

        for (size_t i = 0; i <= line.size(); ++i) {
            if (i == line.size() || line[i] == '\t') {
                if (colIdx >= 9 && sampleIdx < sampleNames.size()) {
                    std::string_view sampleData = line.substr(fieldStart, i - fieldStart);

                    // Extract GT field
                    int gtFieldIdx = 0;
                    size_t gtStart = 0;
                    std::string_view genotype;
                    for (size_t j = 0; j <= sampleData.size(); ++j) {
                        if (j == sampleData.size() || sampleData[j] == ':') {
                            if (gtFieldIdx == gtIndex) {
                                genotype = sampleData.substr(gtStart, j - gtStart);
                                break;
                            }
                            gtFieldIdx++;
                            gtStart = j + 1;
                        }
                    }

                    if (!genotype.empty()) {
                        // Parse allele numbers (e.g., "0/1" or "2|1")
                        size_t alleleStart = 0;
                        for (size_t j = 0; j <= genotype.size(); ++j) {
                            if (j == genotype.size() || genotype[j] == '/' || genotype[j] == '|') {
                                std::string_view alleleStr = genotype.substr(alleleStart, j - alleleStart);
                                alleleStart = j + 1;

                                if (alleleStr.empty() || alleleStr[0] == '.') continue;

                                // Parse allele number
                                int aVal = 0;
                                bool valid = true;
                                for (char c : alleleStr) {
                                    if (c >= '0' && c <= '9') {
                                        aVal = aVal * 10 + (c - '0');
                                    } else {
                                        valid = false;
                                        break;
                                    }
                                }

                                if (!valid || aVal == 0) continue; // REF allele or invalid

                                if (aVal > 0 && static_cast<size_t>(aVal) <= altAlleles.size()) {
                                    std::string_view actualAlt = altAlleles[aVal - 1];

                                    // Build variant key
                                    std::string variantKey;
                                    variantKey.reserve(chrom.size() + posField.size() + ref.size() + actualAlt.size() + 4);
                                    variantKey.append(chrom);
                                    variantKey.push_back(':');
                                    variantKey.append(posField);
                                    variantKey.push_back(':');
                                    variantKey.append(ref);
                                    variantKey.push_back(':');
                                    variantKey.append(actualAlt);

                                    // Lookup in frequency map
                                    auto variantIt = variantPopFreqs.find(variantKey);
                                    if (variantIt != variantPopFreqs.end()) {
                                        const auto &popFreqs = variantIt->second;

                                        // Find population with highest frequency
                                        double bestFreq = -1.0;
                                        std::string bestPop;

                                        for (const auto &popFreq : popFreqs) {
                                            if (popFreq.second > bestFreq) {
                                                bestFreq = popFreq.second;
                                                bestPop = popFreq.first;
                                            }
                                        }

                                        if (!bestPop.empty() && bestFreq >= 0.0) {
                                            sampleScores[sampleNames[sampleIdx]][bestPop] += bestFreq;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    sampleIdx++;
                }
                colIdx++;
                fieldStart = i + 1;
            }
        }

        pos = lineEnd + 1;
    }

    if (!foundChromHeader) {
        std::cerr << "Error: No valid VCF data found in input.\n";
        return false;
    }

    // Output results
    out << "Sample\tInferred_Population\n";
    for (const auto &sName : sampleNames) {
        auto it = sampleScores.find(sName);
        if (it == sampleScores.end() || it->second.empty()) {
            out << sName << "\tUnknown\n";
            continue;
        }

        const auto &popMap = it->second;
        std::string bestPop = "Unknown";
        double bestScore = -1.0;

        for (const auto &ps : popMap) {
            if (ps.second > bestScore) {
                bestScore = ps.second;
                bestPop = ps.first;
            }
        }

        out << sName << "\t" << bestPop << "\n";
    }

    return true;
}

// ----------------------------------------------------
// inferAncestry - stdin version (fallback)
// ----------------------------------------------------
bool VCFXAncestryInferrer::inferAncestry(std::istream &vcfInput, std::ostream &out) {
    std::string line;
    bool foundChromHeader = false;
    std::vector<std::string> headerFields;
    std::vector<std::string> sampleNames;

    std::unordered_map<std::string, std::unordered_map<std::string, double>> sampleScores;

    std::vector<std::string> fields;
    fields.reserve(16);

    if (!vcfInput) {
        std::cerr << "Error: Invalid VCF input stream.\n";
        return false;
    }

    int lineNum = 0;
    while (std::getline(vcfInput, line)) {
        lineNum++;
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
                headerFields.clear();
                vcfx::split_tabs(line, headerFields);

                if (headerFields.size() < 9) {
                    std::cerr << "Error: Invalid VCF header format. Expected at least 9 columns.\n";
                    return false;
                }

                for (size_t c = 9; c < headerFields.size(); ++c) {
                    sampleNames.push_back(headerFields[c]);
                    for (const auto &pop : populations) {
                        sampleScores[headerFields[c]][pop] = 0.0;
                    }
                }
            }
            continue;
        }

        if (!foundChromHeader) {
            std::cerr << "Error: Encountered VCF data before #CHROM header.\n";
            return false;
        }

        fields.clear();
        vcfx::split_tabs(line, fields);

        if (fields.size() < 10) {
            if (!quiet) {
                std::cerr << "Warning: Line " << lineNum << " has fewer than 10 columns, skipping.\n";
            }
            continue;
        }

        const std::string &chrom = fields[0];
        const std::string &pos = fields[1];
        const std::string &ref = fields[3];
        const std::string &altStr = fields[4];
        const std::string &format = fields[8];

        // Split ALT by comma
        std::vector<std::string> altAlleles;
        {
            std::stringstream altSS(altStr);
            std::string altTok;
            while (std::getline(altSS, altTok, ',')) {
                altAlleles.push_back(altTok);
            }
        }

        if (altAlleles.empty()) continue;

        // Find GT index in format
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
        if (gtIndex < 0) continue;

        // For each sample
        for (size_t s = 0; s < sampleNames.size(); ++s) {
            size_t sampleCol = 9 + s;
            if (sampleCol >= fields.size()) continue;

            const std::string &sampleData = fields[sampleCol];

            std::vector<std::string> sampleParts;
            {
                std::stringstream sampSS(sampleData);
                std::string part;
                while (std::getline(sampSS, part, ':')) {
                    sampleParts.push_back(part);
                }
            }
            if (gtIndex >= (int)sampleParts.size()) continue;

            std::string genotype = sampleParts[gtIndex];
            std::replace(genotype.begin(), genotype.end(), '|', '/');

            std::vector<std::string> alleleNums;
            {
                std::stringstream gtSS(genotype);
                std::string a;
                while (std::getline(gtSS, a, '/')) {
                    alleleNums.push_back(a);
                }
            }
            if (alleleNums.empty()) continue;

            for (auto &aStr : alleleNums) {
                if (aStr.empty() || aStr == ".") continue;

                bool numeric = true;
                for (char c : aStr) {
                    if (!isdigit(c)) {
                        numeric = false;
                        break;
                    }
                }
                if (!numeric) continue;

                int aVal = std::stoi(aStr);
                if (aVal == 0) continue; // REF allele

                if (aVal > 0 && (size_t)aVal <= altAlleles.size()) {
                    std::string actualAlt = altAlleles[aVal - 1];
                    std::string variantKey = chrom + ":" + pos + ":" + ref + ":" + actualAlt;

                    auto variantIt = variantPopFreqs.find(variantKey);
                    if (variantIt != variantPopFreqs.end()) {
                        const auto &popFreqs = variantIt->second;

                        double bestFreq = -1.0;
                        std::string bestPop;

                        for (const auto &popFreq : popFreqs) {
                            if (popFreq.second > bestFreq) {
                                bestFreq = popFreq.second;
                                bestPop = popFreq.first;
                            }
                        }

                        if (!bestPop.empty() && bestFreq >= 0.0) {
                            sampleScores[sampleNames[s]][bestPop] += bestFreq;
                        }
                    }
                }
            }
        }
    }

    if (!foundChromHeader) {
        std::cerr << "Error: No valid VCF data found in input.\n";
        return false;
    }

    // Output results
    out << "Sample\tInferred_Population\n";
    for (auto &sName : sampleNames) {
        auto it = sampleScores.find(sName);
        if (it == sampleScores.end() || it->second.empty()) {
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

    return true;
}
