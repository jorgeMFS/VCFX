#include "VCFX_ref_comparator.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string_view>
#include <vector>

// Memory-mapped file support
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

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

// Buffered output for efficiency
class OutputBuffer {
public:
    explicit OutputBuffer(std::ostream& os, size_t bufSize = 1024 * 1024)
        : out(os), buffer(bufSize) {}

    ~OutputBuffer() { flush(); }

    void write(std::string_view sv) {
        if (pos + sv.size() > buffer.size()) {
            flush();
        }
        if (sv.size() > buffer.size()) {
            out.write(sv.data(), sv.size());
        } else {
            std::memcpy(buffer.data() + pos, sv.data(), sv.size());
            pos += sv.size();
        }
    }

    void write(char c) {
        if (pos >= buffer.size()) flush();
        buffer[pos++] = c;
    }

    void flush() {
        if (pos > 0) {
            out.write(buffer.data(), pos);
            pos = 0;
        }
    }

private:
    std::ostream& out;
    std::vector<char> buffer;
    size_t pos = 0;
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

// Fast integer parsing from string_view
static inline int parseIntFast(std::string_view sv) {
    int result = 0;
    for (char c : sv) {
        if (c >= '0' && c <= '9') {
            result = result * 10 + (c - '0');
        }
    }
    return result;
}

// A small helper to remove whitespace
static inline void stripSpaces(std::string &s) {
    s.erase(std::remove_if(s.begin(), s.end(), [](unsigned char c) { return std::isspace(c); }), s.end());
}

// Case-insensitive string comparison
static inline bool caseInsensitiveEqual(std::string_view a, std::string_view b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::toupper(static_cast<unsigned char>(a[i])) != std::toupper(static_cast<unsigned char>(b[i]))) {
            return false;
        }
    }
    return true;
}

int VCFXRefComparator::run(int argc, char *argv[]) {
    if (argc == 1) {
        displayHelp();
        return 0;
    }
    bool showHelp = false;
    bool quiet = false;
    std::string referencePath;
    std::string inputFile;

    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"reference", required_argument, 0, 'r'},
        {"input", required_argument, 0, 'i'},
        {"quiet", no_argument, 0, 'q'},
        {0, 0, 0, 0}
    };

    while (true) {
        int c = ::getopt_long(argc, argv, "hr:i:q", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 'r':
            referencePath = optarg;
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
    if (showHelp) {
        displayHelp();
        return 0;
    }
    if (referencePath.empty()) {
        std::cerr << "Error: must specify --reference <FASTA>.\n";
        displayHelp();
        return 1;
    }
    if (!loadReference(referencePath)) {
        std::cerr << "Error: failed to load reference from " << referencePath << "\n";
        return 1;
    }

    // Use mmap for file input (fast path)
    if (!inputFile.empty()) {
        if (!compareVCFMmap(inputFile.c_str(), std::cout, quiet)) {
            std::cerr << "Error: failed to process input file: " << inputFile << "\n";
            return 1;
        }
    } else {
        compareVCF(std::cin, std::cout);
    }
    return 0;
}

void VCFXRefComparator::displayHelp() {
    std::cout << "VCFX_ref_comparator: Compare VCF REF/ALT with a reference genome.\n\n"
                 "Usage:\n"
                 "  VCFX_ref_comparator --reference ref.fasta -i input.vcf > output.vcf\n"
                 "  VCFX_ref_comparator --reference ref.fasta < input.vcf > output.vcf\n\n"
                 "Options:\n"
                 "  -h, --help             Show this help.\n"
                 "  -r, --reference FILE   Reference FASTA file.\n"
                 "  -i, --input FILE       Input VCF file (uses mmap for better performance).\n"
                 "  -q, --quiet            Suppress warnings.\n\n"
                 "Description:\n"
                 "  Reads a reference FASTA into memory. Then reads each variant line:\n"
                 "   - If chromosome or position is invalid, logs a warning and sets REF_COMPARISON=UNKNOWN_CHROM or "
                 "INVALID_POS.\n"
                 "   - Otherwise, compares the VCF's REF vs the reference substring. Then for each ALT, indicates "
                 "'REF_MATCH' if ALT= reference substring or 'NOVEL'.\n"
                 "  The result is appended to the 'INFO' field as REF_COMPARISON=...\n\n"
                 "Example:\n"
                 "  VCFX_ref_comparator --reference genome.fa -i in.vcf > out.vcf\n";
}

bool VCFXRefComparator::loadReference(const std::string &referenceFastaPath) {
    std::ifstream in(referenceFastaPath);
    if (!in.is_open()) {
        std::cerr << "Error: cannot open reference " << referenceFastaPath << "\n";
        return false;
    }
    referenceGenome.clear();
    std::string line, currentChrom;
    std::ostringstream seq;
    while (true) {
        if (!std::getline(in, line))
            break;
        if (line.empty())
            continue;
        if (line[0] == '>') {
            // store old chrom
            if (!currentChrom.empty()) {
                referenceGenome[currentChrom] = seq.str();
            }
            seq.str("");
            seq.clear();
            // parse new chrom
            currentChrom = line.substr(1);
            // if there's whitespace, strip it after first token
            {
                std::istringstream iss(currentChrom);
                iss >> currentChrom;
            }
            // uppercase
            std::transform(currentChrom.begin(), currentChrom.end(), currentChrom.begin(), ::toupper);
        } else {
            stripSpaces(line);
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            seq << line;
        }
    }
    if (!currentChrom.empty()) {
        referenceGenome[currentChrom] = seq.str();
    }
    return true;
}

bool VCFXRefComparator::compareVCFMmap(const char* filepath, std::ostream& out, bool quiet) {
    MappedFile mf;
    if (!mf.open(filepath)) {
        return false;
    }

    if (mf.size == 0) {
        return true; // Empty file is valid
    }

    OutputBuffer outBuf(out);
    const char* pos = mf.data;
    const char* end = mf.data + mf.size;
    bool foundChromHeader = false;
    bool infoHeaderInserted = false;

    // Temp buffers for string operations
    std::string chromUpper;
    std::string refUpper;
    std::string altUpper;
    std::string outLine;
    outLine.reserve(16384);

    while (pos < end) {
        const char* lineEnd = findNewline(pos, end);
        std::string_view line(pos, lineEnd - pos);

        // Skip empty lines
        if (line.empty()) {
            outBuf.write('\n');
            pos = lineEnd + 1;
            continue;
        }

        // Handle header lines
        if (line[0] == '#') {
            if (line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
                foundChromHeader = true;
                // Insert new INFO line BEFORE we write #CHROM
                if (!infoHeaderInserted) {
                    outBuf.write("##INFO=<ID=REF_COMPARISON,Number=1,Type=String,Description=\"Comparison of REF/ALT vs reference genome substring\">\n");
                    infoHeaderInserted = true;
                }
            }
            outBuf.write(line);
            outBuf.write('\n');
            pos = lineEnd + 1;
            continue;
        }

        // Data line - check if header was found
        if (!foundChromHeader) {
            if (!quiet) {
                std::cerr << "Warning: data line encountered before #CHROM => skipping.\n";
            }
            pos = lineEnd + 1;
            continue;
        }

        // Parse fields
        std::string_view chromField = getNthField(line, 0);
        std::string_view posField = getNthField(line, 1);
        std::string_view idField = getNthField(line, 2);
        std::string_view refField = getNthField(line, 3);
        std::string_view altField = getNthField(line, 4);
        std::string_view qualField = getNthField(line, 5);
        std::string_view filterField = getNthField(line, 6);
        std::string_view infoField = getNthField(line, 7);

        if (chromField.empty() || posField.empty() || refField.empty()) {
            if (!quiet) {
                std::cerr << "Warning: VCF line has insufficient columns => skipping.\n";
            }
            pos = lineEnd + 1;
            continue;
        }

        // Uppercase chrom
        chromUpper.assign(chromField);
        std::transform(chromUpper.begin(), chromUpper.end(), chromUpper.begin(), ::toupper);

        // Parse position
        int varPos = parseIntFast(posField);

        // Build output line with comparison result
        outLine.clear();
        outLine.append(chromField);
        outLine.push_back('\t');
        outLine.append(posField);
        outLine.push_back('\t');
        outLine.append(idField);
        outLine.push_back('\t');
        outLine.append(refField);
        outLine.push_back('\t');
        outLine.append(altField);
        outLine.push_back('\t');
        outLine.append(qualField);
        outLine.push_back('\t');
        outLine.append(filterField);
        outLine.push_back('\t');

        // Build INFO with comparison
        std::string newInfo(infoField);

        // Find reference
        auto it = referenceGenome.find(chromUpper);
        if (it == referenceGenome.end()) {
            if (!newInfo.empty() && newInfo.back() != ';') newInfo += ';';
            newInfo += "REF_COMPARISON=UNKNOWN_CHROM";
        } else {
            const std::string& seq = it->second;
            if (varPos < 1 || varPos > static_cast<int>(seq.size())) {
                if (!newInfo.empty() && newInfo.back() != ';') newInfo += ';';
                newInfo += "REF_COMPARISON=INVALID_POS";
            } else {
                // Get genome reference
                refUpper.assign(refField);
                std::transform(refUpper.begin(), refUpper.end(), refUpper.begin(), ::toupper);

                std::string genomeRef = seq.substr(varPos - 1, refUpper.size());

                // Split alt by comma and compare each
                std::string comparisons;
                size_t altStart = 0;
                for (size_t i = 0; i <= altField.size(); ++i) {
                    if (i == altField.size() || altField[i] == ',') {
                        std::string_view alt = altField.substr(altStart, i - altStart);
                        altUpper.assign(alt);
                        std::transform(altUpper.begin(), altUpper.end(), altUpper.begin(), ::toupper);

                        if (!comparisons.empty()) comparisons += ',';
                        if (altUpper == genomeRef) {
                            comparisons += "REF_MATCH";
                        } else {
                            comparisons += "NOVEL";
                        }
                        altStart = i + 1;
                    }
                }

                if (!newInfo.empty() && newInfo.back() != ';') newInfo += ';';
                newInfo += "REF_COMPARISON=" + comparisons;
            }
        }

        outLine.append(newInfo);

        // Append remaining columns (samples)
        size_t fieldCount = 0;
        size_t fieldStart = 0;
        for (size_t i = 0; i <= line.size(); ++i) {
            if (i == line.size() || line[i] == '\t') {
                if (fieldCount >= 8) {
                    outLine.push_back('\t');
                    outLine.append(line.substr(fieldStart, i - fieldStart));
                }
                fieldCount++;
                fieldStart = i + 1;
            }
        }

        outBuf.write(outLine);
        outBuf.write('\n');
        pos = lineEnd + 1;
    }

    outBuf.flush();
    return true;
}

void VCFXRefComparator::compareVCF(std::istream &vcfIn, std::ostream &vcfOut) {
    bool foundChromHeader = false;
    infoHeaderInserted = false;
    std::string line;
    std::vector<std::string> fields;
    fields.reserve(16);
    while (true) {
        if (!std::getline(vcfIn, line))
            break;
        if (line.empty()) {
            vcfOut << line << "\n";
            continue;
        }
        if (line[0] == '#') {
            // check if #CHROM
            if (line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
                // Insert new INFO line BEFORE we write #CHROM
                if (!infoHeaderInserted) {
                    vcfOut << "##INFO=<ID=REF_COMPARISON,Number=1,Type=String,Description=\"Comparison of REF/ALT vs "
                              "reference genome substring\">\n";
                    infoHeaderInserted = true;
                }
                vcfOut << line << "\n";
            } else {
                vcfOut << line << "\n";
            }
            continue;
        }
        if (!foundChromHeader) {
            std::cerr << "Warning: data line encountered before #CHROM => skipping.\n";
            continue;
        }
        // parse fields
        vcfx::split_tabs(line, fields);
        if (fields.size() < 8) {
            std::cerr << "Warning: VCF line has <8 columns => skipping.\n";
            continue;
        }
        // fields: 0=CHROM,1=POS,2=ID,3=REF,4=ALT,5=QUAL,6=FILTER,7=INFO,...
        std::string &chrom = fields[0];
        std::string &posStr = fields[1];
        std::string &ref = fields[3];
        std::string &alt = fields[4];
        std::string &info = fields[7];

        // uppercase chrom
        std::transform(chrom.begin(), chrom.end(), chrom.begin(), ::toupper);

        int pos = 0;
        try {
            pos = std::stoi(posStr);
        } catch (...) {
            // out of range
            if (!info.empty() && info.back() != ';')
                info += ';';
            info += "REF_COMPARISON=INVALID_POS";
            // rewrite
            std::ostringstream newLine;
            for (int i = 0; i < 7; i++) {
                newLine << fields[i];
                if (i < 7)
                    newLine << "\t";
            }
            newLine << info;
            for (size_t i = 8; i < fields.size(); i++) {
                newLine << "\t" << fields[i];
            }
            vcfOut << newLine.str() << "\n";
            continue;
        }
        // find reference
        auto it = referenceGenome.find(chrom);
        if (it == referenceGenome.end()) {
            // unknown chrom
            if (!info.empty() && info.back() != ';')
                info += ';';
            info += "REF_COMPARISON=UNKNOWN_CHROM";
            // rewrite
            std::ostringstream newLine;
            for (int i = 0; i < 7; i++) {
                newLine << fields[i];
                if (i < 7)
                    newLine << "\t";
            }
            newLine << info;
            for (size_t i = 8; i < fields.size(); i++) {
                newLine << "\t" << fields[i];
            }
            vcfOut << newLine.str() << "\n";
            continue;
        }
        const std::string &seq = it->second;
        if (pos < 1 || pos > (int)seq.size()) {
            // invalid pos
            if (!info.empty() && info.back() != ';')
                info += ';';
            info += "REF_COMPARISON=INVALID_POS";
            // rewrite
            std::ostringstream newLine;
            for (int i = 0; i < 7; i++) {
                newLine << fields[i];
                if (i < 7)
                    newLine << "\t";
            }
            newLine << info;
            for (size_t i = 8; i < fields.size(); i++) {
                newLine << "\t" << fields[i];
            }
            vcfOut << newLine.str() << "\n";
            continue;
        }
        // ref from genome
        std::string genomeRef = seq.substr(pos - 1, ref.size()); // 1-based
        // uppercase
        std::transform(ref.begin(), ref.end(), ref.begin(), ::toupper);

        // compare ref vs genomeRef
        bool refMatch = (ref == genomeRef);

        // split alt by comma
        std::vector<std::string> altAlleles;
        {
            std::stringstream as(alt);
            std::string a;
            while (std::getline(as, a, ',')) {
                std::transform(a.begin(), a.end(), a.begin(), ::toupper);
                altAlleles.push_back(a);
            }
        }
        // for each alt, see if alt== genomeRef => "REF_MATCH" else "NOVEL"
        std::vector<std::string> comparisons;
        for (auto &a : altAlleles) {
            if (a == genomeRef)
                comparisons.push_back("REF_MATCH");
            else
                comparisons.push_back("NOVEL");
        }
        // build comparisonStr
        std::string comparisonStr;
        for (size_t i = 0; i < comparisons.size(); i++) {
            if (i > 0)
                comparisonStr += ",";
            comparisonStr += comparisons[i];
        }
        if (!info.empty() && info.back() != ';')
            info += ';';
        // e.g. "REF_COMPARISON=REF_MATCH,REF_MATCH" or "NOVEL" etc.
        if (!refMatch) {
            // optionally label mismatch => we won't specifically do that
            // the alt status is in the comparisons
        }
        info += "REF_COMPARISON=" + comparisonStr;

        // rebuild line
        std::ostringstream outLine;
        for (int i = 0; i < 7; i++) {
            outLine << fields[i];
            if (i < 6)
                outLine << "\t"; // Changed from i<7 to i<6 to avoid adding an extra tab after FILTER
        }
        outLine << "\t" << info; // Add tab after FILTER, then add INFO without spaces
        // any other columns
        for (size_t i = 8; i < fields.size(); i++) {
            outLine << "\t" << fields[i];
        }
        vcfOut << outLine.str() << "\n";
    }
}

static void show_help() {
    VCFXRefComparator obj;
    char arg0[] = "VCFX_ref_comparator";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_ref_comparator", show_help))
        return 0;
    VCFXRefComparator refComp;
    return refComp.run(argc, argv);
}
