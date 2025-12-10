#include "VCFX_impact_filter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>

// mmap support
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// SIMD support detection
#if defined(__x86_64__) || defined(_M_X64)
#if defined(__AVX2__)
#include <immintrin.h>
#define USE_AVX2 1
#elif defined(__SSE2__)
#include <emmintrin.h>
#define USE_SSE2 1
#endif
#endif

// =============================================================================
// Memory-mapped file wrapper (RAII)
// =============================================================================
struct MappedFile {
    const char *data = nullptr;
    size_t size = 0;
    int fd = -1;

    bool open(const char *path) {
        fd = ::open(path, O_RDONLY);
        if (fd < 0)
            return false;

        struct stat st;
        if (fstat(fd, &st) < 0) {
            close();
            return false;
        }
        size = static_cast<size_t>(st.st_size);
        if (size == 0) {
            return true;
        }

        data = static_cast<const char *>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED) {
            data = nullptr;
            close();
            return false;
        }

        madvise(const_cast<char *>(data), size, MADV_SEQUENTIAL);
        madvise(const_cast<char *>(data), size, MADV_WILLNEED);
        return true;
    }

    void close() {
        if (data && size > 0) {
            munmap(const_cast<char *>(data), size);
            data = nullptr;
        }
        if (fd >= 0) {
            ::close(fd);
            fd = -1;
        }
        size = 0;
    }

    ~MappedFile() { close(); }
};

// =============================================================================
// Output buffer for efficient writing
// =============================================================================
class OutputBuffer {
    static constexpr size_t BUFFER_SIZE = 1024 * 1024;  // 1MB buffer
    char* buffer;
    size_t pos = 0;
    std::ostream& out;

public:
    explicit OutputBuffer(std::ostream& os) : out(os) {
        buffer = new char[BUFFER_SIZE];
    }

    ~OutputBuffer() {
        flush();
        delete[] buffer;
    }

    void write(const char* data, size_t len) {
        if (pos + len > BUFFER_SIZE) {
            flush();
        }
        if (len > BUFFER_SIZE) {
            out.write(data, len);
            return;
        }
        memcpy(buffer + pos, data, len);
        pos += len;
    }

    void write(std::string_view sv) {
        write(sv.data(), sv.size());
    }

    void writeLine(std::string_view sv) {
        if (pos + sv.size() + 1 > BUFFER_SIZE) {
            flush();
        }
        if (sv.size() + 1 > BUFFER_SIZE) {
            out.write(sv.data(), sv.size());
            out.put('\n');
            return;
        }
        memcpy(buffer + pos, sv.data(), sv.size());
        pos += sv.size();
        buffer[pos++] = '\n';
    }

    void writeChar(char c) {
        if (pos + 1 > BUFFER_SIZE) {
            flush();
        }
        buffer[pos++] = c;
    }

    void flush() {
        if (pos > 0) {
            out.write(buffer, pos);
            pos = 0;
        }
    }
};

// =============================================================================
// SIMD-optimized newline detection
// =============================================================================
#if defined(USE_AVX2)
static inline const char *findNewlineSIMD(const char *p, const char *end) {
    const __m256i nl = _mm256_set1_epi8('\n');
    while (p + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(p));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, nl);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask) {
            return p + __builtin_ctz(static_cast<unsigned int>(mask));
        }
        p += 32;
    }
    return static_cast<const char *>(memchr(p, '\n', static_cast<size_t>(end - p)));
}
#elif defined(USE_SSE2)
static inline const char *findNewlineSIMD(const char *p, const char *end) {
    const __m128i nl = _mm_set1_epi8('\n');
    while (p + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i *>(p));
        __m128i cmp = _mm_cmpeq_epi8(chunk, nl);
        int mask = _mm_movemask_epi8(cmp);
        if (mask) {
            return p + __builtin_ctz(static_cast<unsigned int>(mask));
        }
        p += 16;
    }
    return static_cast<const char *>(memchr(p, '\n', static_cast<size_t>(end - p)));
}
#else
static inline const char *findNewlineSIMD(const char *p, const char *end) {
    return static_cast<const char *>(memchr(p, '\n', static_cast<size_t>(end - p)));
}
#endif

// =============================================================================
// Impact level classification
// =============================================================================
enum class ImpactLevel { UNKNOWN = 0, MODIFIER = 1, LOW = 2, MODERATE = 3, HIGH = 4 };

// Case-insensitive substring search (no allocation)
static inline bool containsIgnoreCase(const char* haystack, size_t haystackLen,
                                       const char* needle, size_t needleLen) {
    if (needleLen > haystackLen) return false;

    for (size_t i = 0; i <= haystackLen - needleLen; i++) {
        bool match = true;
        for (size_t j = 0; j < needleLen; j++) {
            char h = haystack[i + j];
            char n = needle[j];
            // Convert to uppercase
            if (h >= 'a' && h <= 'z') h -= 32;
            if (n >= 'a' && n <= 'z') n -= 32;
            if (h != n) {
                match = false;
                break;
            }
        }
        if (match) return true;
    }
    return false;
}

// Classify impact from a string (zero-copy version)
static inline ImpactLevel classifyImpactView(const char* str, size_t len) {
    if (containsIgnoreCase(str, len, "HIGH", 4)) {
        return ImpactLevel::HIGH;
    } else if (containsIgnoreCase(str, len, "MODERATE", 8)) {
        return ImpactLevel::MODERATE;
    } else if (containsIgnoreCase(str, len, "LOW", 3)) {
        return ImpactLevel::LOW;
    } else if (containsIgnoreCase(str, len, "MODIFIER", 8)) {
        return ImpactLevel::MODIFIER;
    }
    return ImpactLevel::UNKNOWN;
}

static inline bool meetsThreshold(ImpactLevel variantLevel, ImpactLevel targetLevel) {
    return static_cast<int>(variantLevel) >= static_cast<int>(targetLevel);
}

// Find IMPACT= in INFO field (case-insensitive, zero-copy)
// Returns pointer to start of value and sets valueLen
static const char* findImpactValue(const char* info, size_t infoLen, size_t& valueLen) {
    // Search for "IMPACT=" case-insensitively
    const char* p = info;
    const char* end = info + infoLen;

    while (p + 7 <= end) {
        // Check for "IMPACT=" case-insensitively
        bool match = true;
        const char* pattern = "IMPACT=";
        for (int i = 0; i < 7; i++) {
            char c = p[i];
            char pat = pattern[i];
            if (c >= 'a' && c <= 'z') c -= 32;
            if (c != pat) {
                match = false;
                break;
            }
        }

        if (match) {
            // Make sure this is at start or preceded by semicolon
            if (p == info || *(p-1) == ';') {
                const char* valueStart = p + 7;
                // Find end of value (semicolon or end of string)
                const char* valueEnd = valueStart;
                while (valueEnd < end && *valueEnd != ';') {
                    valueEnd++;
                }
                valueLen = valueEnd - valueStart;
                return valueStart;
            }
        }
        p++;
    }

    valueLen = 0;
    return nullptr;
}

// =============================================================================
// Find nth tab-delimited field (zero-copy)
// =============================================================================
static inline bool getNthTabField(const char* line, size_t lineLen, int fieldIndex,
                                   const char*& fieldStart, size_t& fieldLen) {
    int currentField = 0;
    size_t start = 0;

    for (size_t i = 0; i <= lineLen; i++) {
        if (i == lineLen || line[i] == '\t') {
            if (currentField == fieldIndex) {
                fieldStart = line + start;
                fieldLen = i - start;
                return true;
            }
            currentField++;
            start = i + 1;
        }
    }
    return false;
}

// =============================================================================
// Memory-mapped file processing (FAST PATH)
// =============================================================================
bool VCFXImpactFilter::filterByImpactMmap(const char* filepath, std::ostream& out,
                                           int targetLevelInt) {
    ImpactLevel targetLevel = static_cast<ImpactLevel>(targetLevelInt);
    MappedFile file;
    if (!file.open(filepath)) {
        std::cerr << "Error: Cannot open file: " << filepath << "\n";
        return false;
    }

    if (file.size == 0) return true;

    const char *ptr = file.data;
    const char *end = file.data + file.size;

    OutputBuffer outBuf(out);
    bool headerFound = false;
    bool wroteInfoMeta = false;

    while (ptr < end) {
        const char *lineEnd = findNewlineSIMD(ptr, end);
        if (!lineEnd) lineEnd = end;

        std::string_view line(ptr, static_cast<size_t>(lineEnd - ptr));

        // Handle Windows line endings
        if (!line.empty() && line.back() == '\r') {
            line.remove_suffix(1);
        }

        // Empty line
        if (line.empty()) {
            outBuf.writeChar('\n');
            ptr = lineEnd + 1;
            continue;
        }

        // Header lines
        if (line[0] == '#') {
            // Insert meta info line before #CHROM
            if (!wroteInfoMeta && line.size() >= 6 && line[1] == 'C' && line[2] == 'H' &&
                line[3] == 'R' && line[4] == 'O' && line[5] == 'M') {
                outBuf.write("##INFO=<ID=EXTRACTED_IMPACT,Number=1,Type=String,Description=\"Extracted from IMPACT=... in info.\">\n");
                wroteInfoMeta = true;
                headerFound = true;
            }
            outBuf.writeLine(line);
            ptr = lineEnd + 1;
            continue;
        }

        // Data line before header
        if (!headerFound) {
            std::cerr << "Error: VCF data encountered before #CHROM line.\n";
            return false;
        }

        // Parse INFO field (column 7, 0-indexed)
        const char* infoStart;
        size_t infoLen;
        if (!getNthTabField(line.data(), line.size(), 7, infoStart, infoLen)) {
            // Invalid line, skip
            ptr = lineEnd + 1;
            continue;
        }

        // Find IMPACT= in INFO
        size_t impactValueLen;
        const char* impactValue = findImpactValue(infoStart, infoLen, impactValueLen);

        std::string_view extracted;
        ImpactLevel varLevel;

        if (impactValue && impactValueLen > 0) {
            extracted = std::string_view(impactValue, impactValueLen);
            varLevel = classifyImpactView(impactValue, impactValueLen);
        } else {
            extracted = std::string_view("UNKNOWN", 7);
            varLevel = ImpactLevel::UNKNOWN;
        }

        // Check if meets threshold
        if (meetsThreshold(varLevel, targetLevel)) {
            // Output line with EXTRACTED_IMPACT appended to INFO
            // We need to reconstruct the line with modified INFO field

            // Find positions of all tabs to locate field boundaries
            const char* lineData = line.data();
            size_t lineLen = line.size();

            // Write fields 0-6 (CHROM through FILTER)
            size_t fieldStart = 0;
            int fieldIdx = 0;
            for (size_t i = 0; i <= lineLen; i++) {
                if (i == lineLen || lineData[i] == '\t') {
                    if (fieldIdx <= 6) {
                        outBuf.write(lineData + fieldStart, i - fieldStart);
                        if (i < lineLen) outBuf.writeChar('\t');
                    } else if (fieldIdx == 7) {
                        // INFO field - append EXTRACTED_IMPACT
                        if (i - fieldStart == 1 && lineData[fieldStart] == '.') {
                            // INFO is just "."
                            outBuf.write("EXTRACTED_IMPACT=");
                            outBuf.write(extracted);
                        } else {
                            outBuf.write(lineData + fieldStart, i - fieldStart);
                            outBuf.write(";EXTRACTED_IMPACT=");
                            outBuf.write(extracted);
                        }
                        // Write remaining fields (FORMAT and samples)
                        if (i < lineLen) {
                            outBuf.write(lineData + i, lineLen - i);
                        }
                        break;
                    }
                    fieldIdx++;
                    fieldStart = i + 1;
                }
            }
            outBuf.writeChar('\n');
        }

        ptr = lineEnd + 1;
    }

    outBuf.flush();
    return true;
}

// =============================================================================
// Main entry point
// =============================================================================
int VCFXImpactFilter::run(int argc, char *argv[]) {
    if (argc == 1) {
        displayHelp();
        return 0;
    }

    int opt;
    bool showHelp = false;
    std::string targetImpact;
    std::string inputFile;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"filter-impact", required_argument, 0, 'i'},
        {"input", required_argument, 0, 'I'},
        {0, 0, 0, 0}
    };

    optind = 1;

    while ((opt = getopt_long(argc, argv, "hi:I:", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'i':
            targetImpact = optarg;
            break;
        case 'I':
            inputFile = optarg;
            break;
        default:
            showHelp = true;
        }
    }

    // Check for positional argument (input file)
    if (inputFile.empty() && optind < argc) {
        inputFile = argv[optind];
    }

    if (showHelp || targetImpact.empty()) {
        displayHelp();
        return targetImpact.empty() ? 1 : 0;
    }

    // Parse target impact level
    ImpactLevel targetLevel = classifyImpactView(targetImpact.c_str(), targetImpact.size());
    if (targetLevel == ImpactLevel::UNKNOWN) {
        std::cerr << "Error: Unrecognized impact level \"" << targetImpact << "\".\n"
                  << "Must be one of HIGH, MODERATE, LOW, MODIFIER.\n";
        return 1;
    }

    if (!inputFile.empty() && inputFile != "-") {
        // Use mmap fast path
        return filterByImpactMmap(inputFile.c_str(), std::cout, static_cast<int>(targetLevel)) ? 0 : 1;
    } else {
        // Fallback to stdin
        filterByImpact(std::cin, std::cout, targetImpact);
        return 0;
    }
}

void VCFXImpactFilter::displayHelp() {
    std::cout << "VCFX_impact_filter: Filter VCF variants based on predicted impact from annotations.\n\n"
              << "Usage:\n"
              << "  VCFX_impact_filter --filter-impact <LEVEL> [options] [input.vcf]\n"
              << "  VCFX_impact_filter --filter-impact <LEVEL> < input.vcf > filtered.vcf\n\n"
              << "Options:\n"
              << "  -h, --help                   Show this help message\n"
              << "  -i, --filter-impact <LEVEL>  One of: HIGH, MODERATE, LOW, MODIFIER\n"
              << "  -I, --input FILE             Input VCF file (uses fast memory-mapped I/O)\n\n"
              << "Performance:\n"
              << "  File input (-I) uses memory-mapped I/O for 10-20x faster processing.\n"
              << "  Features include:\n"
              << "  - SIMD-optimized line scanning (AVX2/SSE2 on x86_64)\n"
              << "  - Zero-copy string parsing (no regex)\n"
              << "  - 1MB output buffering\n\n"
              << "Description:\n"
              << "  Looks in INFO for 'IMPACT=...' (case-insensitive), extracts that string,\n"
              << "  classifies it by whether it contains 'HIGH', 'MODERATE', 'LOW', or 'MODIFIER'.\n"
              << "  Then only outputs lines whose classification is >= the requested level.\n"
              << "  Also appends ';EXTRACTED_IMPACT=Value' to the INFO field.\n\n"
              << "Example:\n"
              << "  VCFX_impact_filter --filter-impact HIGH -I input.vcf > filtered.vcf\n";
}

// =============================================================================
// Stdin-based processing (FALLBACK PATH)
// =============================================================================
void VCFXImpactFilter::filterByImpact(std::istream &in, std::ostream &out, const std::string &targetImpact) {
    ImpactLevel targetLevel = classifyImpactView(targetImpact.c_str(), targetImpact.size());
    if (targetLevel == ImpactLevel::UNKNOWN) {
        std::cerr << "Error: Unrecognized impact level \"" << targetImpact << "\".\n"
                  << "Must be one of HIGH, MODERATE, LOW, MODIFIER.\n";
        return;
    }

    bool wroteHeader = false;
    bool wroteInfoMeta = false;

    std::string line;
    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);

    while (std::getline(in, line)) {
        if (line.empty()) {
            outputBuffer += '\n';
            continue;
        }

        if (line[0] == '#') {
            if (!wroteInfoMeta && line.rfind("#CHROM", 0) == 0) {
                outputBuffer += "##INFO=<ID=EXTRACTED_IMPACT,Number=1,Type=String,Description=\"Extracted from IMPACT=... in info.\">\n";
                wroteInfoMeta = true;
            }
            outputBuffer += line;
            outputBuffer += '\n';
            if (line.rfind("#CHROM", 0) == 0) {
                wroteHeader = true;
            }
            if (outputBuffer.size() > 900000) {
                out << outputBuffer;
                outputBuffer.clear();
            }
            continue;
        }

        if (!wroteHeader) {
            std::cerr << "Error: VCF data encountered before #CHROM line.\n";
            return;
        }

        // Parse INFO field
        const char* infoStart;
        size_t infoLen;
        if (!getNthTabField(line.c_str(), line.size(), 7, infoStart, infoLen)) {
            continue;
        }

        size_t impactValueLen;
        const char* impactValue = findImpactValue(infoStart, infoLen, impactValueLen);

        std::string extracted;
        ImpactLevel varLevel;

        if (impactValue && impactValueLen > 0) {
            extracted = std::string(impactValue, impactValueLen);
            varLevel = classifyImpactView(impactValue, impactValueLen);
        } else {
            extracted = "UNKNOWN";
            varLevel = ImpactLevel::UNKNOWN;
        }

        if (meetsThreshold(varLevel, targetLevel)) {
            // Reconstruct line with EXTRACTED_IMPACT
            const char* lineData = line.c_str();
            size_t lineLen = line.size();

            size_t fieldStart = 0;
            int fieldIdx = 0;
            for (size_t i = 0; i <= lineLen; i++) {
                if (i == lineLen || lineData[i] == '\t') {
                    if (fieldIdx <= 6) {
                        outputBuffer.append(lineData + fieldStart, i - fieldStart);
                        if (i < lineLen) outputBuffer += '\t';
                    } else if (fieldIdx == 7) {
                        if (i - fieldStart == 1 && lineData[fieldStart] == '.') {
                            outputBuffer += "EXTRACTED_IMPACT=";
                            outputBuffer += extracted;
                        } else {
                            outputBuffer.append(lineData + fieldStart, i - fieldStart);
                            outputBuffer += ";EXTRACTED_IMPACT=";
                            outputBuffer += extracted;
                        }
                        if (i < lineLen) {
                            outputBuffer.append(lineData + i, lineLen - i);
                        }
                        break;
                    }
                    fieldIdx++;
                    fieldStart = i + 1;
                }
            }
            outputBuffer += '\n';

            if (outputBuffer.size() > 900000) {
                out << outputBuffer;
                outputBuffer.clear();
            }
        }
    }

    if (!outputBuffer.empty()) {
        out << outputBuffer;
    }
}

static void show_help() {
    VCFXImpactFilter obj;
    char arg0[] = "VCFX_impact_filter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_impact_filter", show_help))
        return 0;
    VCFXImpactFilter filt;
    return filt.run(argc, argv);
}
