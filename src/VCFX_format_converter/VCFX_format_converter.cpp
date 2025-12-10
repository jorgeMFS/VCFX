#include "VCFX_format_converter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <getopt.h>
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

// -----------------------------------------------------------------------
// printHelp
// -----------------------------------------------------------------------
void printHelp() {
    std::cout << "VCFX_format_converter\n"
              << "Usage: VCFX_format_converter [OPTIONS]\n\n"
              << "Options:\n"
              << "  --to-bed             Convert VCF to BED format.\n"
              << "  --to-csv             Convert VCF to CSV format.\n"
              << "  -i, --input FILE     Input VCF file (uses mmap for better performance).\n"
              << "  --help, -h           Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Converts VCF files to specified formats (BED or CSV).\n\n"
              << "Example:\n"
              << "  ./VCFX_format_converter --to-bed -i input.vcf > output.bed\n"
              << "  ./VCFX_format_converter --to-csv < input.vcf > output.csv\n";
}

// -----------------------------------------------------------------------
// parseArguments
// -----------------------------------------------------------------------
bool parseArguments(int argc, char *argv[], OutputFormat &format) {
    format = OutputFormat::UNKNOWN;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--to-bed") {
            format = OutputFormat::BED;
        } else if (arg == "--to-csv") {
            format = OutputFormat::CSV;
        } else if (arg == "--help" || arg == "-h") {
            // We'll handle help outside
        }
    }
    return (format != OutputFormat::UNKNOWN);
}

// -----------------------------------------------------------------------
// convertVCFtoBED
//   - For each variant line, output a single BED line:
//     chrom, start=(pos-1 clamped to >=0), end=(pos-1 + ref.size()), name=id
// -----------------------------------------------------------------------
void convertVCFtoBED(std::istream &in, std::ostream &out) {
    std::string line;
    std::vector<std::string> fields;
    fields.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            // Skip header or empty
            continue;
        }

        vcfx::split_tabs(line, fields);

        if (fields.size() < 5) {
            // Not enough columns to parse properly
            continue;
        }

        const std::string &chrom = fields[0];
        int pos = 0;
        try {
            pos = std::stoi(fields[1]);
        } catch (...) {
            // invalid pos => skip
            continue;
        }
        const std::string &id = fields[2];
        const std::string &ref = fields[3];
        // alt = fields[4] if needed, but not used in basic bed

        // Start coordinate
        // If pos==1 => start=0
        int start = std::max(0, pos - 1);  // clamp to >=0
        int end = start + (int)ref.size(); // simplistic approach if ref is multi-base

        // Format: chrom, start, end, name
        out << chrom << "\t" << start << "\t" << end << "\t" << id << "\n";
    }
}

// -----------------------------------------------------------------------
// convertVCFtoBEDMmap - Memory-mapped version
// -----------------------------------------------------------------------
bool convertVCFtoBEDMmap(const char* filepath, std::ostream& out) {
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

    // Pre-allocated buffer for integer to string conversion
    char numBuf[32];

    while (pos < end) {
        const char* lineEnd = findNewline(pos, end);
        std::string_view line(pos, lineEnd - pos);

        // Skip empty lines and headers
        if (line.empty() || line[0] == '#') {
            pos = lineEnd + 1;
            continue;
        }

        // Extract fields: CHROM, POS, ID, REF
        std::string_view chrom = getNthField(line, 0);
        std::string_view posField = getNthField(line, 1);
        std::string_view id = getNthField(line, 2);
        std::string_view ref = getNthField(line, 3);

        if (chrom.empty() || posField.empty() || ref.empty()) {
            pos = lineEnd + 1;
            continue;
        }

        // Parse position
        int varPos = parseIntFast(posField);
        int bedStart = std::max(0, varPos - 1);
        int bedEnd = bedStart + static_cast<int>(ref.size());

        // Output: chrom, start, end, name
        outBuf.write(chrom);
        outBuf.write('\t');

        // Convert start to string
        int len = snprintf(numBuf, sizeof(numBuf), "%d", bedStart);
        outBuf.write(std::string_view(numBuf, len));
        outBuf.write('\t');

        // Convert end to string
        len = snprintf(numBuf, sizeof(numBuf), "%d", bedEnd);
        outBuf.write(std::string_view(numBuf, len));
        outBuf.write('\t');

        outBuf.write(id.empty() ? std::string_view(".") : id);
        outBuf.write('\n');

        pos = lineEnd + 1;
    }

    outBuf.flush();
    return true;
}

// -----------------------------------------------------------------------
// Helper function: CSV-escape a single field.
//   - If field has comma or quote, enclose in quotes and double the quotes.
// -----------------------------------------------------------------------
static std::string csvEscape(const std::string &field) {
    bool needQuotes = false;
    // check if it contains a comma or a quote
    if (field.find(',') != std::string::npos || field.find('"') != std::string::npos) {
        needQuotes = true;
    }
    // if it has any control chars or whitespace, you might also consider quoting
    // but let's just do commas and quotes

    if (!needQuotes) {
        return field;
    }
    // we do need quotes -> double all existing quotes
    std::string tmp;
    tmp.reserve(field.size() + 2);
    tmp.push_back('"');
    for (char c : field) {
        if (c == '"') {
            // double it by writing two quotes
            tmp.push_back('"');
            tmp.push_back('"');
        } else {
            tmp.push_back(c);
        }
    }
    tmp.push_back('"');
    return tmp;
}

// CSV-escape for string_view - writes directly to output buffer
static inline void csvEscapeToBuffer(std::string_view field, OutputBuffer& out) {
    bool needQuotes = false;
    for (char c : field) {
        if (c == ',' || c == '"') {
            needQuotes = true;
            break;
        }
    }

    if (!needQuotes) {
        out.write(field);
        return;
    }

    out.write('"');
    for (char c : field) {
        if (c == '"') {
            out.write('"');
            out.write('"');
        } else {
            out.write(c);
        }
    }
    out.write('"');
}

// -----------------------------------------------------------------------
// convertVCFtoCSV
//   - For each variant line, produce a CSV line of the same columns
//   - We do minimal CSV: if a field has comma or quote, we enclose in quotes
//     and double any quotes inside.
// -----------------------------------------------------------------------
void convertVCFtoCSV(std::istream &in, std::ostream &out) {
    std::string line;
    bool wroteHeader = false;
    std::vector<std::string> headerCols; // from #CHROM line if we want to produce a CSV header
    std::vector<std::string> hdrTokens;
    hdrTokens.reserve(16);
    std::vector<std::string> fields;
    fields.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            // Optionally parse #CHROM line to produce a CSV header.
            // If you'd like a CSV header row, you can do something like:
            if (!wroteHeader && line.rfind("#CHROM", 0) == 0) {
                // parse the columns
                std::string lineNoChr = line.substr(1); // drop '#'
                vcfx::split_tabs(lineNoChr, hdrTokens);

                // output as CSV header
                for (size_t i = 0; i < hdrTokens.size(); i++) {
                    out << csvEscape(hdrTokens[i]);
                    if (i + 1 < hdrTokens.size())
                        out << ",";
                }
                out << "\n";
                wroteHeader = true;
            }
            // Skip the actual variant lines in the header
            continue;
        }

        // This is a data line
        vcfx::split_tabs(line, fields);

        // Join them as CSV, escaping if needed
        for (size_t i = 0; i < fields.size(); i++) {
            out << csvEscape(fields[i]);
            if (i + 1 < fields.size())
                out << ",";
        }
        out << "\n";
    }
}

// -----------------------------------------------------------------------
// convertVCFtoCSVMmap - Memory-mapped version
// -----------------------------------------------------------------------
bool convertVCFtoCSVMmap(const char* filepath, std::ostream& out) {
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
    bool wroteHeader = false;

    while (pos < end) {
        const char* lineEnd = findNewline(pos, end);
        std::string_view line(pos, lineEnd - pos);

        // Skip empty lines
        if (line.empty()) {
            pos = lineEnd + 1;
            continue;
        }

        // Handle header lines
        if (line[0] == '#') {
            // Check for #CHROM header to produce CSV header
            if (!wroteHeader && line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
                // Parse columns (skip the #)
                std::string_view headerLine = line.substr(1);
                size_t fieldStart = 0;
                bool first = true;

                for (size_t i = 0; i <= headerLine.size(); ++i) {
                    if (i == headerLine.size() || headerLine[i] == '\t') {
                        std::string_view field = headerLine.substr(fieldStart, i - fieldStart);
                        if (!first) {
                            outBuf.write(',');
                        }
                        csvEscapeToBuffer(field, outBuf);
                        first = false;
                        fieldStart = i + 1;
                    }
                }
                outBuf.write('\n');
                wroteHeader = true;
            }
            pos = lineEnd + 1;
            continue;
        }

        // Data line - parse and output as CSV
        size_t fieldStart = 0;
        bool first = true;

        for (size_t i = 0; i <= line.size(); ++i) {
            if (i == line.size() || line[i] == '\t') {
                std::string_view field = line.substr(fieldStart, i - fieldStart);
                if (!first) {
                    outBuf.write(',');
                }
                csvEscapeToBuffer(field, outBuf);
                first = false;
                fieldStart = i + 1;
            }
        }
        outBuf.write('\n');

        pos = lineEnd + 1;
    }

    outBuf.flush();
    return true;
}

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_format_converter", show_help))
        return 0;

    OutputFormat format = OutputFormat::UNKNOWN;
    std::string inputFile;

    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"to-bed", no_argument, 0, 'B'},
        {"to-csv", no_argument, 0, 'C'},
        {"input", required_argument, 0, 'i'},
        {0, 0, 0, 0}
    };

    // Reset getopt
    optind = 1;

    while (true) {
        int c = getopt_long(argc, argv, "hi:BC", long_opts, nullptr);
        if (c == -1) break;

        switch (c) {
        case 'h':
            printHelp();
            return 0;
        case 'B':
            format = OutputFormat::BED;
            break;
        case 'C':
            format = OutputFormat::CSV;
            break;
        case 'i':
            inputFile = optarg;
            break;
        default:
            break;
        }
    }

    // Also check for --to-bed and --to-csv without getopt (backwards compat)
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--to-bed") {
            format = OutputFormat::BED;
        } else if (arg == "--to-csv") {
            format = OutputFormat::CSV;
        }
    }

    if (format == OutputFormat::UNKNOWN) {
        std::cerr << "No valid output format specified (--to-bed or --to-csv).\n";
        printHelp();
        return 1;
    }

    // Use mmap if input file is specified
    if (!inputFile.empty()) {
        bool success = false;
        switch (format) {
        case OutputFormat::BED:
            success = convertVCFtoBEDMmap(inputFile.c_str(), std::cout);
            break;
        case OutputFormat::CSV:
            success = convertVCFtoCSVMmap(inputFile.c_str(), std::cout);
            break;
        default:
            std::cerr << "Unsupported output format.\n";
            return 1;
        }
        if (!success) {
            std::cerr << "Error: failed to process input file: " << inputFile << "\n";
            return 1;
        }
    } else {
        // Use stdin
        switch (format) {
        case OutputFormat::BED:
            convertVCFtoBED(std::cin, std::cout);
            break;
        case OutputFormat::CSV:
            convertVCFtoCSV(std::cin, std::cout);
            break;
        default:
            std::cerr << "Unsupported output format.\n";
            return 1;
        }
    }

    return 0;
}
