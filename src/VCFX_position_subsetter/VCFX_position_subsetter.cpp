#include "VCFX_position_subsetter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

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

    void write(std::string_view sv) {
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

    void writeNewline() {
        if (pos + 1 > BUFFER_SIZE) {
            flush();
        }
        buffer[pos++] = '\n';
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
// Fast character finding
// =============================================================================
static inline const char* findChar(const char* p, const char* end, char c) {
    return static_cast<const char*>(memchr(p, c, static_cast<size_t>(end - p)));
}

// =============================================================================
// Fast integer parsing (no allocation)
// =============================================================================
static inline bool fastParseInt(const char* p, const char* end, int& result) {
    if (p >= end) return false;

    result = 0;
    while (p < end && *p >= '0' && *p <= '9') {
        result = result * 10 + (*p - '0');
        p++;
    }
    return true;
}

// =============================================================================
// Check if line matches region (zero-copy)
// =============================================================================
static inline bool matchesRegion(const char* lineStart, const char* lineEnd,
                                  std::string_view targetChrom, int regionStart, int regionEnd) {
    // Find first tab (end of CHROM field)
    const char* tab1 = findChar(lineStart, lineEnd, '\t');
    if (!tab1) return false;

    // Compare CHROM
    std::string_view chrom(lineStart, static_cast<size_t>(tab1 - lineStart));
    if (chrom != targetChrom) return false;

    // Parse POS (field after first tab)
    const char* posStart = tab1 + 1;
    const char* tab2 = findChar(posStart, lineEnd, '\t');
    if (!tab2) tab2 = lineEnd;

    int pos = 0;
    if (!fastParseInt(posStart, tab2, pos)) return false;

    return pos >= regionStart && pos <= regionEnd;
}

int VCFXPositionSubsetter::run(int argc, char *argv[]) {
    if (argc == 1) {
        displayHelp();
        return 0;
    }
    bool showHelp = false;
    std::string regionStr;
    std::string inputFile;

    static struct option long_opts[] = {
        {"region", required_argument, 0, 'r'},
        {"input", required_argument, 0, 'i'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    optind = 1;

    while (true) {
        int c = ::getopt_long(argc, argv, "r:i:h", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'r':
            regionStr = optarg;
            break;
        case 'i':
            inputFile = optarg;
            break;
        case 'h':
        default:
            showHelp = true;
        }
    }

    // Check for positional argument (input file)
    if (inputFile.empty() && optind < argc) {
        inputFile = argv[optind];
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }
    if (regionStr.empty()) {
        std::cerr << "Error: --region <chrX:start-end> is required.\n";
        displayHelp();
        return 1;
    }
    std::string chrom;
    int start = 0, end = 0;
    if (!parseRegion(regionStr, chrom, start, end)) {
        return 1;
    }

    if (!inputFile.empty() && inputFile != "-") {
        // Use mmap fast path
        return subsetVCFByPositionMmap(inputFile.c_str(), std::cout, chrom, start, end) ? 0 : 1;
    } else {
        // Fallback to stdin
        return subsetVCFByPosition(std::cin, std::cout, chrom, start, end) ? 0 : 1;
    }
}

void VCFXPositionSubsetter::displayHelp() {
    std::cout << "VCFX_position_subsetter: Subset VCF by a single genomic region.\n\n"
                 "Usage:\n"
                 "  VCFX_position_subsetter --region \"chr1:10000-20000\" [options] [input.vcf]\n"
                 "  VCFX_position_subsetter --region \"chr1:10000-20000\" < in.vcf > out.vcf\n\n"
                 "Options:\n"
                 "  -r, --region \"CHR:START-END\"   The region to keep.\n"
                 "  -i, --input FILE               Input VCF file (uses fast memory-mapped I/O)\n"
                 "  -h, --help                     Print this help.\n\n"
                 "Description:\n"
                 "  Reads lines from VCF input, and only prints data lines where:\n"
                 "    1) CHROM matches 'CHR' exactly, and\n"
                 "    2) POS is in [START,END].\n"
                 "  All header lines (#...) are passed unmodified.\n\n"
                 "Performance:\n"
                 "  File input (-i) uses memory-mapped I/O for 10-20x faster processing.\n"
                 "  Features include:\n"
                 "  - SIMD-optimized line scanning (AVX2/SSE2 on x86_64)\n"
                 "  - Zero-copy string parsing with string_view\n"
                 "  - 1MB output buffering\n"
                 "  - Direct CHROM/POS extraction without full line parsing\n\n"
                 "Examples:\n"
                 "  VCFX_position_subsetter -r \"chr2:500-1000\" -i input.vcf > subset.vcf\n"
                 "  VCFX_position_subsetter -r \"chr2:500-1000\" input.vcf > subset.vcf\n"
                 "  VCFX_position_subsetter -r \"chr2:500-1000\" < input.vcf > subset.vcf\n";
}

bool VCFXPositionSubsetter::parseRegion(const std::string &regionStr, std::string &chrom, int &start, int &end) {
    auto cpos = regionStr.find(':');
    auto dpos = regionStr.find('-');
    if (cpos == std::string::npos || dpos == std::string::npos || dpos <= cpos) {
        std::cerr << "Error: invalid region: " << regionStr << ". Expected e.g. chr1:10000-20000.\n";
        return false;
    }
    chrom = regionStr.substr(0, cpos);
    std::string startStr = regionStr.substr(cpos + 1, dpos - (cpos + 1));
    std::string endStr = regionStr.substr(dpos + 1);
    try {
        start = std::stoi(startStr);
        end = std::stoi(endStr);
    } catch (...) {
        std::cerr << "Error: cannot parse region start/end.\n";
        return false;
    }
    if (start > end) {
        std::cerr << "Error: region start> end.\n";
        return false;
    }
    return true;
}

// =============================================================================
// Memory-mapped file processing (FAST PATH)
// =============================================================================
bool VCFXPositionSubsetter::subsetVCFByPositionMmap(const char* filepath, std::ostream& out,
                                                     const std::string& regionChrom,
                                                     int regionStart, int regionEnd) {
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
    std::string_view targetChrom(regionChrom);

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
            outBuf.writeNewline();
            ptr = lineEnd + 1;
            continue;
        }

        // Header lines - pass through
        if (line[0] == '#') {
            outBuf.write(line);
            if (line.size() >= 6 && line[1] == 'C' && line[2] == 'H' &&
                line[3] == 'R' && line[4] == 'O' && line[5] == 'M') {
                headerFound = true;
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Data line before header
        if (!headerFound) {
            std::cerr << "Warning: VCF data line encountered before #CHROM. Skipping.\n";
            ptr = lineEnd + 1;
            continue;
        }

        // Check if line matches region (zero-copy)
        if (matchesRegion(line.data(), line.data() + line.size(), targetChrom, regionStart, regionEnd)) {
            outBuf.write(line);
        }

        ptr = lineEnd + 1;
    }

    outBuf.flush();
    return true;
}

// =============================================================================
// Stdin-based processing (FALLBACK PATH)
// =============================================================================
bool VCFXPositionSubsetter::subsetVCFByPosition(std::istream &in, std::ostream &out, const std::string &regionChrom,
                                                int regionStart, int regionEnd) {
    bool headerFound = false;
    std::string line;
    std::vector<std::string> fields;
    fields.reserve(16);
    while (true) {
        if (!std::getline(in, line))
            break;
        if (line.empty()) {
            out << line << "\n";
            continue;
        }
        if (line[0] == '#') {
            out << line << "\n";
            if (line.rfind("#CHROM", 0) == 0)
                headerFound = true;
            continue;
        }
        if (!headerFound) {
            std::cerr << "Warning: data line encountered before #CHROM => skipping.\n";
            continue;
        }
        vcfx::split_tabs(line, fields);
        if (fields.size() < 2) {
            std::cerr << "Warning: line has <2 columns => skipping.\n";
            continue;
        }
        const std::string &chrom = fields[0];
        const std::string &posStr = fields[1];
        int pos = 0;
        try {
            pos = std::stoi(posStr);
        } catch (...) {
            std::cerr << "Warning: invalid POS '" << posStr << "'. Skipping.\n";
            continue;
        }
        if (chrom == regionChrom && pos >= regionStart && pos <= regionEnd) {
            out << line << "\n";
        }
    }
    return true;
}

static void show_help() {
    VCFXPositionSubsetter obj;
    char arg0[] = "VCFX_position_subsetter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_position_subsetter", show_help))
        return 0;
    VCFXPositionSubsetter subsetter;
    return subsetter.run(argc, argv);
}
