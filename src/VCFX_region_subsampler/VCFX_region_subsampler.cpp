#include "VCFX_region_subsampler.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
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

// Binary search for position in sorted intervals
static inline bool inRegionsBinary(const std::vector<Region>& ivs, int pos) {
    int left = 0, right = static_cast<int>(ivs.size()) - 1;
    while (left <= right) {
        int mid = (left + right) / 2;
        if (pos < ivs[mid].start) {
            right = mid - 1;
        } else if (pos > ivs[mid].end) {
            left = mid + 1;
        } else {
            return true;
        }
    }
    return false;
}

static void mergeIntervals(std::vector<Region> &ivs) {
    if (ivs.empty())
        return;
    std::sort(ivs.begin(), ivs.end(), [](const Region &a, const Region &b) { return a.start < b.start; });
    std::vector<Region> merged;
    merged.push_back(ivs[0]);
    for (size_t i = 1; i < ivs.size(); i++) {
        Region &last = merged.back();
        const Region &curr = ivs[i];
        if (curr.start <= last.end + 1) {
            // Overlapping or contiguous
            if (curr.end > last.end) {
                last.end = curr.end;
            }
        } else {
            merged.push_back(curr);
        }
    }
    ivs = merged;
}

// For a sorted list of intervals, check whether pos is in any (binary search)
static bool inRegions(const std::vector<Region> &ivs, int pos) {
    int left = 0, right = (int)ivs.size() - 1;
    while (left <= right) {
        int mid = (left + right) / 2;
        if (pos < ivs[mid].start) {
            right = mid - 1;
        } else if (pos > ivs[mid].end) {
            left = mid + 1;
        } else {
            // pos in [ivs[mid].start..ivs[mid].end]
            return true;
        }
    }
    return false;
}

int VCFXRegionSubsampler::run(int argc, char *argv[]) {
    bool showHelp = false;
    bool quiet = false;
    std::string bedFile;
    std::string inputFile;

    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"region-bed", required_argument, 0, 'b'},
        {"input", required_argument, 0, 'i'},
        {"quiet", no_argument, 0, 'q'},
        {0, 0, 0, 0}
    };

    // Parse arguments
    while (true) {
        int c = ::getopt_long(argc, argv, "hb:i:q", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 'b':
            bedFile = optarg;
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

    // If user asked for help or used an invalid short/long option
    if (showHelp) {
        displayHelp();
        return 0;
    }

    // If no region-bed was specified
    if (bedFile.empty()) {
        std::cerr << "Error: Must specify --region-bed <FILE>.\n";
        displayHelp();
        return 1;
    }

    // Load BED file => store intervals
    if (!loadRegions(bedFile, regions)) {
        std::cerr << "Error: failed to load regions from " << bedFile << "\n";
        return 1;
    }

    // Sort & merge intervals for each chrom
    sortAndMergeIntervals(regions);

    // Use mmap for file input (fast path)
    if (!inputFile.empty()) {
        if (!processVCFMmap(inputFile.c_str(), std::cout, quiet)) {
            std::cerr << "Error: failed to process input file: " << inputFile << "\n";
            return 1;
        }
    } else {
        // Process stdin VCF (fallback)
        processVCF(std::cin, std::cout);
    }
    return 0;
}

void VCFXRegionSubsampler::displayHelp() {
    std::cout << "VCFX_region_subsampler: Keep only variants whose (CHROM,POS) is in a set of regions.\n\n"
              << "Usage:\n"
              << "  VCFX_region_subsampler -b FILE -i input.vcf > out.vcf\n"
              << "  VCFX_region_subsampler --region-bed FILE < input.vcf > out.vcf\n\n"
              << "Options:\n"
              << "  -h, --help             Show help.\n"
              << "  -b, --region-bed FILE  BED file listing multiple regions.\n"
              << "  -i, --input FILE       Input VCF file (uses mmap for better performance).\n"
              << "  -q, --quiet            Suppress warnings.\n\n"
              << "Description:\n"
              << "  Reads the BED, which is <chrom> <start> <end> in 0-based. This tool converts\n"
              << "  them to 1-based [start+1 .. end]. Then merges intervals per chrom.\n"
              << "  Then only lines in the VCF that fall in those intervals for that CHROM are printed.\n\n"
              << "Example:\n"
              << "  VCFX_region_subsampler --region-bed myregions.bed -i input.vcf > out.vcf\n";
}

bool VCFXRegionSubsampler::loadRegions(const std::string &bedFilePath,
                                       std::unordered_map<std::string, std::vector<Region>> &chromRegions) {
    std::ifstream in(bedFilePath);
    if (!in.is_open()) {
        std::cerr << "Error: cannot open BED " << bedFilePath << "\n";
        return false;
    }
    std::string line;
    int lineCount = 0;

    while (true) {
        if (!std::getline(in, line))
            break;
        lineCount++;
        if (line.empty() || line[0] == '#')
            continue;

        std::stringstream ss(line);
        std::string chrom;
        int start = 0, end = 0;
        if (!(ss >> chrom >> start >> end)) {
            std::cerr << "Warning: skipping invalid bed line " << lineCount << ": " << line << "\n";
            continue;
        }
        if (start < 0)
            start = 0;
        Region r;
        r.start = start + 1; // 1-based
        r.end = end;

        if (r.end < r.start) {
            // ignore negative or zero-length intervals
            continue;
        }
        chromRegions[chrom].push_back(r);
    }
    return true;
}

void VCFXRegionSubsampler::sortAndMergeIntervals(std::unordered_map<std::string, std::vector<Region>> &chromRegions) {
    for (auto &kv : chromRegions) {
        auto &ivs = kv.second;
        std::sort(ivs.begin(), ivs.end(), [](const Region &a, const Region &b) { return a.start < b.start; });

        std::vector<Region> merged;
        merged.reserve(ivs.size());
        merged.push_back(ivs[0]);
        for (size_t i = 1; i < ivs.size(); i++) {
            Region &last = merged.back();
            const Region &curr = ivs[i];
            if (curr.start <= last.end + 1) {
                if (curr.end > last.end) {
                    last.end = curr.end;
                }
            } else {
                merged.push_back(curr);
            }
        }
        ivs = merged;
    }
}

bool VCFXRegionSubsampler::isInAnyRegion(const std::string &chrom, int pos) const {
    auto it = regions.find(chrom);
    if (it == regions.end())
        return false;
    const auto &ivs = it->second;

    int left = 0, right = (int)ivs.size() - 1;
    while (left <= right) {
        int mid = (left + right) / 2;
        if (pos < ivs[mid].start) {
            right = mid - 1;
        } else if (pos > ivs[mid].end) {
            left = mid + 1;
        } else {
            // pos is in [ivs[mid].start, ivs[mid].end]
            return true;
        }
    }
    return false;
}

bool VCFXRegionSubsampler::processVCFMmap(const char* filepath, std::ostream& out, bool quiet) {
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

    // Cache for current chromosome's regions (optimization for sorted VCFs)
    std::string currentChrom;
    const std::vector<Region>* currentRegions = nullptr;

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
            outBuf.write(line);
            outBuf.write('\n');
            if (line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
                foundChromHeader = true;
            }
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

        // Extract CHROM (field 0) and POS (field 1)
        std::string_view chromField = getNthField(line, 0);
        std::string_view posField = getNthField(line, 1);

        if (chromField.empty() || posField.empty()) {
            if (!quiet) {
                std::cerr << "Warning: line has insufficient columns => skipping.\n";
            }
            pos = lineEnd + 1;
            continue;
        }

        // Parse position
        int varPos = parseIntFast(posField);

        // Check if chromosome changed (cache lookup)
        if (currentRegions == nullptr ||
            chromField.size() != currentChrom.size() ||
            chromField != currentChrom) {
            currentChrom = std::string(chromField);
            auto it = regions.find(currentChrom);
            currentRegions = (it != regions.end()) ? &(it->second) : nullptr;
        }

        // Check if position is in any region
        bool inRegion = false;
        if (currentRegions != nullptr) {
            inRegion = inRegionsBinary(*currentRegions, varPos);
        }

        if (inRegion) {
            outBuf.write(line);
            outBuf.write('\n');
        }

        pos = lineEnd + 1;
    }

    outBuf.flush();
    return true;
}

void VCFXRegionSubsampler::processVCF(std::istream &in, std::ostream &out) {
    bool foundChromHeader = false;
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
            if (line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
            }
            continue;
        }
        if (!foundChromHeader) {
            std::cerr << "Warning: data line encountered before #CHROM => skipping.\n";
            continue;
        }

        vcfx::split_tabs(line, fields);
        if (fields.size() < 8) {
            std::cerr << "Warning: line has <8 columns => skipping.\n";
            continue;
        }

        std::string &chrom = fields[0];
        std::string &posStr = fields[1];
        int pos = 0;
        try {
            pos = std::stoi(posStr);
        } catch (...) {
            std::cerr << "Warning: invalid POS => skipping.\n";
            continue;
        }

        if (isInAnyRegion(chrom, pos)) {
            out << line << "\n";
        }
    }
}

static void show_help() {
    VCFXRegionSubsampler obj;
    char arg0[] = "VCFX_region_subsampler";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_region_subsampler", show_help))
        return 0;
    VCFXRegionSubsampler app;
    return app.run(argc, argv);
}
