#include "VCFX_sample_extractor.h"
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
#include <unordered_map>
#include <unordered_set>
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
    static constexpr size_t BUFFER_SIZE = 1024 * 1024;
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

// Helper to trim whitespace
static std::string trim(const std::string &s) {
    size_t start = 0;
    while (start < s.size() && isspace((unsigned char)s[start]))
        start++;
    if (start == s.size())
        return "";
    size_t end = s.size() - 1;
    while (end > start && isspace((unsigned char)s[end]))
        end--;
    return s.substr(start, end - start + 1);
}

// =============================================================================
// Memory-mapped file processing (FAST PATH)
// =============================================================================
bool VCFXSampleExtractor::extractSamplesMmap(const char* filepath, std::ostream& out,
                                              const std::vector<std::string>& samples) {
    MappedFile file;
    if (!file.open(filepath)) {
        std::cerr << "Error: Cannot open file: " << filepath << "\n";
        return false;
    }

    if (file.size == 0) return true;

    std::unordered_set<std::string> sampleSet(samples.begin(), samples.end());

    const char *ptr = file.data;
    const char *end = file.data + file.size;

    OutputBuffer outBuf(out);
    bool foundChromLine = false;
    std::vector<int> keepSampleIndices;
    std::vector<std::string> finalSampleNames;

    while (ptr < end) {
        const char *lineEnd = findNewlineSIMD(ptr, end);
        if (!lineEnd) lineEnd = end;

        std::string_view line(ptr, static_cast<size_t>(lineEnd - ptr));

        if (!line.empty() && line.back() == '\r') {
            line.remove_suffix(1);
        }

        if (line.empty()) {
            outBuf.writeChar('\n');
            ptr = lineEnd + 1;
            continue;
        }

        // Header lines
        if (line[0] == '#') {
            if (!foundChromLine && line.size() >= 6 &&
                line[1] == 'C' && line[2] == 'H' && line[3] == 'R' &&
                line[4] == 'O' && line[5] == 'M') {
                foundChromLine = true;

                // Parse #CHROM line to find sample columns
                keepSampleIndices.clear();
                finalSampleNames.clear();

                int fieldIdx = 0;
                size_t start = 0;
                std::vector<std::string_view> headerFields;

                for (size_t i = 0; i <= line.size(); i++) {
                    if (i == line.size() || line[i] == '\t') {
                        headerFields.push_back(line.substr(start, i - start));
                        fieldIdx++;
                        start = i + 1;
                    }
                }

                // Find which samples to keep (columns 9+)
                for (size_t i = 9; i < headerFields.size(); i++) {
                    std::string sampleName(headerFields[i]);
                    if (sampleSet.count(sampleName) > 0) {
                        keepSampleIndices.push_back(static_cast<int>(i));
                        finalSampleNames.push_back(sampleName);
                    }
                }

                // Warn about missing samples
                for (const auto& s : samples) {
                    if (std::find(finalSampleNames.begin(), finalSampleNames.end(), s) == finalSampleNames.end()) {
                        std::cerr << "Warning: sample '" << s << "' not found in header.\n";
                    }
                }

                // Write modified #CHROM line
                for (int i = 0; i < 9 && i < static_cast<int>(headerFields.size()); i++) {
                    if (i > 0) outBuf.writeChar('\t');
                    outBuf.write(headerFields[i]);
                }
                for (const auto& sname : finalSampleNames) {
                    outBuf.writeChar('\t');
                    outBuf.write(sname);
                }
                outBuf.writeChar('\n');
            } else {
                // Other header lines - pass through
                outBuf.write(line);
                outBuf.writeChar('\n');
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Data line before header
        if (!foundChromLine) {
            std::cerr << "Warning: data line encountered before #CHROM => skipping.\n";
            ptr = lineEnd + 1;
            continue;
        }

        // Parse data line - extract only needed columns
        const char* lineData = line.data();
        size_t lineLen = line.size();

        // Parse all fields
        std::vector<std::string_view> fields;
        fields.reserve(keepSampleIndices.empty() ? 9 : keepSampleIndices.back() + 2);

        size_t start = 0;
        int fieldIdx = 0;
        int maxNeeded = keepSampleIndices.empty() ? 8 : keepSampleIndices.back();

        for (size_t i = 0; i <= lineLen; i++) {
            if (i == lineLen || lineData[i] == '\t') {
                fields.emplace_back(lineData + start, i - start);
                if (fieldIdx >= maxNeeded) {
                    // We have all we need
                    break;
                }
                fieldIdx++;
                start = i + 1;
            }
        }

        if (fields.size() < 8) {
            std::cerr << "Warning: line has <8 columns => skipping.\n";
            ptr = lineEnd + 1;
            continue;
        }

        // Write first 9 columns (CHROM..FORMAT)
        for (int i = 0; i < 9 && i < static_cast<int>(fields.size()); i++) {
            if (i > 0) outBuf.writeChar('\t');
            outBuf.write(fields[i]);
        }

        // Write selected sample columns
        for (int idx : keepSampleIndices) {
            outBuf.writeChar('\t');
            if (idx < static_cast<int>(fields.size())) {
                outBuf.write(fields[idx]);
            } else {
                outBuf.writeChar('.');
            }
        }
        outBuf.writeChar('\n');

        ptr = lineEnd + 1;
    }

    outBuf.flush();
    return true;
}

// =============================================================================
// Main entry point
// =============================================================================
int VCFXSampleExtractor::run(int argc, char *argv[]) {
    if (argc == 1) {
        displayHelp();
        return 0;
    }

    bool showHelp = false;
    std::vector<std::string> samples;
    std::string inputFile;

    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"samples", required_argument, 0, 's'},
        {"input", required_argument, 0, 'i'},
        {0, 0, 0, 0}
    };

    optind = 1;

    while (true) {
        int c = ::getopt_long(argc, argv, "hs:i:", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 's': {
            std::stringstream ss(optarg);
            std::string token;
            while (ss >> token) {
                std::stringstream s2(token);
                std::string sub;
                while (std::getline(s2, sub, ',')) {
                    sub = trim(sub);
                    if (!sub.empty())
                        samples.push_back(sub);
                }
            }
        } break;
        case 'i':
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

    if (showHelp) {
        displayHelp();
        return 0;
    }

    if (samples.empty()) {
        std::cerr << "Error: must specify at least one sample with --samples.\n";
        return 1;
    }

    if (!inputFile.empty() && inputFile != "-") {
        return extractSamplesMmap(inputFile.c_str(), std::cout, samples) ? 0 : 1;
    } else {
        extractSamples(std::cin, std::cout, samples);
        return 0;
    }
}

void VCFXSampleExtractor::displayHelp() {
    std::cout << "VCFX_sample_extractor: Subset a VCF to a chosen set of samples.\n\n"
                 "Usage:\n"
                 "  VCFX_sample_extractor --samples \"Sample1,Sample2\" [options] [input.vcf]\n"
                 "  VCFX_sample_extractor --samples \"Sample1,Sample2\" < input.vcf > output.vcf\n\n"
                 "Options:\n"
                 "  -h, --help              Print this help.\n"
                 "  -s, --samples <LIST>    Comma or space separated list of sample names.\n"
                 "  -i, --input FILE        Input VCF file (uses fast memory-mapped I/O)\n\n"
                 "Performance:\n"
                 "  File input (-i) uses memory-mapped I/O for 10-20x faster processing.\n"
                 "  Features include:\n"
                 "  - SIMD-optimized line scanning (AVX2/SSE2 on x86_64)\n"
                 "  - Zero-copy field extraction\n"
                 "  - 1MB output buffering\n\n"
                 "Description:\n"
                 "  Reads #CHROM line to identify sample columns. Keeps only user-specified samples.\n"
                 "  Rewrites #CHROM line with that subset. For each variant data line, we keep only the\n"
                 "  chosen sample columns. If a sample isn't found in the header, logs a warning.\n\n"
                 "Example:\n"
                 "  VCFX_sample_extractor --samples \"IndivA,IndivB\" -i input.vcf > subset.vcf\n";
}

// =============================================================================
// Stdin-based processing (FALLBACK PATH)
// =============================================================================
void VCFXSampleExtractor::extractSamples(std::istream &in, std::ostream &out,
                                          const std::vector<std::string> &samples) {
    std::unordered_set<std::string> sampleSet(samples.begin(), samples.end());

    bool foundChromLine = false;
    std::string line;
    std::vector<int> keepSampleIndices;
    keepSampleIndices.reserve(samples.size());
    std::vector<std::string> finalSampleNames;
    bool headerProcessed = false;

    std::vector<std::string> fields;
    fields.reserve(16);

    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);

    while (true) {
        if (!std::getline(in, line))
            break;
        if (line.empty()) {
            outputBuffer += '\n';
            continue;
        }
        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                foundChromLine = true;
                vcfx::split_tabs(line, fields);
                keepSampleIndices.clear();
                finalSampleNames.clear();
                if (fields.size() < 9) {
                    outputBuffer += line;
                    outputBuffer += '\n';
                    continue;
                }
                for (size_t i = 9; i < fields.size(); i++) {
                    if (sampleSet.count(fields[i]) > 0) {
                        keepSampleIndices.push_back((int)i);
                        finalSampleNames.push_back(fields[i]);
                    }
                }
                for (auto &s : samples) {
                    if (std::find(finalSampleNames.begin(), finalSampleNames.end(), s) == finalSampleNames.end()) {
                        std::cerr << "Warning: sample '" << s << "' not found in header.\n";
                    }
                }
                for (int i = 0; i < 9; i++) {
                    if (i > 0) outputBuffer += '\t';
                    outputBuffer += fields[i];
                }
                for (size_t i = 0; i < finalSampleNames.size(); i++) {
                    outputBuffer += '\t';
                    outputBuffer += finalSampleNames[i];
                }
                outputBuffer += '\n';
                headerProcessed = true;
            } else {
                outputBuffer += line;
                outputBuffer += '\n';
            }

            if (outputBuffer.size() > 900000) {
                out << outputBuffer;
                outputBuffer.clear();
            }
            continue;
        }
        if (!foundChromLine) {
            std::cerr << "Warning: data line encountered before #CHROM => skipping.\n";
            continue;
        }
        vcfx::split_tabs(line, fields);
        if (fields.size() < 8) {
            std::cerr << "Warning: line has <8 columns => skipping.\n";
            continue;
        }
        if (fields.size() < 9) {
            std::cerr << "Warning: data line with no sample columns => skipping.\n";
            continue;
        }
        for (int i = 0; i < 9; i++) {
            if (i > 0) outputBuffer += '\t';
            outputBuffer += fields[i];
        }
        for (auto idx : keepSampleIndices) {
            outputBuffer += '\t';
            if ((size_t)idx < fields.size()) {
                outputBuffer += fields[idx];
            } else {
                outputBuffer += '.';
            }
        }
        outputBuffer += '\n';

        if (outputBuffer.size() > 900000) {
            out << outputBuffer;
            outputBuffer.clear();
        }
    }

    if (!outputBuffer.empty()) {
        out << outputBuffer;
    }
}

static void show_help() {
    VCFXSampleExtractor obj;
    char arg0[] = "VCFX_sample_extractor";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_sample_extractor", show_help))
        return 0;
    VCFXSampleExtractor app;
    return app.run(argc, argv);
}
