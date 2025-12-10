#include "VCFX_population_filter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#ifdef __unix__
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __AVX2__
#include <immintrin.h>
#endif

int VCFXPopulationFilter::run(int argc, char *argv[]) {
    if (argc == 1) {
        displayHelp();
        return 0;
    }
    bool showHelp = false;
    std::string populationTag;
    std::string popMapFile;
    const char* inputFile = nullptr;
    bool quiet = false;

    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"population", required_argument, 0, 'p'},
        {"pop-map", required_argument, 0, 'm'},
        {"input", required_argument, 0, 'i'},
        {"quiet", no_argument, 0, 'q'},
        {0, 0, 0, 0}
    };

    optind = 1;
    while (true) {
        int c = ::getopt_long(argc, argv, "hp:m:i:q", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 'p':
            populationTag = optarg;
            break;
        case 'm':
            popMapFile = optarg;
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
    if (populationTag.empty() || popMapFile.empty()) {
        std::cerr << "Error: Must specify --population <TAG> and --pop-map <file>.\n";
        displayHelp();
        return 1;
    }

    std::unordered_set<std::string> samplesToInclude;
    if (!loadPopulationMap(popMapFile, populationTag, samplesToInclude)) {
        std::cerr << "Error: Unable to load or parse pop map.\n";
        return 1;
    }
    if (samplesToInclude.empty() && !quiet) {
        std::cerr << "Warning: No samples found for population tag: " << populationTag << "\n";
    }

    if (inputFile) {
        return filterPopulationMmap(inputFile, std::cout, samplesToInclude, populationTag, quiet) ? 0 : 1;
    } else {
        filterPopulation(std::cin, std::cout, samplesToInclude, populationTag);
        return 0;
    }
}

void VCFXPopulationFilter::displayHelp() {
    std::cout << "VCFX_population_filter: Subset VCF to samples in specified population.\n\n"
                 "Usage:\n"
                 "  VCFX_population_filter [options] < input.vcf > output.vcf\n"
                 "  VCFX_population_filter -p TAG -m pops.txt -i input.vcf > output.vcf\n\n"
                 "Options:\n"
                 "  -h, --help               Print this help.\n"
                 "  -p, --population <TAG>   Population tag to keep (e.g. 'EUR','AFR', etc.)\n"
                 "  -m, --pop-map <FILE>     Tab-delimited file: 'SampleName <tab> Population'\n"
                 "  -i, --input FILE         Input VCF file (default: stdin)\n"
                 "  -q, --quiet              Suppress warnings\n\n"
                 "Description:\n"
                 "  Reads the pop map, finds samples that match the chosen population.\n"
                 "  Then reads the VCF from stdin and prints lines with only those sample columns.\n"
                 "  If a sample is not in that population, it's dropped from the #CHROM header and data columns.\n\n"
                 "Example:\n"
                 "  VCFX_population_filter --population AFR --pop-map pops.txt < input.vcf > out.vcf\n"
                 "  VCFX_population_filter -p AFR -m pops.txt -i input.vcf > out.vcf\n";
}

bool VCFXPopulationFilter::loadPopulationMap(const std::string &popMapFile, const std::string &popTag,
                                             std::unordered_set<std::string> &samplesToInclude) {
    std::ifstream f(popMapFile);
    if (!f.is_open()) {
        std::cerr << "Error: cannot open " << popMapFile << "\n";
        return false;
    }
    std::string line;
    while (true) {
        if (!std::getline(f, line))
            break;
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::string samp, pop;
        if (!(ss >> samp >> pop)) {
            std::cerr << "Warning: popmap line invalid: " << line << "\n";
            continue;
        }
        if (pop == popTag) {
            samplesToInclude.insert(samp);
        }
    }
    return true;
}

void VCFXPopulationFilter::filterPopulation(std::istream &in, std::ostream &out,
                                            const std::unordered_set<std::string> &samplesToInclude,
                                            const std::string &popTag) {
    bool foundChromLine = false;
    std::string line;
    std::vector<std::string> finalHeader;
    std::vector<int> sampleIndices;
    std::vector<std::string> fields;
    fields.reserve(16);
    std::vector<std::string> columns;
    columns.reserve(16);

    while (true) {
        if (!std::getline(in, line))
            break;
        if (line.empty()) {
            out << line << "\n";
            continue;
        }
        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                foundChromLine = true;
                vcfx::split_tabs(line, fields);
                // fields[0..8] => fixed, fields[9..] => samples
                sampleIndices.clear();
                for (size_t i = 9; i < fields.size(); i++) {
                    if (samplesToInclude.count(fields[i]) > 0) {
                        sampleIndices.push_back((int)i);
                    }
                }
                // build new #CHROM line
                std::ostringstream newChrom;
                for (int i = 0; i < 9; i++) {
                    newChrom << fields[i];
                    if (i < 8 || !sampleIndices.empty())
                        newChrom << "\t";
                }
                for (size_t i = 0; i < sampleIndices.size(); i++) {
                    newChrom << fields[sampleIndices[i]];
                    if (i + 1 < sampleIndices.size())
                        newChrom << "\t";
                }
                out << newChrom.str() << "\n";
            } else {
                out << line << "\n";
            }
            continue;
        }
        if (!foundChromLine) {
            std::cerr << "Warning: data line before #CHROM => skipping.\n";
            continue;
        }
        vcfx::split_tabs(line, columns);
        if (columns.size() < 9) {
            std::cerr << "Warning: line with fewer than 9 columns => skipping.\n";
            continue;
        }
        // build new line
        std::ostringstream newLine;
        for (int i = 0; i < 9; i++) {
            newLine << columns[i];
            if (i < 8 || !sampleIndices.empty())
                newLine << "\t";
        }
        for (size_t i = 0; i < sampleIndices.size(); i++) {
            newLine << columns[sampleIndices[i]];
            if (i + 1 < sampleIndices.size())
                newLine << "\t";
        }
        out << newLine.str() << "\n";
    }
    if (!foundChromLine) {
        std::cerr << "Error: No #CHROM header found in VCF.\n";
    }
}

// ============================================================================
// MMAP-based high-performance implementation
// ============================================================================

#ifdef __unix__

namespace {

// RAII wrapper for memory-mapped files
struct MappedFile {
    const char* data = nullptr;
    size_t size = 0;
    int fd = -1;

    MappedFile() = default;
    ~MappedFile() {
        if (data && data != MAP_FAILED) {
            munmap(const_cast<char*>(data), size);
        }
        if (fd >= 0) {
            close(fd);
        }
    }
    MappedFile(const MappedFile&) = delete;
    MappedFile& operator=(const MappedFile&) = delete;

    bool open(const char* filepath) {
        fd = ::open(filepath, O_RDONLY);
        if (fd < 0) return false;

        struct stat st;
        if (fstat(fd, &st) < 0) return false;
        size = static_cast<size_t>(st.st_size);
        if (size == 0) return true;

        data = static_cast<const char*>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED) {
            data = nullptr;
            return false;
        }

        madvise(const_cast<char*>(data), size, MADV_SEQUENTIAL | MADV_WILLNEED);
        return true;
    }
};

// Output buffer for efficient writes
class OutputBuffer {
public:
    explicit OutputBuffer(std::ostream& os, size_t bufSize = 1024 * 1024)
        : out_(os), buffer_(bufSize), pos_(0) {}

    ~OutputBuffer() { flush(); }

    void append(std::string_view sv) {
        if (pos_ + sv.size() > buffer_.size()) {
            flush();
        }
        if (sv.size() > buffer_.size()) {
            out_.write(sv.data(), sv.size());
        } else {
            memcpy(buffer_.data() + pos_, sv.data(), sv.size());
            pos_ += sv.size();
        }
    }

    void append(char c) {
        if (pos_ >= buffer_.size()) {
            flush();
        }
        buffer_[pos_++] = c;
    }

    void flush() {
        if (pos_ > 0) {
            out_.write(buffer_.data(), pos_);
            pos_ = 0;
        }
    }

private:
    std::ostream& out_;
    std::vector<char> buffer_;
    size_t pos_;
};

// SIMD-accelerated newline search
inline const char* findNewline(const char* start, const char* end) {
#ifdef __AVX2__
    const __m256i newline = _mm256_set1_epi8('\n');
    while (start + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(start));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, newline);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask != 0) {
            return start + __builtin_ctz(mask);
        }
        start += 32;
    }
#endif
#ifdef __SSE2__
    const __m128i newline = _mm_set1_epi8('\n');
    while (start + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i*>(start));
        __m128i cmp = _mm_cmpeq_epi8(chunk, newline);
        int mask = _mm_movemask_epi8(cmp);
        if (mask != 0) {
            return start + __builtin_ctz(mask);
        }
        start += 16;
    }
#endif
    const char* p = static_cast<const char*>(memchr(start, '\n', end - start));
    return p ? p : end;
}

// Get nth tab-separated field using zero-copy
inline std::string_view getNthField(std::string_view line, size_t n) {
    const char* p = line.data();
    const char* end = p + line.size();
    size_t idx = 0;

    while (idx < n && p < end) {
        const char* tab = static_cast<const char*>(memchr(p, '\t', end - p));
        if (!tab) return std::string_view();
        p = tab + 1;
        ++idx;
    }

    if (p >= end) return std::string_view();

    const char* tab = static_cast<const char*>(memchr(p, '\t', end - p));
    return std::string_view(p, (tab ? tab : end) - p);
}

} // anonymous namespace

bool VCFXPopulationFilter::filterPopulationMmap(const char* filepath, std::ostream& out,
                                                 const std::unordered_set<std::string>& samplesToInclude,
                                                 const std::string& popTag, bool quiet) {
    MappedFile file;
    if (!file.open(filepath)) {
        if (!quiet) {
            std::cerr << "Error: Cannot open file: " << filepath << "\n";
        }
        return false;
    }

    if (file.size == 0) {
        return true;
    }

    const char* data = file.data;
    const char* end = data + file.size;
    const char* pos = data;

    OutputBuffer buf(out);
    bool foundChromLine = false;
    std::vector<size_t> sampleIndices;

    while (pos < end) {
        const char* lineEnd = findNewline(pos, end);
        std::string_view line(pos, lineEnd - pos);
        pos = (lineEnd < end) ? lineEnd + 1 : end;

        if (line.empty()) {
            buf.append('\n');
            continue;
        }

        if (line[0] == '#') {
            if (line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
                foundChromLine = true;

                // Parse header to find sample indices
                sampleIndices.clear();
                const char* p = line.data();
                const char* lend = p + line.size();
                size_t colIdx = 0;

                while (p < lend) {
                    const char* tab = static_cast<const char*>(memchr(p, '\t', lend - p));
                    const char* fieldEnd = tab ? tab : lend;

                    if (colIdx >= 9) {
                        std::string_view sampleName(p, fieldEnd - p);
                        std::string sampleStr(sampleName);
                        if (samplesToInclude.count(sampleStr) > 0) {
                            sampleIndices.push_back(colIdx);
                        }
                    }

                    if (!tab) break;
                    p = tab + 1;
                    ++colIdx;
                }

                // Output new header line
                p = line.data();
                colIdx = 0;
                while (p < lend && colIdx < 9) {
                    const char* tab = static_cast<const char*>(memchr(p, '\t', lend - p));
                    const char* fieldEnd = tab ? tab : lend;
                    buf.append(std::string_view(p, fieldEnd - p));
                    if (colIdx < 8 || !sampleIndices.empty()) {
                        buf.append('\t');
                    }
                    if (!tab) break;
                    p = tab + 1;
                    ++colIdx;
                }

                // Output selected sample names
                for (size_t i = 0; i < sampleIndices.size(); ++i) {
                    std::string_view field = getNthField(line, sampleIndices[i]);
                    buf.append(field);
                    if (i + 1 < sampleIndices.size()) {
                        buf.append('\t');
                    }
                }
                buf.append('\n');
            } else {
                buf.append(line);
                buf.append('\n');
            }
            continue;
        }

        if (!foundChromLine) {
            if (!quiet) {
                std::cerr << "Warning: data line before #CHROM => skipping.\n";
            }
            continue;
        }

        // Data line - output first 9 columns + selected samples
        const char* p = line.data();
        const char* lend = p + line.size();
        size_t colIdx = 0;

        // Output first 9 columns
        while (p < lend && colIdx < 9) {
            const char* tab = static_cast<const char*>(memchr(p, '\t', lend - p));
            const char* fieldEnd = tab ? tab : lend;
            buf.append(std::string_view(p, fieldEnd - p));
            if (colIdx < 8 || !sampleIndices.empty()) {
                buf.append('\t');
            }
            if (!tab) break;
            p = tab + 1;
            ++colIdx;
        }

        // Output selected sample columns
        for (size_t i = 0; i < sampleIndices.size(); ++i) {
            std::string_view field = getNthField(line, sampleIndices[i]);
            buf.append(field);
            if (i + 1 < sampleIndices.size()) {
                buf.append('\t');
            }
        }
        buf.append('\n');
    }

    if (!foundChromLine && !quiet) {
        std::cerr << "Error: No #CHROM header found in VCF.\n";
    }

    return true;
}

#else
// Non-Unix fallback
bool VCFXPopulationFilter::filterPopulationMmap(const char* filepath, std::ostream& out,
                                                 const std::unordered_set<std::string>& samplesToInclude,
                                                 const std::string& popTag, bool quiet) {
    (void)quiet;
    std::ifstream in(filepath);
    if (!in) {
        std::cerr << "Error: Cannot open file: " << filepath << "\n";
        return false;
    }
    filterPopulation(in, out, samplesToInclude, popTag);
    return true;
}
#endif

static void show_help() {
    VCFXPopulationFilter obj;
    char arg0[] = "VCFX_population_filter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_population_filter", show_help))
        return 0;
    VCFXPopulationFilter pf;
    return pf.run(argc, argv);
}
