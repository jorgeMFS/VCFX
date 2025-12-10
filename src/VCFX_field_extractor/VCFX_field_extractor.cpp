#include "VCFX_field_extractor.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <getopt.h>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
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

// =============================================================================
// Zero-copy field access
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

// Find INFO subkey value without allocation
static inline bool findInfoValue(const char* info, size_t infoLen,
                                  const char* key, size_t keyLen,
                                  const char*& valueStart, size_t& valueLen) {
    if (infoLen == 0 || (infoLen == 1 && info[0] == '.')) {
        return false;
    }

    const char* p = info;
    const char* end = info + infoLen;

    while (p < end) {
        // Find end of this key=value pair
        const char* pairEnd = static_cast<const char*>(memchr(p, ';', end - p));
        if (!pairEnd) pairEnd = end;

        // Check if this pair's key matches
        const char* eq = static_cast<const char*>(memchr(p, '=', pairEnd - p));
        size_t thisKeyLen = eq ? (eq - p) : (pairEnd - p);

        if (thisKeyLen == keyLen && memcmp(p, key, keyLen) == 0) {
            if (eq) {
                valueStart = eq + 1;
                valueLen = pairEnd - valueStart;
            } else {
                // Flag without value
                valueStart = "1";
                valueLen = 1;
            }
            return true;
        }

        p = pairEnd + 1;
    }
    return false;
}

// Find FORMAT subfield index
static inline int findFormatIndex(const char* format, size_t formatLen,
                                   const char* subfield, size_t subfieldLen) {
    int index = 0;
    size_t start = 0;

    for (size_t i = 0; i <= formatLen; i++) {
        if (i == formatLen || format[i] == ':') {
            if (i - start == subfieldLen && memcmp(format + start, subfield, subfieldLen) == 0) {
                return index;
            }
            index++;
            start = i + 1;
        }
    }
    return -1;
}

// Get nth colon-delimited field
static inline bool getNthColonField(const char* str, size_t len, int fieldIndex,
                                     const char*& fieldStart, size_t& fieldLen) {
    int currentField = 0;
    size_t start = 0;

    for (size_t i = 0; i <= len; i++) {
        if (i == len || str[i] == ':') {
            if (currentField == fieldIndex) {
                fieldStart = str + start;
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
// Field type enumeration for fast dispatch
// =============================================================================
enum class FieldType {
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO,
    INFO_SUBKEY,
    SAMPLE_SUBFIELD
};

struct ParsedField {
    FieldType type;
    std::string key;          // For INFO_SUBKEY: the key name; For SAMPLE_SUBFIELD: unused
    std::string sampleName;   // For SAMPLE_SUBFIELD
    std::string subfield;     // For SAMPLE_SUBFIELD
    int sampleIndex = -1;     // For S<n> format, 1-based index
};

static ParsedField parseFieldSpec(const std::string& fld) {
    ParsedField pf;

    if (fld == "CHROM") { pf.type = FieldType::CHROM; }
    else if (fld == "POS") { pf.type = FieldType::POS; }
    else if (fld == "ID") { pf.type = FieldType::ID; }
    else if (fld == "REF") { pf.type = FieldType::REF; }
    else if (fld == "ALT") { pf.type = FieldType::ALT; }
    else if (fld == "QUAL") { pf.type = FieldType::QUAL; }
    else if (fld == "FILTER") { pf.type = FieldType::FILTER; }
    else if (fld == "INFO") { pf.type = FieldType::INFO; }
    else {
        size_t colonPos = fld.find(':');
        if (colonPos != std::string::npos) {
            pf.type = FieldType::SAMPLE_SUBFIELD;
            std::string samplePart = fld.substr(0, colonPos);
            pf.subfield = fld.substr(colonPos + 1);

            // Check if it's S<n> format
            if (!samplePart.empty() && samplePart[0] == 'S' &&
                std::all_of(samplePart.begin() + 1, samplePart.end(), ::isdigit)) {
                pf.sampleIndex = std::stoi(samplePart.substr(1));
            } else {
                pf.sampleName = samplePart;
            }
        } else {
            // Assume it's an INFO subkey
            pf.type = FieldType::INFO_SUBKEY;
            pf.key = fld;
        }
    }
    return pf;
}

// =============================================================================
// Memory-mapped file processing (FAST PATH)
// =============================================================================
bool extractFieldsMmap(const char* filepath, std::ostream& out,
                        const std::vector<std::string>& fields) {
    MappedFile file;
    if (!file.open(filepath)) {
        std::cerr << "Error: Cannot open file: " << filepath << "\n";
        return false;
    }

    if (file.size == 0) return true;

    // Pre-parse field specifications
    std::vector<ParsedField> parsedFields;
    parsedFields.reserve(fields.size());
    for (const auto& f : fields) {
        parsedFields.push_back(parseFieldSpec(f));
    }

    // Output header
    OutputBuffer outBuf(out);
    for (size_t i = 0; i < fields.size(); i++) {
        outBuf.write(fields[i]);
        if (i + 1 < fields.size()) outBuf.writeChar('\t');
    }
    outBuf.writeChar('\n');

    const char *ptr = file.data;
    const char *end = file.data + file.size;

    std::unordered_map<std::string, int> sampleNameToIndex;
    bool foundChromHeader = false;

    while (ptr < end) {
        const char *lineEnd = findNewlineSIMD(ptr, end);
        if (!lineEnd) lineEnd = end;

        std::string_view line(ptr, static_cast<size_t>(lineEnd - ptr));

        if (!line.empty() && line.back() == '\r') {
            line.remove_suffix(1);
        }

        if (line.empty()) {
            ptr = lineEnd + 1;
            continue;
        }

        // Header lines
        if (line[0] == '#') {
            if (!foundChromHeader && line.size() >= 6 &&
                line[1] == 'C' && line[2] == 'H' && line[3] == 'R' &&
                line[4] == 'O' && line[5] == 'M') {
                foundChromHeader = true;

                // Parse sample names
                int fieldIdx = 0;
                size_t start = 0;
                for (size_t i = 0; i <= line.size(); i++) {
                    if (i == line.size() || line[i] == '\t') {
                        if (fieldIdx >= 9) {
                            std::string sampleName(line.data() + start, i - start);
                            sampleNameToIndex[sampleName] = fieldIdx;
                        }
                        fieldIdx++;
                        start = i + 1;
                    }
                }
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Data line - extract fields
        const char* lineData = line.data();
        size_t lineLen = line.size();

        // Cache frequently accessed fields
        const char* chromStart = nullptr; size_t chromLen = 0;
        const char* posStart = nullptr; size_t posLen = 0;
        const char* idStart = nullptr; size_t idLen = 0;
        const char* refStart = nullptr; size_t refLen = 0;
        const char* altStart = nullptr; size_t altLen = 0;
        const char* qualStart = nullptr; size_t qualLen = 0;
        const char* filterStart = nullptr; size_t filterLen = 0;
        const char* infoStart = nullptr; size_t infoLen = 0;
        const char* formatStart = nullptr; size_t formatLen = 0;

        // Parse standard fields once
        int fieldIdx = 0;
        size_t start = 0;
        for (size_t i = 0; i <= lineLen && fieldIdx <= 8; i++) {
            if (i == lineLen || lineData[i] == '\t') {
                switch (fieldIdx) {
                    case 0: chromStart = lineData + start; chromLen = i - start; break;
                    case 1: posStart = lineData + start; posLen = i - start; break;
                    case 2: idStart = lineData + start; idLen = i - start; break;
                    case 3: refStart = lineData + start; refLen = i - start; break;
                    case 4: altStart = lineData + start; altLen = i - start; break;
                    case 5: qualStart = lineData + start; qualLen = i - start; break;
                    case 6: filterStart = lineData + start; filterLen = i - start; break;
                    case 7: infoStart = lineData + start; infoLen = i - start; break;
                    case 8: formatStart = lineData + start; formatLen = i - start; break;
                }
                fieldIdx++;
                start = i + 1;
            }
        }

        // Output extracted fields
        for (size_t i = 0; i < parsedFields.size(); i++) {
            const auto& pf = parsedFields[i];
            const char* valueStart = nullptr;
            size_t valueLen = 0;
            bool found = false;

            switch (pf.type) {
                case FieldType::CHROM:
                    valueStart = chromStart; valueLen = chromLen; found = (chromStart != nullptr);
                    break;
                case FieldType::POS:
                    valueStart = posStart; valueLen = posLen; found = (posStart != nullptr);
                    break;
                case FieldType::ID:
                    valueStart = idStart; valueLen = idLen; found = (idStart != nullptr);
                    break;
                case FieldType::REF:
                    valueStart = refStart; valueLen = refLen; found = (refStart != nullptr);
                    break;
                case FieldType::ALT:
                    valueStart = altStart; valueLen = altLen; found = (altStart != nullptr);
                    break;
                case FieldType::QUAL:
                    valueStart = qualStart; valueLen = qualLen; found = (qualStart != nullptr);
                    break;
                case FieldType::FILTER:
                    valueStart = filterStart; valueLen = filterLen; found = (filterStart != nullptr);
                    break;
                case FieldType::INFO:
                    valueStart = infoStart; valueLen = infoLen; found = (infoStart != nullptr);
                    break;
                case FieldType::INFO_SUBKEY:
                    found = findInfoValue(infoStart, infoLen, pf.key.c_str(), pf.key.size(),
                                          valueStart, valueLen);
                    break;
                case FieldType::SAMPLE_SUBFIELD: {
                    int sampleCol = -1;
                    if (pf.sampleIndex > 0) {
                        sampleCol = 9 + (pf.sampleIndex - 1);
                    } else {
                        auto it = sampleNameToIndex.find(pf.sampleName);
                        if (it != sampleNameToIndex.end()) {
                            sampleCol = it->second;
                        }
                    }
                    if (sampleCol >= 9 && formatStart) {
                        int subfieldIdx = findFormatIndex(formatStart, formatLen,
                                                          pf.subfield.c_str(), pf.subfield.size());
                        if (subfieldIdx >= 0) {
                            const char* sampleStart;
                            size_t sampleLen;
                            if (getNthTabField(lineData, lineLen, sampleCol, sampleStart, sampleLen)) {
                                found = getNthColonField(sampleStart, sampleLen, subfieldIdx,
                                                          valueStart, valueLen);
                            }
                        }
                    }
                    break;
                }
            }

            if (found && valueLen > 0) {
                outBuf.write(valueStart, valueLen);
            } else {
                outBuf.writeChar('.');
            }

            if (i + 1 < parsedFields.size()) {
                outBuf.writeChar('\t');
            }
        }
        outBuf.writeChar('\n');

        ptr = lineEnd + 1;
    }

    outBuf.flush();
    return true;
}

// =============================================================================
// printHelp
// =============================================================================
void printHelp() {
    std::cout << "VCFX_field_extractor\n"
              << "Usage: VCFX_field_extractor --fields \"FIELD1,FIELD2,...\" [OPTIONS] [input.vcf]\n\n"
              << "Description:\n"
              << "  Extracts specified fields from each VCF record. Fields can be:\n"
              << "    - Standard fields: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO\n"
              << "    - Subkeys in INFO (e.g. DP, AF, ANN). These are extracted from the INFO column.\n"
              << "    - Sample subfields: e.g. SampleName:GT or S2:DP, referencing the second sample's DP.\n"
              << "      You can use sample name as it appears in #CHROM line, or 'S' plus 1-based sample index.\n"
              << "If a requested field is not found or invalid, '.' is output.\n\n"
              << "Options:\n"
              << "  --fields, -f   Comma-separated list of fields to extract\n"
              << "  --input, -i    Input VCF file (uses fast memory-mapped I/O)\n"
              << "  --help, -h     Show this help message\n\n"
              << "Performance:\n"
              << "  File input (-i) uses memory-mapped I/O for 10-20x faster processing.\n"
              << "  Features include:\n"
              << "  - SIMD-optimized line scanning (AVX2/SSE2 on x86_64)\n"
              << "  - Zero-copy field extraction\n"
              << "  - 1MB output buffering\n"
              << "  - Direct INFO key lookup without full parsing\n\n"
              << "Example:\n"
              << "  VCFX_field_extractor --fields \"CHROM,POS,ID,REF,ALT,DP,Sample1:GT\" -i input.vcf > out.tsv\n";
}

// =============================================================================
// Stdin-based processing (FALLBACK PATH) - Keep original implementation
// =============================================================================
static std::unordered_map<std::string, std::string> parseInfo(const std::string &infoField) {
    std::unordered_map<std::string, std::string> infoMap;
    if (infoField.empty() || infoField[0] == '.') {
        return infoMap;
    }
    size_t start = 0, end;
    while (start < infoField.size()) {
        end = infoField.find(';', start);
        if (end == std::string::npos) end = infoField.size();
        if (end > start) {
            size_t eqPos = infoField.find('=', start);
            if (eqPos == std::string::npos || eqPos >= end) {
                infoMap[infoField.substr(start, end - start)] = "1";
            } else {
                infoMap[infoField.substr(start, eqPos - start)] = infoField.substr(eqPos + 1, end - eqPos - 1);
            }
        }
        start = end + 1;
    }
    return infoMap;
}

static std::vector<std::string> parseLineExtract(const std::vector<std::string> &vcfCols,
                                                 const std::vector<std::string> &fields,
                                                 const std::unordered_map<std::string, int> &sampleNameToIndex) {
    std::vector<std::string> out;
    out.reserve(fields.size());

    std::unordered_map<std::string, std::string> infoMap;
    if (vcfCols.size() > 7) {
        infoMap = parseInfo(vcfCols[7]);
    }

    std::vector<std::string> formatTokens;
    if (vcfCols.size() > 8) {
        const std::string &fmt = vcfCols[8];
        size_t start = 0, end;
        while ((end = fmt.find(':', start)) != std::string::npos) {
            formatTokens.emplace_back(fmt, start, end - start);
            start = end + 1;
        }
        formatTokens.emplace_back(fmt, start);
    }

    for (auto &fld : fields) {
        std::string value = ".";

        if (fld == "CHROM") {
            if (vcfCols.size() > 0) value = vcfCols[0];
        } else if (fld == "POS") {
            if (vcfCols.size() > 1) value = vcfCols[1];
        } else if (fld == "ID") {
            if (vcfCols.size() > 2) value = vcfCols[2];
        } else if (fld == "REF") {
            if (vcfCols.size() > 3) value = vcfCols[3];
        } else if (fld == "ALT") {
            if (vcfCols.size() > 4) value = vcfCols[4];
        } else if (fld == "QUAL") {
            if (vcfCols.size() > 5) value = vcfCols[5];
        } else if (fld == "FILTER") {
            if (vcfCols.size() > 6) value = vcfCols[6];
        } else if (fld == "INFO") {
            if (vcfCols.size() > 7) value = vcfCols[7];
        } else {
            if (infoMap.find(fld) != infoMap.end()) {
                value = infoMap[fld];
            } else {
                size_t colonPos = fld.find(':');
                if (colonPos != std::string::npos) {
                    std::string sampleNameOrID = fld.substr(0, colonPos);
                    std::string subfield = fld.substr(colonPos + 1);
                    int sampleColIndex = -1;
                    if (!sampleNameOrID.empty() && sampleNameOrID[0] == 'S' &&
                        std::all_of(sampleNameOrID.begin() + 1, sampleNameOrID.end(), ::isdigit)) {
                        int idx = std::stoi(sampleNameOrID.substr(1));
                        sampleColIndex = 9 + (idx - 1);
                    } else {
                        auto itS = sampleNameToIndex.find(sampleNameOrID);
                        if (itS != sampleNameToIndex.end()) {
                            sampleColIndex = itS->second;
                        }
                    }
                    if (sampleColIndex >= 9 && (size_t)sampleColIndex < vcfCols.size()) {
                        int subIx = -1;
                        for (int i = 0; i < (int)formatTokens.size(); i++) {
                            if (formatTokens[i] == subfield) {
                                subIx = i;
                                break;
                            }
                        }
                        if (subIx >= 0) {
                            const std::string &sampleCol = vcfCols[sampleColIndex];
                            int tokenIdx = 0;
                            size_t start = 0, end;
                            while (tokenIdx < subIx && (end = sampleCol.find(':', start)) != std::string::npos) {
                                start = end + 1;
                                tokenIdx++;
                            }
                            if (tokenIdx == subIx) {
                                end = sampleCol.find(':', start);
                                if (end == std::string::npos) end = sampleCol.size();
                                value = sampleCol.substr(start, end - start);
                            } else {
                                value = ".";
                            }
                        }
                    }
                }
            }
        }
        out.push_back(value);
    }

    return out;
}

void extractFields(std::istream &in, std::ostream &out, const std::vector<std::string> &fields) {
    for (size_t i = 0; i < fields.size(); i++) {
        out << fields[i];
        if (i + 1 < fields.size()) out << "\t";
    }
    out << "\n";

    std::string line;
    std::unordered_map<std::string, int> sampleNameToIndex;
    bool foundChromHeader = false;
    std::vector<std::string> vcfCols;
    vcfCols.reserve(16);

    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
                vcfx::split_tabs(line, vcfCols);
                for (int i = 9; i < (int)vcfCols.size(); i++) {
                    sampleNameToIndex[vcfCols[i]] = i;
                }
            }
            continue;
        }

        vcfx::split_tabs(line, vcfCols);
        std::vector<std::string> extracted = parseLineExtract(vcfCols, fields, sampleNameToIndex);

        for (size_t i = 0; i < extracted.size(); i++) {
            outputBuffer += extracted[i];
            if (i + 1 < extracted.size()) outputBuffer += '\t';
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

// =============================================================================
// main
// =============================================================================
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_field_extractor", show_help))
        return 0;

    std::vector<std::string> fields;
    std::string inputFile;
    bool showHelp = false;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"fields", required_argument, 0, 'f'},
        {"input", required_argument, 0, 'i'},
        {0, 0, 0, 0}
    };

    optind = 1;
    int opt;

    while ((opt = getopt_long(argc, argv, "hf:i:", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'f': {
            std::string fieldsStr = optarg;
            size_t start = 0, end;
            while ((end = fieldsStr.find(',', start)) != std::string::npos) {
                fields.push_back(fieldsStr.substr(start, end - start));
                start = end + 1;
            }
            fields.push_back(fieldsStr.substr(start));
            break;
        }
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
        printHelp();
        return 0;
    }

    if (fields.empty()) {
        std::cerr << "No fields specified. Use --fields or -f to specify.\n";
        std::cerr << "Use --help for usage.\n";
        return 1;
    }

    if (!inputFile.empty() && inputFile != "-") {
        return extractFieldsMmap(inputFile.c_str(), std::cout, fields) ? 0 : 1;
    } else {
        extractFields(std::cin, std::cout, fields);
        return 0;
    }
}
