#include "vcfx_core.h"
#include <algorithm>
#include <cctype>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <cstring>

namespace vcfx {

std::string trim(const std::string& str) {
    auto first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) {
        return "";
    }
    auto last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> result;
    std::istringstream iss(str);
    std::string item;
    while (std::getline(iss, item, delimiter)) {
        result.push_back(item);
    }
    return result;
}

bool flag_present(int argc, char* argv[], const char* long_flag,
                  const char* short_flag) {
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], long_flag) == 0 ||
            (short_flag && std::strcmp(argv[i], short_flag) == 0)) {
            return true;
        }
    }
    return false;
}

void print_error(const std::string& msg, std::ostream& os) {
    os << "Error: " << msg << '\n';
}

void print_version(const std::string& tool, const std::string& version,
                   std::ostream& os) {
    os << tool << " version " << version << '\n';
}

// ------------------------------------------------------------
// Internal helper: decompress gzip/BGZF data from 'in' into 'out'
// ------------------------------------------------------------
static bool decompress_gzip_stream(std::istream& in, std::string& out) {
    constexpr int CHUNK = 16384;
    char inBuf[CHUNK];
    char outBuf[CHUNK];

    z_stream strm;
    std::memset(&strm, 0, sizeof(strm));
    if (inflateInit2(&strm, 15 + 32) != Z_OK) {
        return false;
    }

    int ret = Z_OK;
    do {
        in.read(inBuf, CHUNK);
        strm.avail_in = static_cast<uInt>(in.gcount());
        if (strm.avail_in == 0 && in.eof()) {
            break;
        }
        strm.next_in = reinterpret_cast<Bytef*>(inBuf);

        do {
            strm.avail_out = CHUNK;
            strm.next_out = reinterpret_cast<Bytef*>(outBuf);
            ret = inflate(&strm, Z_NO_FLUSH);
            if (ret == Z_STREAM_ERROR || ret == Z_NEED_DICT ||
                ret == Z_DATA_ERROR || ret == Z_MEM_ERROR) {
                inflateEnd(&strm);
                return false;
            }
            size_t have = CHUNK - strm.avail_out;
            if (have > 0) {
                out.append(outBuf, have);
            }
        } while (strm.avail_out == 0);
    } while (ret != Z_STREAM_END);

    inflateEnd(&strm);
    return ret == Z_STREAM_END;
}

// ------------------------------------------------------------
// Detect gzip magic numbers on a stream without consuming them
// ------------------------------------------------------------
static bool stream_has_gzip_magic(std::istream& in) {
    int c1 = in.get();
    if (c1 == EOF) {
        return false;
    }
    int c2 = in.get();
    if (c2 == EOF) {
        in.unget();
        return false;
    }
    bool isGz = (static_cast<unsigned char>(c1) == 0x1f &&
                 static_cast<unsigned char>(c2) == 0x8b);
    in.putback(static_cast<char>(c2));
    in.putback(static_cast<char>(c1));
    return isGz;
}

bool read_maybe_compressed(std::istream& in, std::string& out) {
    out.clear();
    if (stream_has_gzip_magic(in)) {
        return decompress_gzip_stream(in, out);
    }
    std::ostringstream oss;
    oss << in.rdbuf();
    out = oss.str();
    return true;
}

bool read_file_maybe_compressed(const std::string& path, std::string& out) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }
    bool isGz = false;
    if (path.size() >= 3 &&
        (path.compare(path.size() - 3, 3, ".gz") == 0)) {
        isGz = true;
    } else if (path.size() >= 4 &&
               (path.compare(path.size() - 4, 4, ".bgz") == 0)) {
        isGz = true;
    } else if (path.size() >= 5 &&
               (path.compare(path.size() - 5, 5, ".bgzf") == 0)) {
        isGz = true;
    }
    if (isGz || stream_has_gzip_magic(file)) {
        return decompress_gzip_stream(file, out);
    }
    std::ostringstream oss;
    oss << file.rdbuf();
    out = oss.str();
    return true;
}

bool handle_common_flags(int argc, char* argv[],
                         const std::string& tool,
                         void (*print_help)(),
                         std::ostream& os) {
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--version") == 0 ||
            std::strcmp(argv[i], "-v") == 0) {
            print_version(tool, get_version(), os);
            return true;
        }
        if (std::strcmp(argv[i], "--help") == 0 ||
            std::strcmp(argv[i], "-h") == 0) {
            if (print_help) {
                print_help();
            }
            return true;
        }
    }
    return false;
}

}  // namespace vcfx
