#include "VCFX_compressor.h"
#include <zlib.h>
#include <cstring>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_compressor\n"
              << "Usage: VCFX_compressor [OPTIONS]\n\n"
              << "Options:\n"
              << "  --compress, -c         Compress the input VCF file.\n"
              << "  --decompress, -d       Decompress the input VCF file.\n"
              << "  --help, -h             Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Compresses or decompresses VCF files using gzip compression.\n\n"
              << "Examples:\n"
              << "  Compress:\n"
              << "    ./VCFX_compressor --compress < input.vcf > output.vcf.gz\n\n"
              << "  Decompress:\n"
              << "    ./VCFX_compressor --decompress < input.vcf.gz > output.vcf\n";
}

// Function to perform compression or decompression
bool compressDecompressVCF(std::istream& in, std::ostream& out, bool compress) {
    const int bufferSize = 16384;
    char buffer[bufferSize];

    if (compress) {
        // Initialize zlib for compression
        z_stream strm;
        std::memset(&strm, 0, sizeof(strm));
        if (deflateInit(&strm, Z_DEFAULT_COMPRESSION) != Z_OK) {
            std::cerr << "Error: deflateInit failed.\n";
            return false;
        }

        while (in) {
            in.read(buffer, bufferSize);
            std::streamsize bytesRead = in.gcount();
            if (bytesRead > 0) {
                strm.avail_in = static_cast<uInt>(bytesRead);
                strm.next_in = reinterpret_cast<Bytef*>(buffer);

                do {
                    strm.avail_out = bufferSize;
                    strm.next_out = reinterpret_cast<Bytef*>(buffer);
                    deflate(&strm, in.eof() ? Z_FINISH : Z_NO_FLUSH);
                    std::streamsize have = bufferSize - strm.avail_out;
                    out.write(buffer, have);
                } while (strm.avail_out == 0);
            }
        }

        deflateEnd(&strm);
    } else {
        // Initialize zlib for decompression
        z_stream strm;
        std::memset(&strm, 0, sizeof(strm));
        if (inflateInit(&strm) != Z_OK) {
            std::cerr << "Error: inflateInit failed.\n";
            return false;
        }

        while (in) {
            in.read(buffer, bufferSize);
            std::streamsize bytesRead = in.gcount();
            if (bytesRead > 0) {
                strm.avail_in = static_cast<uInt>(bytesRead);
                strm.next_in = reinterpret_cast<Bytef*>(buffer);

                do {
                    strm.avail_out = bufferSize;
                    strm.next_out = reinterpret_cast<Bytef*>(buffer);
                    int ret = inflate(&strm, in.eof() ? Z_FINISH : Z_NO_FLUSH);
                    if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR) {
                        std::cerr << "Error: inflate failed with code " << ret << ".\n";
                        inflateEnd(&strm);
                        return false;
                    }
                    std::streamsize have = bufferSize - strm.avail_out;
                    out.write(buffer, have);
                } while (strm.avail_out == 0);
            }
        }

        inflateEnd(&strm);
    }

    return true;
}

int main(int argc, char* argv[]) {
    bool compress = false;
    bool decompress = false;

    // Argument parsing for help and operation mode
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--compress" || arg == "-c") {
            compress = true;
        } else if (arg == "--decompress" || arg == "-d") {
            decompress = true;
        } else if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    if ((compress && decompress) || (!compress && !decompress)) {
        std::cerr << "Error: Specify either --compress or --decompress.\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    if (!compress && !decompress) {
        printHelp();
        return 1;
    }

    bool success = compressDecompressVCF(std::cin, std::cout, compress);
    return success ? 0 : 1;
}
