#include "VCFX_indexer.h"
#include <getopt.h>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstdio>

// A small utility: splits a string by delimiter
static std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string t;
    while(std::getline(ss,t,delim)) {
        tokens.push_back(t);
    }
    return tokens;
}

void VCFXIndexer::displayHelp() {
    std::cout 
        << "VCFX_indexer\n"
        << "Usage: VCFX_indexer [--help]\n\n"
        << "Description:\n"
        << "  Reads a VCF from stdin and writes an index to stdout. The index is a TSV with\n"
        << "  columns: CHROM, POS, FILE_OFFSET. The FILE_OFFSET is the byte offset from the\n"
        << "  start of the file to the beginning of that variant line.\n\n"
        << "Example:\n"
        << "  VCFX_indexer < input.vcf > index.tsv\n";
}

int VCFXIndexer::run(int argc, char* argv[]) {
    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };

    while(true) {
        int c= getopt_long(argc, argv, "h", long_opts, nullptr);
        if (c==-1) break;
        switch(c){
            case 'h':
                displayHelp();
                return 0;
            default:
                displayHelp();
                return 1;
        }
    }

    // Now we do the main
    createVCFIndex(std::cin, std::cout);
    return 0;
}

// The main function to create the index
void VCFXIndexer::createVCFIndex(std::istream &in, std::ostream &out) {
    // We'll read from 'in' using a method that tracks the byte offset
    // so we can't just do std::getline alone. We can do a manual read approach or we use 'tellg()' before reading each line
    out << "CHROM\tPOS\tFILE_OFFSET\n";

    bool foundChromHeader=false;

    // We track the offset to each line by calling tellg() at the start of reading each line
    // caution: if the input is not seekable (like a pipe), tellg() might return -1
    // We'll fallback to a manual offset if tellg() doesn't work
    bool canSeek = true;
    std::streampos lastPos = in.tellg();
    if(lastPos < 0) {
        // not seekable => we do a manual counting approach
        canSeek=false;
    }

    // We'll do a manual offset count
    // for each line we read, we add line.size()+1 for the newline
    // but for cross-platform newlines we might mismatch. We'll do a simpler approach
    // If canSeek is false, we do manual counting
    // If canSeek is true, we rely on tellg()
    long manualOffset= 0;

    while(true) {
        std::streampos fpos= in.tellg(); // get current position
        if(fpos<0 && !canSeek) {
            // fallback => use manualOffset
        } else if(fpos<0 && canSeek) {
            // error
        }
        if(in.eof()) break;
        if(!in.good()) break;

        long lineOffset = canSeek ? (long)fpos : manualOffset;

        std::string line;
        if(!std::getline(in,line)) {
            break; 
        }
        // If the line is empty, continue
        if(line.empty()) {
            // add length
            if(!canSeek) {
                manualOffset += (long)line.size()+1; 
            }
            continue;
        }
        // If it's a header line
        if(line[0]=='#') {
            // check if it's #CHROM => set foundChromHeader
            if(!foundChromHeader && line.rfind("#CHROM",0)==0) {
                foundChromHeader= true;
            }
            // either way no index line
            if(!canSeek) {
                manualOffset += (long)line.size()+1;
            }
            continue;
        }

        // If we have not found #CHROM, that's an error
        if(!foundChromHeader) {
            std::cerr<< "Error: no #CHROM header found before variant lines.\n";
            return;
        }

        // parse
        auto fields = split(line,'\t');
        if(fields.size()<2) {
            // skip
            if(!canSeek) {
                manualOffset += (long)line.size()+1;
            }
            continue;
        }
        std::string &chrom = fields[0];
        std::string &posStr= fields[1];
        int posVal=0;
        try {
            posVal = std::stoi(posStr);
        } catch(...) {
            // skip
            if(!canSeek) {
                manualOffset += (long)line.size()+1;
            }
            continue;
        }
        // output index line
        out << chrom << "\t" << posVal << "\t" << lineOffset << "\n";

        // update offset
        if(!canSeek) {
            manualOffset += (long)line.size()+1; 
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXIndexer idx;
    return idx.run(argc, argv);
}
