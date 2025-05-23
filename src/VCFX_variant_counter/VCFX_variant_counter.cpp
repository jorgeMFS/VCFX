#include "VCFX_variant_counter.h"
#include <getopt.h>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include "vcfx_core.h"
#include <zlib.h>
#include <cstring>

void VCFXVariantCounter::displayHelp(){
    std::cout <<
"VCFX_variant_counter: Counts the total number of valid variants in a VCF.\n\n"
"Usage:\n"
"  VCFX_variant_counter [options] < input.vcf\n\n"
"Options:\n"
"  -h, --help        Show this help.\n"
"  -s, --strict      Fail on any data line with <8 columns.\n\n"
"Description:\n"
"  Reads a VCF from stdin, ignores all header lines (#). For each data line,\n"
"  we check if it has >=8 columns; if it does, we count it; if fewer columns:\n"
"   * if --strict => we exit with error,\n"
"   * otherwise => we skip with a warning.\n"
"  Finally, we print 'Total Variants: X'.\n\n"
"Example:\n"
"  VCFX_variant_counter < input.vcf\n"
"  VCFX_variant_counter --strict < input.vcf\n";
}

int VCFXVariantCounter::run(int argc, char* argv[]){
    bool showHelp=false;
    static struct option long_opts[]={
        {"help", no_argument, 0,'h'},
        {"strict", no_argument, 0,'s'},
        {0,0,0,0}
    };
    
    // Check for options/flags
    if (argc > 1) {
        while(true){
            int c= ::getopt_long(argc, argv,"hs", long_opts,nullptr);
            if(c==-1) break;
            switch(c){
                case 'h':
                    showHelp=true;
                    break;
                case 's':
                    strictMode= true;
                    break;
                default:
                    showHelp=true;
            }
        }
    }
    
    if(showHelp){
        displayHelp();
        return 0;
    }
    
    auto peek1 = std::cin.peek();
    bool isEmpty = (peek1 == EOF);
    bool isGzip = false;
    if(!isEmpty){
        int c1 = std::cin.get();
        int c2 = std::cin.get();
        if(c2 != EOF){
            isGzip = (static_cast<unsigned char>(c1) == 0x1f &&
                      static_cast<unsigned char>(c2) == 0x8b);
            std::cin.putback(static_cast<char>(c2));
        }
        std::cin.putback(static_cast<char>(c1));
    }

    int total = -1;
    if(isEmpty){
        total = 0;
    } else if(isGzip){
        total = countVariantsGzip(std::cin);
    } else {
        total = countVariants(std::cin);
    }
    if(total<0){
        // indicates an error if strict
        return 1;
    }
    std::cout<<"Total Variants: "<< total <<"\n";
    return 0;
}

bool VCFXVariantCounter::processLine(const std::string &line, int lineNumber, int &count){
    if(line.empty()) return true;
    if(line[0]=='#') return true;
    std::stringstream ss(line);
    std::vector<std::string> fields;
    {
        std::string col;
        while(std::getline(ss,col,'\t')){
            fields.push_back(col);
        }
    }
    if(fields.size()<8){
        if(strictMode){
            std::cerr<<"Error: line "<< lineNumber <<" has <8 columns.\n";
            return false;
        } else {
            std::cerr<<"Warning: skipping line "<<lineNumber<<" with <8 columns.\n";
            return true;
        }
    }
    count++;
    return true;
}

int VCFXVariantCounter::countVariants(std::istream &in){
    int count=0;
    int lineNumber=0;
    std::string line;
    while(std::getline(in,line)){
        lineNumber++;
        if(!processLine(line, lineNumber, count)) return -1;
    }
    return count;
}

int VCFXVariantCounter::countVariantsGzip(std::istream &in){
    constexpr int CHUNK = 16384;
    char inBuf[CHUNK];
    char outBuf[CHUNK];
    z_stream strm; std::memset(&strm,0,sizeof(strm));
    if(inflateInit2(&strm,15+32)!=Z_OK){
        std::cerr<<"Error: inflateInit2 failed.\n";
        return -1;
    }
    int count=0; int lineNumber=0; std::string buffer; int ret=Z_OK;
    do {
        in.read(inBuf, CHUNK);
        strm.avail_in = static_cast<uInt>(in.gcount());
        if(strm.avail_in==0 && in.eof()) break;
        strm.next_in = reinterpret_cast<Bytef*>(inBuf);
        do {
            strm.avail_out = CHUNK;
            strm.next_out = reinterpret_cast<Bytef*>(outBuf);
            ret = inflate(&strm, Z_NO_FLUSH);
            if (ret == Z_STREAM_ERROR || ret == Z_NEED_DICT || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR){
                std::cerr<<"Error: decompression failed.\n";
                inflateEnd(&strm);
                return -1;
            }
            size_t have = CHUNK - strm.avail_out;
            if(have>0){
                buffer.append(outBuf, have);
                size_t pos;
                while((pos = buffer.find('\n')) != std::string::npos){
                    std::string line = buffer.substr(0,pos);
                    buffer.erase(0,pos+1);
                    lineNumber++;
                    if(!processLine(line,lineNumber,count)){
                        inflateEnd(&strm);
                        return -1;
                    }
                }
            }
        } while(strm.avail_out==0);
    } while(ret != Z_STREAM_END);

    if(!buffer.empty()){
        lineNumber++;
        if(!processLine(buffer,lineNumber,count)){
            inflateEnd(&strm);
            return -1;
        }
    }
    inflateEnd(&strm);
    return count;
}

int main(int argc, char* argv[]) {
    if (vcfx::handle_version_flag(argc, argv, "VCFX_variant_counter")) return 0;
    VCFXVariantCounter app;
    return app.run(argc, argv);
}
