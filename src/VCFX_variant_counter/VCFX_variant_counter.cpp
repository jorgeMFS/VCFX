#include "VCFX_variant_counter.h"
#include <getopt.h>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include "vcfx_core.h"

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
    
    std::string plainInput;
    if(!vcfx::read_maybe_compressed(std::cin, plainInput)){
        std::cerr << "Error: failed to read input" << std::endl;
        return 1;
    }
    std::istringstream inStream(plainInput);
    int total= countVariants(inStream);
    if(total<0){
        // indicates an error if strict
        return 1;
    }
    std::cout<<"Total Variants: "<< total <<"\n";
    return 0;
}

int VCFXVariantCounter::countVariants(std::istream &in){
    int count=0;
    int lineNumber=0;
    std::string line;
    while(true){
        if(!std::getline(in,line)) break;
        lineNumber++;
        if(line.empty()) continue;
        if(line[0]=='#') continue; 
        // parse columns
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
                return -1; // indicates error
            } else {
                std::cerr<<"Warning: skipping line "<<lineNumber<<" with <8 columns.\n";
                continue;
            }
        }
        // if we get here => count it
        count++;
    }
    return count;
}

int main(int argc, char* argv[]) {
    if (vcfx::handle_version_flag(argc, argv, "VCFX_variant_counter")) return 0;
    VCFXVariantCounter app;
    return app.run(argc, argv);
}
