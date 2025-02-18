#include "VCFX_subsampler.h"
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <cstdlib>

void VCFXSubsampler::displayHelp(){
    std::cout <<
"VCFX_subsampler: Randomly pick N lines from a VCF data section.\n\n"
"Usage:\n"
"  VCFX_subsampler [options] < input.vcf > output.vcf\n\n"
"Options:\n"
"  -s, --subsample <N>   Required: number of data lines (variants) to keep.\n"
"  --seed <INT>          Use a reproducible random seed.\n"
"  -h, --help            Show this help.\n\n"
"Description:\n"
"  We read all header lines (#...) first and output them as-is. Then we do\n"
"  reservoir sampling on subsequent lines (the data lines). If the file has\n"
"  fewer than N lines, we keep them all. We skip lines with <8 columns.\n\n"
"Example:\n"
"  VCFX_subsampler --subsample 1000 < big.vcf > subset.vcf\n"
"  VCFX_subsampler --subsample 1000 --seed 1234 < big.vcf > subset2.vcf\n";
}

int VCFXSubsampler::run(int argc, char* argv[]){
    if(argc==1){
        displayHelp();
        return 0;
    }
    bool showHelp=false;
    int sampleSize=0;
    unsigned int seed= (unsigned int)std::time(nullptr);
    bool seedSpecified= false;

    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"subsample", required_argument, 0, 's'},
        {"seed", required_argument, 0, 1000},
        {0,0,0,0}
    };
    while(true){
        int c= ::getopt_long(argc, argv, "hs:", long_opts,nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 's': {
                try {
                    sampleSize= std::stoi(optarg);
                    if(sampleSize<=0){
                        throw std::invalid_argument("must be >0");
                    }
                } catch(...){
                    std::cerr<<"Error: invalid subsample size.\n";
                    return 1;
                }
            } break;
            case 1000: { // --seed
                seedSpecified= true;
                try {
                    long long val= std::stoll(optarg);
                    seed= (unsigned int)val;
                } catch(...){
                    std::cerr<<"Error: invalid seed.\n";
                    return 1;
                }
            }break;
            default:
                showHelp=true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }
    if(sampleSize<=0){
        std::cerr<<"Error: must specify --subsample <N> with N>0.\n";
        return 1;
    }
    subsampleLines(std::cin, std::cout, sampleSize, seed);
    return 0;
}

void VCFXSubsampler::subsampleLines(std::istream &in, std::ostream &out,
                                    int sampleSize, unsigned int seed)
{
    std::string line;
    // store all # lines first
    while(true){
        std::streampos p= in.tellg();
        if(!std::getline(in,line)) break;
        if(line.empty()){
            out<< line<<"\n";
            continue;
        }
        if(line[0]=='#'){
            out<< line<<"\n";
            continue;
        } else {
            // not a header
            in.seekg(p);
            break;
        }
    }
    // now do reservoir sampling on data lines
    std::vector<std::string> reservoir;
    reservoir.reserve(sampleSize);

    int count=0;
    // random gen
    std::default_random_engine rng(seed);
    // read data lines
    while(true){
        if(!std::getline(in, line)) break;
        if(line.empty()) continue;
        // skip lines with <8 columns
        {
            std::stringstream s2(line);
            std::vector<std::string> fields;
            std::string tmp;
            while(std::getline(s2,tmp,'\t')){
                fields.push_back(tmp);
            }
            if(fields.size()<8){
                std::cerr<<"Warning: skipping line with <8 columns.\n";
                continue;
            }
        }

        if(count< sampleSize){
            reservoir.push_back(line);
        } else {
            // pick random int in [0..count]
            std::uniform_int_distribution<int> dist(0, count);
            int j= dist(rng);
            if(j< sampleSize){
                reservoir[j]= line;
            }
        }
        count++;
    }
    // output reservoir
    for(auto &r: reservoir){
        out<< r<<"\n";
    }
}


int main(int argc, char* argv[]){
    VCFXSubsampler app;
    return app.run(argc, argv);
}
