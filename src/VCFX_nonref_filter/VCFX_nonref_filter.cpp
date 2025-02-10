#include "VCFX_nonref_filter.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <cctype>
#include <cstdlib>

int VCFXNonRefFilter::run(int argc, char* argv[]) {
    int opt; bool showHelp=false;
    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };
    while(true){
        int c= getopt_long(argc, argv, "h", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
            default:
                showHelp=true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }
    filterNonRef(std::cin, std::cout);
    return 0;
}

void VCFXNonRefFilter::displayHelp(){
    std::cout <<
"VCFX_nonref_filter: Filter out variants where all samples are homozygous reference.\n\n"
"Usage:\n"
"  VCFX_nonref_filter [options] < input.vcf > output.vcf\n\n"
"Description:\n"
"  Reads a VCF and checks each variant. If in every sample genotype\n"
"  all alleles are '0' (homozygous reference), that variant is skipped.\n"
"  Otherwise, the variant is printed.\n"
"  Missing genotypes or partial calls are considered non-hom-ref.\n\n"
"Options:\n"
"  --help, -h   Print this help.\n\n"
"Example:\n"
"  VCFX_nonref_filter < input.vcf > nohomref.vcf\n";
}

bool VCFXNonRefFilter::isHomRef(const std::string &genotypeString) const {
    if(genotypeString.empty()|| genotypeString=="."|| genotypeString=="./."|| genotypeString==".|.") return false;
    std::string g(genotypeString);
    for(auto &c:g) if(c=='|') c='/';
    // split by '/'
    auto slashPos= g.find('/');
    if(slashPos== std::string::npos){
        // Could be haploid or partial
        // if "0" => maybe homRef haploid. If "." => missing => not hom-ref
        // If "1" => nonRef
        // We'll do: if g=="0" => homRef, else false
        return (g=="0");
    }
    // if multiple slashes => e.g. polyploid => parse them all
    // let's do a quick approach: split by '/'
    std::vector<std::string> alleles;
    {
        std::stringstream ss(g);
        std::string token;
        while(std::getline(ss, token, '/')){
            alleles.push_back(token);
        }
    }
    // if any allele is not "0" => not hom-ref
    for(auto &a: alleles){
        if(a!="0") return false;
    }
    return true;
}

void VCFXNonRefFilter::filterNonRef(std::istream& in, std::ostream& out) {
    bool headerFound=false;
    std::string line;
    while(true){
        if(!std::getline(in,line)) break;
        if(line.empty()){
            out<< line <<"\n";
            continue;
        }
        if(line[0]=='#'){
            out<< line <<"\n";
            if(line.rfind("#CHROM",0)==0) headerFound=true;
            continue;
        }
        if(!headerFound){
            std::cerr<<"Error: data line encountered before #CHROM header.\n";
            // we can break or skip
            // We'll just continue
            continue;
        }
        std::stringstream ss(line);
        std::vector<std::string> fields; {
            std::string t;
            while(std::getline(ss,t,'\t')) fields.push_back(t);
        }
        if(fields.size()<10){
            // fewer than CHROM..INFO,FORMAT + 1 sample => pass or skip?
            out<< line <<"\n";
            continue;
        }
        std::string &format= fields[8];
        auto fmtParts= std::vector<std::string>();
        {
            std::stringstream fs(format);
            std::string ff; 
            while(std::getline(fs, ff, ':')){
                fmtParts.push_back(ff);
            }
        }
        int gtIndex=-1;
        for(size_t i=0;i<fmtParts.size();i++){
            if(fmtParts[i]=="GT"){ gtIndex=(int)i; break;}
        }
        if(gtIndex<0){
            // No GT => we can't confirm homRef => we keep
            out<< line <<"\n";
            continue;
        }
        bool allHomRef=true;
        for(size_t s=9; s< fields.size(); s++){
            // parse genotype
            auto sampleCol= fields[s];
            std::vector<std::string> subf;
            {
                std::stringstream st(sampleCol);
                std::string token; 
                while(std::getline(st, token, ':')){
                    subf.push_back(token);
                }
            }
            if(gtIndex>=(int)subf.size()){
                // missing => not hom ref
                allHomRef=false; 
                break;
            }
            if(!isHomRef(subf[gtIndex])){
                allHomRef=false;
                break;
            }
        }
        if(!allHomRef) out<< line <<"\n";
    }
}
