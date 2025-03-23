#include "VCFX_sorter.h"
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cstdlib>

void VCFXSorter::displayHelp(){
    std::cout <<
"VCFX_sorter: Sort a VCF by chromosome and position.\n\n"
"Usage:\n"
"  VCFX_sorter [options] < input.vcf > output.vcf\n\n"
"Options:\n"
"  -h, --help          Show help.\n"
"  -n, --natural-chr   Use a natural chromosome sort (chr1 < chr2 < chr10) instead of lexicographic.\n\n"
"Description:\n"
"  Reads all data lines into memory, sorts by (CHROM,POS). Preserves all header lines\n"
"  in original order, and outputs them first, then prints sorted data lines.\n\n"
"Examples:\n"
"  1) Lexicographic:\n"
"     VCFX_sorter < unsorted.vcf > sorted.vcf\n"
"  2) Natural order:\n"
"     VCFX_sorter --natural-chr < unsorted.vcf > sorted.vcf\n";
}

// A function to parse chromosome name in a natural manner if possible
// e.g. "chr10" => (chr,10) so that chr2 < chr10 in numeric sense
// We'll remove "chr" prefix if present, then parse leading digits
static bool parseChromNat(const std::string &chrom,
                          std::string &prefix,
                          long &num,
                          std::string &suffix)
{
    // e.g. "chr10_gl000" => prefix="chr", then 10, then suffix="_gl000"
    // We'll do a naive approach: remove leading "chr" or "Chr" or "CHR"
    std::string c= chrom;
    // uppercase to check
    std::string up;
    up.reserve(c.size());
    for(char ch: c) up.push_back(std::toupper(ch));
    if(up.rfind("CHR",0)==0) {
        // remove "chr" prefix
        prefix= c.substr(0,3); 
        c= c.substr(3);
    } else {
        prefix="";
    }
    // now parse leading digits
    // e.g. c="10_gl000"
    // we parse digits until non-digit
    size_t idx=0;
    while(idx< c.size() && std::isdigit((unsigned char)c[idx])) idx++;
    if(idx==0){
        // means no leading digit
        num=-1;
        suffix= c;
        return true;
    }
    // parse c[0..idx-1] as a number
    try {
        num= std::stol(c.substr(0,idx));
    } catch(...){
        return false;
    }
    suffix= c.substr(idx);
    return true;
}

// We'll define the lexCompare
bool VCFRecord::lexCompare(const VCFRecord &a, const VCFRecord &b){
    if(a.chrom != b.chrom) return a.chrom < b.chrom;
    return a.pos < b.pos;
}

// We'll define the naturalCompare
bool VCFRecord::naturalCompare(const VCFRecord &a, const VCFRecord &b){
    // parse each
    std::string apfx, asuf, bpfx, bsuf;
    long anum, bnum;
    if(!parseChromNat(a.chrom, apfx, anum, asuf) || 
       !parseChromNat(b.chrom, bpfx, bnum, bsuf)) {
        // fallback to lex
        if(a.chrom != b.chrom) return a.chrom < b.chrom;
        return a.pos< b.pos;
    }
    // if prefix differs, lex compare them
    if(apfx != bpfx) return (apfx< bpfx);
    // if both have num>=0
    if(anum>=0 && bnum>=0){
        if(anum!= bnum) return (anum< bnum);
        // if suffix differs => lex compare
        if(asuf!= bsuf) return asuf< bsuf;
        // else compare pos
        if(a.pos!= b.pos) return a.pos< b.pos;
        // all tie
        return false;
    } else if(anum>=0 && bnum<0){
        // numeric < no numeric
        return true;
    } else if(anum<0 && bnum>=0){
        return false;
    } else {
        // both have no numeric => fallback to lex compare
        if(a.chrom!= b.chrom) return a.chrom< b.chrom;
        return a.pos< b.pos;
    }
}

int VCFXSorter::run(int argc, char* argv[]){
    bool showHelp=false;
    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"natural-chr", no_argument, 0, 'n'},
        {0,0,0,0}
    };
    while(true){
        int c= ::getopt_long(argc, argv, "hn", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 'n':
                naturalChromOrder= true;
                break;
            default:
                showHelp=true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }
    
    // Check if stdin has data
    if (std::cin.peek() == EOF) {
        displayHelp();
        return 0;
    }
    
    loadVCF(std::cin);
    sortRecords();
    outputVCF(std::cout);
    return 0;
}

// loads lines into memory
void VCFXSorter::loadVCF(std::istream &in){
    bool foundChrom=false;
    std::string line;
    while(true){
        if(!std::getline(in,line)) break;
        if(line.empty()){
            // we keep it in the header or skip? We'll keep it in header
            headerLines.push_back(line);
            continue;
        }
        if(line[0]=='#'){
            headerLines.push_back(line);
            if(line.rfind("#CHROM",0)==0) foundChrom=true;
            continue;
        }
        // parse data
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string col;
            while(std::getline(ss,col,'\t')){
                fields.push_back(col);
            }
        }
        if(fields.size()<8){
            std::cerr<<"Warning: skipping line with <8 columns.\n";
            continue;
        }
        VCFRecord rec;
        rec.chrom= fields[0];
        try {
            rec.pos= std::stoi(fields[1]);
        } catch(...){
            std::cerr<<"Warning: invalid POS => skipping line.\n";
            continue;
        }
        rec.fields= std::move(fields);
        records.push_back(rec);
    }
    if(!foundChrom){
        std::cerr<<"Warning: no #CHROM line found in input.\n";
    }
}

// sort them in memory
void VCFXSorter::sortRecords(){
    if(naturalChromOrder){
        std::sort(records.begin(), records.end(), VCFRecord::naturalCompare);
    } else {
        std::sort(records.begin(), records.end(), VCFRecord::lexCompare);
    }
}

void VCFXSorter::outputVCF(std::ostream &out){
    // print all header lines
    for(auto &l : headerLines){
        out<< l <<"\n";
    }
    // then print sorted data
    for(auto &rec: records){
        // rebuild line from rec.fields
        bool first=true;
        for(size_t i=0;i< rec.fields.size(); i++){
            if(!first) out<<"\t";
            else first=false;
            out<< rec.fields[i];
        }
        out<<"\n";
    }
}

int main(int argc, char* argv[]){
    VCFXSorter app;
    return app.run(argc, argv);
}