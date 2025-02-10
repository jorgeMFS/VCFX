#include "VCFX_phase_quality_filter.h"
#include <getopt.h>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cctype>

static void split(const std::string &s, char delim, std::vector<std::string> &tokens) {
    tokens.clear();
    std::stringstream ss(s);
    std::string t;
    while(std::getline(ss,t,delim)) {
        tokens.push_back(t);
    }
}

int VCFXPhaseQualityFilter::run(int argc, char* argv[]) {
    if(argc==1) {
        displayHelp();
        return 0;
    }
    bool showHelp=false;
    std::string condition;
    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"filter-pq", required_argument, 0, 'f'},
        {0,0,0,0}
    };
    while(true) {
        int c= ::getopt_long(argc, argv, "hf:", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 'f':
                condition= optarg;
                break;
            default:
                showHelp=true;
        }
    }
    if(showHelp) {
        displayHelp();
        return 0;
    }
    if(condition.empty()) {
        std::cerr<<"Error: Must specify condition with --filter-pq\n";
        displayHelp();
        return 1;
    }
    std::string op;
    double threshold=0.0;
    if(!parseCondition(condition, op, threshold)) {
        std::cerr<<"Error: Unable to parse condition '"<< condition<<"'. e.g. PQ>=30\n";
        return 1;
    }
    filterByPQ(std::cin, std::cout, op, threshold);
    return 0;
}

void VCFXPhaseQualityFilter::displayHelp() {
    std::cout <<
"VCFX_phase_quality_filter: Filter variants by phasing quality (PQ) in the INFO field.\n\n"
"Usage:\n"
"  VCFX_phase_quality_filter --filter-pq \"PQ<OP><THRESHOLD>\" < input.vcf > output.vcf\n\n"
"Options:\n"
"  -h, --help             Print this help message.\n"
"  -f, --filter-pq <COND> Condition like 'PQ>30', 'PQ>=20', 'PQ!=10', etc.\n\n"
"Description:\n"
"  Reads each variant line, extracts 'PQ=' from INFO. If missing or invalid, PQ=0.\n"
"  Keeps lines if 'PQ <OP> THRESHOLD' is true. Otherwise, discards.\n\n"
"Supported operators: >, >=, <, <=, ==, !=\n\n"
"Examples:\n"
"  1) Keep variants with PQ>30:\n"
"     VCFX_phase_quality_filter --filter-pq \"PQ>30\" < in.vcf > out.vcf\n"
"  2) Keep PQ<=15:\n"
"     VCFX_phase_quality_filter --filter-pq \"PQ<=15\" < in.vcf > out.vcf\n";
}

bool VCFXPhaseQualityFilter::parseCondition(const std::string &condition,
                                            std::string &op,
                                            double &threshold) {
    if(condition.size()<3) return false;
    if(condition.rfind("PQ",0)!=0) {
        return false;
    }
    auto sub= condition.substr(2); 
    std::vector<std::string> ops={">=","<=","==","!=",">","<"};
    std::string matched;
    size_t foundPos=std::string::npos;
    for(auto &o : ops){
        auto p= sub.find(o);
        if(p==0) { matched=o; break; }
    }
    if(matched.empty()) {
        return false;
    }
    op= matched;
    auto valStr= sub.substr(op.size());
    if(valStr.empty()) return false;
    try {
        threshold= std::stod(valStr);
    } catch(...){
        return false;
    }
    return true;
}

void VCFXPhaseQualityFilter::filterByPQ(std::istream &in, std::ostream &out,
                                        const std::string &op,
                                        double threshold) {
    bool headerFound=false;
    std::string line;
    while(true) {
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
            std::cerr<<"Warning: data line before #CHROM => skipping\n";
            continue;
        }
        std::vector<std::string> fields;
        {
            std::stringstream ss(line);
            std::string f;
            while(std::getline(ss,f,'\t')) fields.push_back(f);
        }
        if(fields.size()<8){
            std::cerr<<"Warning: VCF line with fewer than 8 columns => skipping.\n";
            continue;
        }
        // parse info
        double pq= parsePQScore(fields[7]);
        bool keep=false;
        if(op==">") {
            if(pq> threshold) keep=true;
        } else if(op==">=") {
            if(pq>=threshold) keep=true;
        } else if(op=="<") {
            if(pq< threshold) keep=true;
        } else if(op=="<=") {
            if(pq<=threshold) keep=true;
        } else if(op=="==") {
            if(pq==threshold) keep=true;
        } else if(op=="!=") {
            if(pq!=threshold) keep=true;
        }
        if(keep) {
            out<< line <<"\n";
        }
    }
}

double VCFXPhaseQualityFilter::parsePQScore(const std::string &info) {
    if(info.empty() || info=="." ) return 0.0;
    std::stringstream ss(info);
    std::string kv;
    while(std::getline(ss, kv, ';')) {
        if(kv.rfind("PQ=",0)==0) {
            auto valStr= kv.substr(3);
            try {
                return std::stod(valStr);
            } catch(...) {
                std::cerr<<"Warning: invalid PQ= '"<<valStr<<"'. Using 0.0.\n";
                return 0.0;
            }
        }
    }
    return 0.0;
}

int main(int argc, char* argv[]){
    VCFXPhaseQualityFilter f;
    return f.run(argc, argv);
}
