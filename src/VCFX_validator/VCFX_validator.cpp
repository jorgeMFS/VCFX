#include "VCFX_validator.h"
#include <getopt.h>
#include <sstream>
#include <vector>
#include <cctype>
#include <cstdlib>

static std::string trim(const std::string &s){
    size_t start=0; 
    while(start<s.size() && std::isspace((unsigned char)s[start])) start++;
    if(start== s.size()) return "";
    size_t end= s.size()-1;
    while(end> start && std::isspace((unsigned char)s[end])) end--;
    return s.substr(start, end-start+1);
}

int VCFXValidator::run(int argc, char* argv[]){
    if(argc==1){
        displayHelp();
        return 0;
    }
    bool showHelp=false;
    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"strict", no_argument, 0, 's'},
        {0,0,0,0}
    };
    while(true){
        int c= ::getopt_long(argc, argv, "hs", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h': showHelp=true; break;
            case 's': strictMode= true; break;
            default: showHelp= true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }
    bool ok= validateVCF(std::cin);
    return (ok ? 0 : 1);
}

void VCFXValidator::displayHelp(){
    std::cout <<
"VCFX_validator: Checks basic validity of a VCF.\n\n"
"Usage:\n"
"  VCFX_validator [options] < input.vcf\n\n"
"Options:\n"
"  -h, --help     Show this help.\n"
"  -s, --strict   Enable stricter checks (not fully implemented, but reserved).\n\n"
"Description:\n"
"  Validates:\n"
"   * All '##' lines are recognized as meta lines.\n"
"   * #CHROM line is present, has at least 8 columns.\n"
"   * Each data line has >=8 columns, checks CHROM non-empty, POS>0,\n"
"     REF/ALT non-empty, QUAL is '.' or non-negative float, FILTER non-empty,\n"
"     INFO is minimal check. Logs errors/warnings.\n"
"  Exits 0 if pass, 1 if fail.\n";
}

bool VCFXValidator::validateMetaLine(const std::string &line, int lineNumber){
    // minimal check: must start with "##"
    if(line.size()<2) return false;
    if(line[0]=='#' && line[1]=='#'){
        return true;
    }
    std::cerr<<"Error: line "<< lineNumber <<" is a header line but doesn't start with '##'.\n";
    return false;
}

bool VCFXValidator::validateChromHeader(const std::string &line, int lineNumber){
    // parse by tab
    std::vector<std::string> f;
    {
        std::stringstream ss(line);
        std::string col;
        while(std::getline(ss,col,'\t')){
            f.push_back(col);
        }
    }
    if(f.size()<8){
        std::cerr<<"Error: #CHROM line at "<< lineNumber <<" has <8 columns.\n";
        return false;
    }
    // typically #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, (FORMAT?), ...
    if(f[0]!="#CHROM"){
        std::cerr<<"Error: #CHROM line doesn't start with '#CHROM' at line "<< lineNumber <<".\n";
        return false;
    }
    return true;
}

bool VCFXValidator::validateDataLine(const std::string &line, int lineNumber){
    // parse
    std::vector<std::string> f;
    {
        std::stringstream ss(line);
        std::string col;
        while(std::getline(ss,col,'\t')){
            f.push_back(trim(col));
        }
    }
    if(f.size()<8){
        std::cerr<<"Error: line "<< lineNumber <<" has <8 columns.\n";
        return false;
    }
    // CHROM
    if(f[0].empty()){
        std::cerr<<"Error: line "<<lineNumber<<" CHROM is empty.\n";
        return false;
    }
    // POS
    {
        int pos=0;
        try{
            pos= std::stoi(f[1]);
            if(pos<=0){
                std::cerr<<"Error: line "<<lineNumber<<" POS must be >0.\n";
                return false;
            }
        } catch(...){
            std::cerr<<"Error: line "<<lineNumber<<" POS not parseable.\n";
            return false;
        }
    }
    // ID can be empty => that's fine
    // REF
    if(f[3].empty()){
        std::cerr<<"Error: line "<< lineNumber<<" REF is empty.\n";
        return false;
    }
    // ALT
    if(f[4].empty()){
        std::cerr<<"Error: line "<< lineNumber<<" ALT is empty.\n";
        return false;
    }
    // QUAL => '.' or non-neg float
    if(f[5]!="."){
        try{
            double q= std::stod(f[5]);
            if(q<0){
                std::cerr<<"Error: line "<<lineNumber<<" negative QUAL.\n";
                return false;
            }
        } catch(...){
            std::cerr<<"Error: line "<<lineNumber<<" invalid QUAL.\n";
            return false;
        }
    }
    // FILTER => must not be empty
    if(f[6].empty()){
        std::cerr<<"Error: line "<<lineNumber<<" FILTER is empty.\n";
        return false;
    }
    // INFO => minimal check
    // could be '.' => ok
    // else splitted by ';' => each either 'key=val' or 'flag'
    if(f[7]!="."){
        std::stringstream infoSS(f[7]);
        std::string token;
        bool anyValid= false;
        while(std::getline(infoSS, token, ';')){
            token= trim(token);
            if(token.empty()) continue;
            auto eq= token.find('=');
            if(eq== std::string::npos){
                // treat as flag
                anyValid= true;
            } else {
                // key=val
                std::string k= token.substr(0, eq);
                // std::string v= token.substr(eq+1);
                if(k.empty()){
                    std::cerr<<"Error: line "<<lineNumber<<" has INFO with empty key.\n";
                    return false;
                }
                anyValid= true;
            }
        }
        // if not anyValid => error
        if(!anyValid){
            std::cerr<<"Error: line "<<lineNumber<<" INFO not valid.\n";
            return false;
        }
    }
    return true;
}

bool VCFXValidator::validateVCF(std::istream &in){
    std::string line;
    int lineNum=0;
    bool foundChromLine= false;

    while(true){
        if(!std::getline(in, line)) break;
        lineNum++;
        if(line.empty()){
            // skip
            continue;
        }
        if(line[0]=='#'){
            // could be ## or #CHROM
            if(line.rfind("##",0)==0){
                // validate meta
                if(!validateMetaLine(line, lineNum)) return false;
            } else if(line.rfind("#CHROM",0)==0){
                // validate #CHROM
                if(!validateChromHeader(line, lineNum)) return false;
                foundChromLine= true;
            } else {
                // might be some weird # line? We'll allow # as well but let's do a basic check
                // We'll treat lines that start with '#' but not "##" or "#CHROM" as okay
            }
        } else {
            // data line
            if(!foundChromLine){
                std::cerr<<"Error: data line encountered before #CHROM at line "<<lineNum<<".\n";
                return false;
            }
            if(!validateDataLine(line, lineNum)) return false;
        }
    }
    if(!foundChromLine){
        std::cerr<<"Error: no #CHROM line found in file.\n";
        return false;
    }
    std::cout<<"VCF file is valid.\n";
    return true;
}

int main(int argc, char* argv[]){
    VCFXValidator validator;
    return validator.run(argc, argv);
}
