#include "vcfx_core.h"
#include "VCFX_validator.h"
#include <getopt.h>
#include <sstream>
#include <vector>
#include <cctype>
#include <cstdlib>
#include <unistd.h>

static std::string trim(const std::string &s){
    size_t start=0;
    while(start<s.size() && std::isspace((unsigned char)s[start])) start++;
    if(start== s.size()) return "";
    size_t end= s.size()-1;
    while(end> start && std::isspace((unsigned char)s[end])) end--;
    return s.substr(start, end-start+1);
}

static std::vector<std::string> split(const std::string &s, char delim){
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) out.push_back(item);
    return out;
}

int VCFXValidator::run(int argc, char* argv[]){
    bool hasStdin = !isatty(fileno(stdin));
    if(argc==1 && !hasStdin){
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
"  -s, --strict   Enable stricter checks.\n\n"
"Description:\n"
"  Validates:\n"
"   * All '##' lines are recognized as meta lines.\n"
"   * #CHROM line is present and well formed.\n"
"   * Each data line has >=8 columns, checks CHROM non-empty, POS>0,\n"
"     REF/ALT non-empty, QUAL is '.' or non-negative float, FILTER non-empty,\n"
"     INFO is minimally checked.\n"
"  In strict mode additional checks are performed:\n"
"   * Data line column count must match the #CHROM header.\n"
"   * Sample columns must match the FORMAT field structure.\n"
"   * Any warning is treated as an error.\n"
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
    if(f[0]!="#CHROM"){
        std::cerr<<"Error: #CHROM line doesn't start with '#CHROM' at line "<< lineNumber <<".\n";
        return false;
    }

    headerColumnCount = static_cast<int>(f.size());
    headerHasFormat = (headerColumnCount > 8);
    sampleCount = headerHasFormat ? headerColumnCount - 9 : 0;

    if(headerHasFormat && f[8] != "FORMAT"){
        std::string msg = "Warning: column 9 of #CHROM header is not 'FORMAT'.";
        if(strictMode){
            std::cerr << "Error: " << msg << "\n";
            return false;
        } else {
            std::cerr << msg << "\n";
        }
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
    if(headerColumnCount>0){
        if(strictMode && static_cast<int>(f.size()) != headerColumnCount){
            std::cerr << "Error: line "<<lineNumber<<" has "<<f.size()
                      <<" columns but header specifies "<<headerColumnCount<<".\n";
            return false;
        } else if(static_cast<int>(f.size()) != headerColumnCount){
            std::cerr << "Warning: line "<<lineNumber<<" column count "<<f.size()
                      <<" differs from header ("<<headerColumnCount<<").\n";
        }
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

    if(headerHasFormat){
        if(f.size()<9){
            std::cerr<<"Error: line "<<lineNumber<<" missing FORMAT/sample columns."<<"\n";
            return false;
        }
        std::vector<std::string> formatParts = split(f[8], ':');
        for(size_t i=9;i<f.size();++i){
            if(f[i]=="." || f[i].empty()) continue;
            std::vector<std::string> sampleParts = split(f[i], ':');
            if(sampleParts.size()!=formatParts.size()){
                std::string msg = "Warning: sample column " + std::to_string(i-8) +
                    " does not match FORMAT field";
                if(strictMode){
                    std::cerr<<"Error: "<<msg<<" on line "<<lineNumber<<".\n";
                    return false;
                } else {
                    std::cerr<<msg<<" on line "<<lineNumber<<".\n";
                }
            }
        }
    } else if(f.size()>8){
        std::string msg = "Warning: data line has sample columns but header lacks FORMAT";
        if(strictMode){
            std::cerr<<"Error: "<<msg<<" on line "<<lineNumber<<".\n";
            return false;
        } else {
            std::cerr<<msg<<" on line "<<lineNumber<<".\n";
        }
    }

    return true;
}

bool VCFXValidator::validateVCF(std::istream &in){
    std::string line;
    int lineNum=0;
    bool foundChromLine= false;
    std::vector<std::string> lines;

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
                // Any line starting with '#' but not "##" or "#CHROM" is invalid
                std::cerr<<"Error: line "<<lineNum<<" is a header line but neither starts with '##' nor is a #CHROM header line.\n";
                return false;
            }
        } else {
            // data line
            if(!foundChromLine){
                std::cerr<<"Error: data line encountered before #CHROM at line "<<lineNum<<".\n";
                return false;
            }
            if(!validateDataLine(line, lineNum)) return false;
        }
        lines.push_back(line);
    }
    if(!foundChromLine){
        std::cerr<<"Error: no #CHROM line found in file.\n";
        return false;
    }
    for(const auto &l : lines){
        std::cout << l << '\n';
    }
    std::cerr<<"VCF file is valid.\n";
    return true;
}

int main(int argc, char* argv[]){
    if (vcfx::handle_version_flag(argc, argv, "VCFX_validator")) return 0;
    VCFXValidator validator;
    return validator.run(argc, argv);
}
