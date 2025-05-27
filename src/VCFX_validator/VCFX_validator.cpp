#include "vcfx_core.h"
#include "VCFX_validator.h"
#include <getopt.h>
#include <sstream>
#include <vector>
#include <cctype>
#include <cstdlib>
#include <unistd.h>
#include <regex>

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

static bool valid_dna(const std::string& s){
    if(s.empty()) return false;
    for(char c: s){
        char u = std::toupper(static_cast<unsigned char>(c));
        if(u!='A' && u!='C' && u!='G' && u!='T' && u!='N') return false;
    }
    return true;
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
        {"report-dups", no_argument, 0, 'd'},
        {0,0,0,0}
    };
    while(true){
        int c= ::getopt_long(argc, argv, "hsd", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h': showHelp=true; break;
            case 's': strictMode= true; break;
            case 'd': reportDuplicates = true; break;
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
"  -s, --strict   Enable stricter checks.\n"
"  -d, --report-dups Report duplicate records.\n\n"
"Description:\n"
"  Validates:\n"
"   * All '##' lines are recognized as meta lines.\n"
"   * #CHROM line is present and well formed.\n"
"   * Each data line has >=8 columns, checks CHROM non-empty, POS>0,\n"
"     REF/ALT non-empty, QUAL is '.' or non-negative float, FILTER non-empty,\n"
"     INFO fields are checked against header definitions.\n"
"   * FORMAT fields and genotype values are validated.\n"
"   * REF and ALT sequences must contain only A,C,G,T,N.\n"
"   * Duplicate variants are detected when --report-dups is used.\n"
"  In strict mode additional checks are performed:\n"
"   * Data line column count must match the #CHROM header.\n"
"   * Sample columns must match the FORMAT field structure.\n"
"   * Any warning is treated as an error.\n"
"  Exits 0 if pass, 1 if fail.\n";
}

bool VCFXValidator::validateMetaLine(const std::string &line, int lineNumber){
    if(line.size()<2) return false;
    if(line.rfind("##INFO=",0)==0 || line.rfind("##FORMAT=",0)==0){
        size_t start=line.find('<');
        size_t end=line.rfind('>');
        if(start==std::string::npos || end==std::string::npos || end<=start){
            std::cerr<<"Error: malformed header at line "<<lineNumber<<".\n";
            return false;
        }
        std::string inside=line.substr(start+1,end-start-1);
        auto fields=split(inside,',');
        std::string id,number,type;
        for(auto &f:fields){
            auto eq=f.find('=');
            if(eq==std::string::npos) continue;
            std::string key=trim(f.substr(0,eq));
            std::string val=trim(f.substr(eq+1));
            if(key=="ID") id=val;
            else if(key=="Number") number=val;
            else if(key=="Type") type=val;
        }
        if(id.empty()){
            std::cerr<<"Error: header line missing ID at "<<lineNumber<<".\n";
            return false;
        }
        FieldDef def{number,type};
        if(line.rfind("##INFO=",0)==0){
            infoDefs[id]=def;
        }else{
            formatDefs[id]=def;
        }
        return true;
    }
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
    if(!valid_dna(f[3])){
        std::cerr<<"Error: line "<<lineNumber<<" REF has invalid characters.\n";
        return false;
    }
    // ALT
    if(f[4].empty()){
        std::cerr<<"Error: line "<< lineNumber<<" ALT is empty.\n";
        return false;
    }
    if(!valid_dna(f[4])){
        std::cerr<<"Error: line "<<lineNumber<<" ALT has invalid characters.\n";
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
    // INFO field validation
    if(f[7]!="."){
        std::stringstream infoSS(f[7]);
        std::string token;
        bool anyValid=false;
        while(std::getline(infoSS, token, ';')){
            token=trim(token);
            if(token.empty()) continue;
            auto eq=token.find('=');
            if(eq==std::string::npos){
                std::string k=token;
                auto it=infoDefs.find(k);
                if(it==infoDefs.end()){
                    std::string msg="INFO field " + k + " not defined in header";
                    if(strictMode){
                        std::cerr<<"Error: "<<msg<<" on line "<<lineNumber<<".\n";
                        return false;
                    }else{
                        std::cerr<<"Warning: "<<msg<<" on line "<<lineNumber<<".\n";
                    }
                }
                anyValid=true;
            }else{
                std::string k=token.substr(0,eq);
                std::string v=token.substr(eq+1);
                if(k.empty()){
                    std::cerr<<"Error: line "<<lineNumber<<" has INFO with empty key.\n";
                    return false;
                }
                auto it=infoDefs.find(k);
                if(it==infoDefs.end()){
                    std::string msg="INFO field " + k + " not defined in header";
                    if(strictMode){
                        std::cerr<<"Error: "<<msg<<" on line "<<lineNumber<<".\n";
                        return false;
                    }else{
                        std::cerr<<"Warning: "<<msg<<" on line "<<lineNumber<<".\n";
                    }
                }else if(!it->second.number.empty() &&
                         std::all_of(it->second.number.begin(), it->second.number.end(), ::isdigit)){
                    size_t expected=static_cast<size_t>(std::stoi(it->second.number));
                    size_t have=split(v, ',').size();
                    if(have!=expected){
                        std::string msg="INFO field " + k + " expected " + it->second.number +
                                         " values";
                        if(strictMode){
                            std::cerr<<"Error: "<<msg<<" on line "<<lineNumber<<".\n";
                            return false;
                        }else{
                            std::cerr<<"Warning: "<<msg<<" on line "<<lineNumber<<".\n";
                        }
                    }
                }
                anyValid=true;
            }
        }
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
        for(const auto &fp : formatParts){
            auto it = formatDefs.find(fp);
            if(it==formatDefs.end()){
                std::string msg = "FORMAT field " + fp + " not defined in header";
                if(strictMode){
                    std::cerr<<"Error: "<<msg<<" on line "<<lineNumber<<".\n";
                    return false;
                } else {
                    std::cerr<<"Warning: "<<msg<<" on line "<<lineNumber<<".\n";
                }
            }
        }
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
            for(size_t j=0;j<sampleParts.size() && j<formatParts.size();++j){
                const std::string &key = formatParts[j];
                const std::string &val = sampleParts[j];
                auto it = formatDefs.find(key);
                if(it!=formatDefs.end() && !it->second.number.empty() &&
                   std::all_of(it->second.number.begin(), it->second.number.end(), ::isdigit)){
                    size_t expected = static_cast<size_t>(std::stoi(it->second.number));
                    size_t have = split(val, ',').size();
                    if(have!=expected){
                        std::string msg = "FORMAT field " + key + " expected " + it->second.number + " values";
                        if(strictMode){
                            std::cerr<<"Error: "<<msg<<" on line "<<lineNumber<<".\n";
                            return false;
                        } else {
                            std::cerr<<"Warning: "<<msg<<" on line "<<lineNumber<<".\n";
                        }
                    }
                }
                if(key=="GT" && !val.empty()){
                    static const std::regex gtRegex("^(\\.|(\\d+([/|]\\d+)*))$");
                    if(!std::regex_match(val, gtRegex)){
                        std::string msg = "invalid genotype";
                        if(strictMode){
                            std::cerr<<"Error: "<<msg<<" on line "<<lineNumber<<".\n";
                            return false;
                        } else {
                            std::cerr<<"Warning: "<<msg<<" on line "<<lineNumber<<".\n";
                        }
                    }
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

    // duplicate detection
    std::string key = f[0] + ":" + f[1] + ":" + f[3] + ":" + f[4];
    if(seenVariants.count(key)){
        std::string msg = "duplicate variant";
        if(reportDuplicates){
            std::cerr<<"Duplicate at line "<<lineNumber<<"\n";
        }
        if(strictMode){
            std::cerr<<"Error: "<<msg<<" on line "<<lineNumber<<".\n";
            return false;
        } else {
            std::cerr<<"Warning: "<<msg<<" on line "<<lineNumber<<".\n";
        }
    } else {
        seenVariants.insert(key);
    }

    return true;
}

bool VCFXValidator::validateVCF(std::istream &in){
    std::string data;
    if(!vcfx::read_maybe_compressed(in, data)){
        std::cerr<<"Error: failed to read input stream.\n";
        return false;
    }
    std::istringstream iss(data);
    std::string line;
    int lineNum=0;
    bool foundChromLine=false;
    std::vector<std::string> lines;

    while(true){
        if(!std::getline(iss, line)) break;
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

static void show_help() { VCFXValidator obj; char arg0[] = "VCFX_validator"; char arg1[] = "--help"; char* argv2[] = {arg0, arg1, nullptr}; obj.run(2, argv2); }

int main(int argc, char* argv[]){
    if (vcfx::handle_common_flags(argc, argv, "VCFX_validator", show_help)) return 0;
    VCFXValidator validator;
    return validator.run(argc, argv);
}
