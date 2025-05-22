#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <unistd.h>
#include <set>

static void print_usage(){
    std::cout << "vcfx - unified interface for VCFX tools\n"
              << "Usage: vcfx [--help] [--list] <subcommand> [args]\n\n"
              << "  <subcommand>  Name of a VCFX tool without the 'VCFX_' prefix\n"
              << "  --list        List available subcommands found in PATH\n"
              << "  --help        Show this help message\n";
}

static void list_commands(){
    const char* path_env = std::getenv("PATH");
    if(!path_env) return;
    std::string paths(path_env);
    std::set<std::string> cmds;
    size_t start=0;
    while(true){
        size_t end = paths.find(':', start);
        std::string dir = paths.substr(start, end - start);
        DIR* d = opendir(dir.c_str());
        if(d){
            struct dirent* e;
            while((e = readdir(d)) != nullptr){
                if(std::strncmp(e->d_name, "VCFX_", 5)==0){
                    std::string name = e->d_name + 5;
                    std::string full = dir + "/" + e->d_name;
                    if(access(full.c_str(), X_OK)==0){
                        cmds.insert(name);
                    }
                }
            }
            closedir(d);
        }
        if(end == std::string::npos) break;
        start = end + 1;
    }
    for(const auto& c : cmds){
        std::cout << c << '\n';
    }
}

int main(int argc, char* argv[]){
    bool show_help = false;
    bool show_list = false;
    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"list", no_argument, 0, 'l'},
        {0,0,0,0}
    };

    int opt;
    while((opt = getopt_long(argc, argv, "hl", long_opts, nullptr)) != -1){
        if(opt == 'h') show_help = true;
        else if(opt == 'l') show_list = true;
        else {
            print_usage();
            return 1;
        }
    }

    if(show_help){
        print_usage();
        return 0;
    }
    if(show_list){
        list_commands();
        return 0;
    }

    if(optind >= argc){
        print_usage();
        return 1;
    }

    std::string sub = argv[optind];
    std::string exec_name = "VCFX_" + sub;

    std::vector<char*> exec_args;
    exec_args.push_back(const_cast<char*>(exec_name.c_str()));
    for(int i = optind + 1; i < argc; ++i){
        exec_args.push_back(argv[i]);
    }
    exec_args.push_back(nullptr);

    execvp(exec_name.c_str(), exec_args.data());
    std::perror(exec_name.c_str());
    return 1;
}

