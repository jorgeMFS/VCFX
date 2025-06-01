#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <getopt.h>
#include <iostream>
#include <limits.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif
#include <filesystem>
#include <fstream>
#include <set>

static const char *argv0_global = nullptr;

static void print_usage() {
    std::cout << "vcfx - unified interface for VCFX tools\n"
              << "Usage: vcfx [--help] [--list] <subcommand> [args]\n\n"
              << "  <subcommand>  Name of a VCFX tool without the 'VCFX_' prefix\n"
              << "  list          Alias for --list\n"
              << "  help <tool>   Show Markdown documentation for a tool if available\n"
              << "  --list        List available subcommands found in PATH\n"
              << "  --help        Show this help message\n";
}

static void list_commands() {
    const char *path_env = std::getenv("PATH");
    if (!path_env)
        return;
    std::string paths(path_env);
    std::set<std::string> cmds;
    size_t start = 0;
    while (true) {
        size_t end = paths.find(':', start);
        std::string dir = paths.substr(start, end - start);
        DIR *d = opendir(dir.c_str());
        if (d) {
            struct dirent *e;
            while ((e = readdir(d)) != nullptr) {
                if (std::strncmp(e->d_name, "VCFX_", 5) == 0) {
                    std::string name = e->d_name + 5;
                    std::string full = dir + "/" + e->d_name;
                    if (access(full.c_str(), X_OK) == 0) {
                        cmds.insert(name);
                    }
                }
            }
            closedir(d);
        }
        if (end == std::string::npos)
            break;
        start = end + 1;
    }
    for (const auto &c : cmds) {
        std::cout << c << '\n';
    }
}

static std::vector<std::string> get_doc_dirs() {
    std::vector<std::string> dirs;
    const char *env = std::getenv("VCFX_DOCS_DIR");
    if (env)
        dirs.emplace_back(env);

    std::string exe;
    char buf[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (len > 0) {
        buf[len] = '\0';
        exe = buf;
    }
#ifdef __APPLE__
    if (exe.empty()) {
        uint32_t size = sizeof(buf);
        if (_NSGetExecutablePath(buf, &size) == 0) {
            exe = buf;
        }
    }
#endif
    if (exe.empty() && argv0_global) {
        exe = argv0_global;
    }

    if (!exe.empty()) {
        auto pos = exe.find_last_of('/');
        if (pos != std::string::npos) {
            std::string base = exe.substr(0, pos);
            dirs.push_back(base + "/../share/doc/VCFX");
            dirs.push_back(base + "/../share/vcfx/docs");
            dirs.push_back(base + "/../docs");
            dirs.push_back(base + "/../../docs");
            dirs.push_back(base + "/../../../docs"); // when run from build/src/vcfx_wrapper
        }
    }
    dirs.push_back("docs");
    return dirs;
}

static int print_tool_doc(const std::string &tool) {
    std::string fname = "VCFX_" + tool + ".md";
    namespace fs = std::filesystem;
    for (const auto &dir : get_doc_dirs()) {
        fs::path base(dir);
        fs::path direct = base / fname;
        std::ifstream in(direct);
        if (in) {
            std::cout << in.rdbuf();
            return 0;
        }
        if (fs::exists(base) && fs::is_directory(base)) {
            for (const auto &entry : fs::recursive_directory_iterator(base)) {
                if (entry.path().filename() == fname && entry.is_regular_file()) {
                    std::ifstream rin(entry.path());
                    if (rin) {
                        std::cout << rin.rdbuf();
                        return 0;
                    }
                }
            }
        }
    }
    std::cerr << "Documentation for '" << tool << "' not found." << std::endl;
    return 1;
}

int main(int argc, char *argv[]) {
    argv0_global = argv[0];
    bool show_help = false;
    bool show_list = false;
    static struct option long_opts[] = {{"help", no_argument, 0, 'h'}, {"list", no_argument, 0, 'l'}, {0, 0, 0, 0}};

    int opt;
    while ((opt = getopt_long(argc, argv, "hl", long_opts, nullptr)) != -1) {
        if (opt == 'h')
            show_help = true;
        else if (opt == 'l')
            show_list = true;
        else {
            print_usage();
            return 1;
        }
    }

    if (show_help) {
        print_usage();
        return 0;
    }
    if (show_list) {
        list_commands();
        return 0;
    }

    if (optind >= argc) {
        print_usage();
        return 1;
    }

    std::string sub = argv[optind];

    if (sub == "list") {
        list_commands();
        return 0;
    }

    if (sub == "help") {
        if (optind + 1 >= argc) {
            print_usage();
            return 0;
        }
        return print_tool_doc(argv[optind + 1]);
    }

    std::string exec_name = "VCFX_" + sub;

    std::vector<char *> exec_args;
    exec_args.push_back(const_cast<char *>(exec_name.c_str()));
    for (int i = optind + 1; i < argc; ++i) {
        exec_args.push_back(argv[i]);
    }
    exec_args.push_back(nullptr);

    execvp(exec_name.c_str(), exec_args.data());
    std::perror(exec_name.c_str());
    return 1;
}
