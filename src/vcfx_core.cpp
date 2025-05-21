#include "vcfx_core.h"
#include <algorithm>
#include <cctype>
#include <sstream>

namespace vcfx {

std::string trim(const std::string& str) {
    auto first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) {
        return "";
    }
    auto last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> result;
    std::istringstream iss(str);
    std::string item;
    while (std::getline(iss, item, delimiter)) {
        result.push_back(item);
    }
    return result;
}

void print_error(const std::string& msg, std::ostream& os) {
    os << "Error: " << msg << '\n';
}

void print_version(const std::string& tool, const std::string& version,
                   std::ostream& os) {
    os << tool << " version " << version << '\n';
}

}  // namespace vcfx
