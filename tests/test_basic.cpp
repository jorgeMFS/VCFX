#include <cstdlib>
#include <iostream>

int main(){
    // We do a minimal check that we can run a tool's --help
    int ret = system("VCFX_variant_counter --help");
    if(ret != 0){
        std::cerr << "VCFX_variant_counter --help returned non-zero!\n";
        return 1;
    }

    std::cout << "test_basic passed.\n";
    return 0;
}
