#include <gtest/gtest.h>
#include <sstream>
#include "../src/VCFX_header_parser/VCFX_header_parser.h"

TEST(HeaderParserTest, ExtractsHeaderLines) {
    std::stringstream input(
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tT\t.\t.\t.\n"
    );
    
    std::stringstream output;
    processHeader(input, output);
    
    std::string expected = 
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    
    EXPECT_EQ(output.str(), expected);
}

TEST(HeaderParserTest, HandlesEmptyInput) {
    std::stringstream input("");
    std::stringstream output;
    processHeader(input, output);
    EXPECT_EQ(output.str(), "");
}