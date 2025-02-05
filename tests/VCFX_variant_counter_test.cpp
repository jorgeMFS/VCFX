#include <gtest/gtest.h>
#include <sstream>
#include "../src/VCFX_variant_counter/VCFX_variant_counter.h"

class VariantCounterTest : public ::testing::Test {
protected:
    void SetUp() override {
        input.str("");
        output.str("");
    }

    std::stringstream input;
    std::stringstream output;
};

TEST_F(VariantCounterTest, CountsVariants) {
    input.str(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tT\t.\t.\t.\n"
        "chr1\t200\t.\tG\tC\t.\t.\t.\n"
    );
    
    EXPECT_EQ(countVariants(input), 2);
}

TEST_F(VariantCounterTest, HandlesEmptyInput) {
    input.str("");
    EXPECT_EQ(countVariants(input), 0);
}

TEST_F(VariantCounterTest, HandlesHeaderOnly) {
    input.str(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

    );
    EXPECT_EQ(countVariants(input), 0);
} 