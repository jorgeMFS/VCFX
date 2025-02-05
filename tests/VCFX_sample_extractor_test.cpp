#include <gtest/gtest.h>
#include <sstream>
#include "../src/VCFX_sample_extractor/VCFX_sample_extractor.h"

class SampleExtractorTest : public ::testing::Test {
protected:
    void SetUp() override {
        input.str("");
        output.str("");
    }

    std::stringstream input;
    std::stringstream output;
};

TEST_F(SampleExtractorTest, ExtractsSampleData) {
    input.str(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n"
        "chr1\t100\trs123\tA\tT\t40\tPASS\t.\tGT\t0/1\t1/1\n"
    );

    extractSampleData(input, output, "Sample1");
    
    std::string expected = 
        "CHROM\tPOS\tID\tREF\tALT\tSample1\n"
        "chr1\t100\trs123\tA\tT\t0/1\n";
    EXPECT_EQ(output.str(), expected);
}

TEST_F(SampleExtractorTest, HandlesNonexistentSample) {
    input.str(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
    );

    extractSampleData(input, output, "NonexistentSample");
    EXPECT_EQ(output.str(), "");
}

TEST_F(SampleExtractorTest, HandlesEmptyInput) {
    input.str("");
    extractSampleData(input, output, "Sample1");
    EXPECT_EQ(output.str(), "");
}

TEST_F(SampleExtractorTest, ParsesArguments) {
    std::string sample_name;
    const char* valid_args[] = {"program", "--sample", "Sample1"};
    EXPECT_TRUE(parseArguments(3, const_cast<char**>(valid_args), sample_name));
    EXPECT_EQ(sample_name, "Sample1");

    const char* help_args[] = {"program", "--help"};
    EXPECT_FALSE(parseArguments(2, const_cast<char**>(help_args), sample_name));
} 