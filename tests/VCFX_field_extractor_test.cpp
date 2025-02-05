#include <gtest/gtest.h>
#include "../src/VCFX_field_extractor/VCFX_field_extractor.h"

TEST(FieldExtractorTest, ExtractsBasicFields) {
    std::string record = "chr1\t100\trs123\tA\tT\t40\tPASS\tDP=150;AF=0.1";
    std::vector<std::string> fields = {"CHROM", "POS", "ID", "QUAL"};
    
    auto extracted = parseFields(record, fields);
    EXPECT_EQ(extracted.size(), 4);
    EXPECT_EQ(extracted[0], "chr1");
    EXPECT_EQ(extracted[1], "100");
    EXPECT_EQ(extracted[2], "rs123");
    EXPECT_EQ(extracted[3], "40");
}

TEST(FieldExtractorTest, ExtractsInfoFields) {
    std::string record = "chr1\t100\trs123\tA\tT\t40\tPASS\tDP=150;AF=0.1";
    std::vector<std::string> fields = {"DP", "AF"};
    
    auto extracted = parseFields(record, fields);
    EXPECT_EQ(extracted.size(), 2);
    EXPECT_EQ(extracted[0], "150");
    EXPECT_EQ(extracted[1], "0.1");
}

TEST(FieldExtractorTest, HandlesNonExistentFields) {
    std::string record = "chr1\t100\trs123\tA\tT\t40\tPASS\tDP=150";
    std::vector<std::string> fields = {"CHROM", "NONEXISTENT", "DP"};
    
    auto extracted = parseFields(record, fields);
    EXPECT_EQ(extracted.size(), 3);
    EXPECT_EQ(extracted[0], "chr1");
    EXPECT_EQ(extracted[1], ".");
    EXPECT_EQ(extracted[2], "150");
}