#include <gtest/gtest.h>
#include "individual.hpp"

TEST(BHIndividualTest, CopyConstructorTest) {
    Individual ind1;

    ind1.m[0] = 0.5;

    Individual ind2(ind1);

    EXPECT_EQ(ind2.m[0],ind1.m[0]);
}

TEST(BHIndividualTest, AssignmentTest) {
    Individual ind1;
    Individual ind2;

    ind1.m[0] = 0.5;

    EXPECT_TRUE(ind1.m[0] != ind2.m[0]);
}
