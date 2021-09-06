#include <functional>
#include <gtest/gtest.h>

#include "../bisection.cpp"

TEST(BisectionTest, BasicAssertions) {

    std::function<float(const float&)> f = [](float x) { return x; };

    float result = bisection(f, -10, 10);

    EXPECT_EQ(result, 0.0);
}


TEST(NewtonsTest, BasicAssertions) {

}

TEST(RegularFalsiTest, BasicAssertions) {

}

TEST(SecantTest, BasicAssertions) {
    
}