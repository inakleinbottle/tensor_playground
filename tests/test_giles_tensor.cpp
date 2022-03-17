//
// Created by sam on 20/02/2022.
//
#include <stdio.h>

#include <gtest/gtest.h>

#include <simple_template_tensor.h>
#include <giles_template_tensor.h>

using namespace playground;

////////////////////// IDENTITY MULTIPLICATION TESTS ///////////////////////

TEST(test_giles_tensor, identity_mul_fixed) {
    constexpr unsigned width = 5;
    constexpr unsigned depth = 5;

    giles_template_tensor<width, depth> identity {1.0};
    giles_template_tensor<width, depth> other;

    for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) {
        other[i] = double(i);
    }

//    EXPECT_EQ(other*identity, other);
//    EXPECT_EQ(identity*other, other);

    giles_template_tensor<width, depth> zero;
    EXPECT_EQ(identity*other - other, zero);
}

TEST(test_giles_tensor, fixed_mul_identity) {
    constexpr unsigned width = 5;
    constexpr unsigned depth = 5;

    giles_template_tensor<width, depth> identity {1.0};
    giles_template_tensor<width, depth> other;

    for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) {
        other[i] = double(i);
    }

//    EXPECT_EQ(other*identity, other);
//    EXPECT_EQ(identity*other, other);

    giles_template_tensor<width, depth> zero;
    EXPECT_EQ(other*identity - other, zero);
}

////////////////////// INCREASING DEPTH TESTS ///////////////////////

class GilesTensorIntegerCoefficientsTests : public ::testing::Test
{
protected:

    static constexpr unsigned width = 5;
    static constexpr unsigned depth = 5;

    giles_template_tensor<width, depth> giles_result;
    simple_template_tensor<width, depth> simple_result;

    virtual void SetUp()
    {
        // giles tensor lhs

        giles_template_tensor<width, depth> giles_lhs;

        for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) {
            giles_lhs[i] = double(i);
        }

        // simple tensor lhs

        simple_template_tensor<width, depth> simple_lhs;

        for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) {
            simple_lhs[i] = double(i);
        }

        // giles tensor rhs

        giles_template_tensor<width, depth> giles_rhs;

        for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) {
            giles_rhs[i] = double(i);
        }

        // simple tensor rhs

        simple_template_tensor<width, depth> simple_rhs;

        for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) {
            simple_rhs[i] = double(i);
        }

        // ------- Multipy -------- //

        giles_result = giles_lhs*giles_rhs;

        simple_result = simple_lhs*simple_rhs;
    }
};

TEST_F(GilesTensorIntegerCoefficientsTests, DepthOneTest)
{
    constexpr unsigned depth_to_check = 1;

    for (size_t i = tensor_alg_size(width, depth_to_check-1); i < tensor_alg_size(width, depth_to_check); ++i) 
    {

        EXPECT_EQ(giles_result[i], simple_result[i]) << "Multiplication result differs at index " << i;
    }

}

TEST_F(GilesTensorIntegerCoefficientsTests, DepthTwoTest)
{
    constexpr unsigned depth_to_check = 2;

    for (size_t i = tensor_alg_size(width, depth_to_check-1); i < tensor_alg_size(width, depth_to_check); ++i) 
    {

        EXPECT_EQ(giles_result[i], simple_result[i]) << "Multiplication result differs at index " << i;
    }

}

TEST_F(GilesTensorIntegerCoefficientsTests, DepthThreeTest)
{
    constexpr unsigned depth_to_check = 3;

    for (size_t i = tensor_alg_size(width, depth_to_check-1); i < tensor_alg_size(width, depth_to_check); ++i) 
    {

        EXPECT_EQ(giles_result[i], simple_result[i]) << "Multiplication result differs at index " << i;
    }

}

TEST_F(GilesTensorIntegerCoefficientsTests, DepthFourTest)
{
    constexpr unsigned depth_to_check = 4;

    for (size_t i = tensor_alg_size(width, depth_to_check-1); i < tensor_alg_size(width, depth_to_check); ++i) 
    {

        EXPECT_EQ(giles_result[i], simple_result[i]) << "Multiplication result differs at index " << i;
    }

}

TEST_F(GilesTensorIntegerCoefficientsTests, DepthFiveTest)
{
    constexpr unsigned depth_to_check = 5;

    for (size_t i = tensor_alg_size(width, depth_to_check-1); i < tensor_alg_size(width, depth_to_check); ++i) 
    {

        EXPECT_EQ(giles_result[i], simple_result[i]) << "Multiplication result differs at index " << i;
    }

}

TEST_F(GilesTensorIntegerCoefficientsTests, FullDepthTest)
{
    for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) 
    {

        EXPECT_EQ(giles_result[i], simple_result[i]) << "Multiplication result differs at index " << i;
    }

}

////////////////////// TILING TESTS ///////////////////////

class GilesTensorDepthFourTilingTests : public ::testing::Test
{
protected:

    static constexpr unsigned width = 5;
    static constexpr unsigned depth = 5;

    void get_contributions_from_degrees(
        unsigned lhs_degree, 
        unsigned rhs_degree, 
        giles_template_tensor<width, depth> &giles_result,
        simple_template_tensor<width, depth> &simple_result
        )
    {
        giles_template_tensor<width, depth> giles_lhs;
        giles_template_tensor<width, depth> giles_rhs;
        simple_template_tensor<width, depth> simple_lhs;
        simple_template_tensor<width, depth> simple_rhs;

        giles_lhs[tensor_alg_size(width, lhs_degree-1)] = 1.0;
        simple_lhs[tensor_alg_size(width, lhs_degree-1)] = 1.0;
        giles_rhs[tensor_alg_size(width, rhs_degree-1)] = 1.0;
        simple_rhs[tensor_alg_size(width, rhs_degree-1)] = 1.0;

        giles_result = giles_lhs*giles_rhs;
        simple_result = simple_lhs*simple_rhs;
    }

    virtual void SetUp()
    {

    }
};

TEST_F(GilesTensorDepthFourTilingTests, OneThreeTest)
{
    // lhs = {0 1{1} 0 ... 0}
    // rhs = {0 ... 0 1{111} 0 ... 0}
    // ans = {0 ... 0 1{1111} 0 ... 0}

    giles_template_tensor<width, depth> giles_result;
    simple_template_tensor<width, depth> simple_result;

    get_contributions_from_degrees(1, 3, giles_result, simple_result);
    
    for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) 
    {
        EXPECT_EQ(giles_result[i], simple_result[i]) << "Multiplication result differs at index " << i;
    }
}

TEST_F(GilesTensorDepthFourTilingTests, TwoTwoTest)
{
    // lhs = {0 ... 0 1{11} 0 ... 0}
    // rhs = {0 ... 0 1{11} 0 ... 0}
    // ans = {0 ... 0 1{1111} 0 ... 0}

    giles_template_tensor<width, depth> giles_result;
    simple_template_tensor<width, depth> simple_result;

    get_contributions_from_degrees(2, 2, giles_result, simple_result);
    
    for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) 
    {
        EXPECT_EQ(giles_result[i], simple_result[i]) << "Multiplication result differs at index " << i;
    }

}

TEST_F(GilesTensorDepthFourTilingTests, ThreeOneTest)
{
    // lhs = {0 ... 0 1{111} 0 ... 0}
    // rhs = {0 ... 0 1{1} 0 ... 0}
    // ans = {0 ... 0 1{1111} 0 ... 0}

    giles_template_tensor<width, depth> giles_result;
    simple_template_tensor<width, depth> simple_result;

    get_contributions_from_degrees(3, 1, giles_result, simple_result);

    for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) 
    {
        EXPECT_EQ(giles_result[i], simple_result[i]) << "Multiplication result differs at index " << i;
    }

}