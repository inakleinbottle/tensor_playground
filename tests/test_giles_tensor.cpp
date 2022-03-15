//
// Created by sam on 20/02/2022.
//
#include <stdio.h>

#include <gtest/gtest.h>

#include <simple_template_tensor.h>
#include <giles_template_tensor.h>


using namespace playground;

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

TEST(test_giles_tensor, integer_giles_tensor_simple_tensor_compare_multiplication)
{
    constexpr unsigned width = 5;
    constexpr unsigned depth = 5;

    // -- lhs -- //

    // giles tensor lhs

    giles_template_tensor<width, depth> giles_lhs;

    for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) {
        giles_lhs[i] = double(i);
    }

    // std::cout << "giles_lhs=" << giles_lhs << std::endl;

    // simple tensor lhs

    simple_template_tensor<width, depth> simple_lhs;

    for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) {
        simple_lhs[i] = double(i);
    }

    // std::cout << "simple_lhs=" << simple_lhs << std::endl;

    // -- rhs -- //

    // giles tensor rhs

    giles_template_tensor<width, depth> giles_rhs;

    for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) {
        giles_rhs[i] = double(i);
    }

    // std::cout << "giles_rhs=" << giles_rhs << std::endl;

    // simple tensor rhs

    simple_template_tensor<width, depth> simple_rhs;

    for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) {
        simple_rhs[i] = double(i);
    }

    // std::cout << "simple_rhs=" << simple_rhs << std::endl;

    // ------- Multipy -------- //

    giles_template_tensor<width, depth> giles_result = giles_lhs*giles_rhs;

    simple_template_tensor<width, depth> simple_result = simple_lhs*simple_rhs;

    // ------- Compare --------- //

    for (size_t i = 0; i < tensor_alg_size(width, depth); ++i) 
    {

        EXPECT_EQ(giles_result[i], simple_result[i]) << "Multiplication result differs at index " << i;
    }

}


