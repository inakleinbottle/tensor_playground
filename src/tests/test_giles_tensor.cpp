//
// Created by sam on 20/02/2022.
//

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


