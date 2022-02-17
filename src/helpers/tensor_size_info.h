//
// Created by sam on 07/11/2021.
//

#ifndef TENSOR_PLAYGROUND_TENSOR_SIZE_INFO_H
#define TENSOR_PLAYGROUND_TENSOR_SIZE_INFO_H


namespace playground {
    constexpr size_t power(size_t base, unsigned exponent) noexcept {
        return (exponent == 0) ? 1 :
               (exponent == 1) ? base :
               ((exponent % 2 == 0) ? 1 : base) * power(base, exponent / 2) * power(base, exponent / 2);
    }

    constexpr size_t tensor_alg_size(size_t Width, unsigned degree) noexcept {
        return (power(Width, degree + 1) - 1) / (Width - 1);
    }

    constexpr size_t tensor_start_of_degree(size_t width, unsigned degree) noexcept {
        return (power(width, degree) - 1) / (width - 1);
    }

}

#endif //TENSOR_PLAYGROUND_TENSOR_SIZE_INFO_H
