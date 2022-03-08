//
// Created by sam on 14/02/2022.
//

#ifndef TENSOR_PLAYGROUND_GILES_TEMPLATE_TENSOR_H
#define TENSOR_PLAYGROUND_GILES_TEMPLATE_TENSOR_H

//#define TENSOR_PLAYGROUND_INLINE_DEF

#include "index_word.h"
#include "mul_level_func.h"
#include <simple_template_tensor.h>

namespace playground {

template<unsigned Width, unsigned Depth, typename Coeffs = double>
class giles_template_tensor : public simple_template_tensor<Width, Depth, Coeffs> {
    using base_type = simple_template_tensor<Width, Depth, Coeffs>;
    using coeff_container_t = typename base_type::base_type;

   public:
    using base_type::operator[];
    using base_type::assign;
    using base_type::begin;
    using base_type::data;
    using base_type::emplace_back;
    using base_type::end;
    using base_type::size;
    using typename base_type::size_type;
    using degree_type = unsigned;
    using base_type::size_array;

    static constexpr degree_type max_depth = Depth;

    using index_key = index_word<Width, size_type>;

    static constexpr degree_type tile_letters = 2;
    static constexpr size_type tile_width = power(Width, tile_letters);
    static constexpr size_type tile_size = tile_width * tile_width;

    using base_type::operator=;

    coeff_container_t reverse_data;

   private:
    void fill_reverse_data() {
        reverse_data.resize(tensor_alg_size(Width, Depth-1));
        for (auto k = index_key(0); k < size_array[Depth - 1]; ++k) {
            reverse_data[k.reverse()] = operator[](k);
        }
    }

   public:
    giles_template_tensor() : base_type() , reverse_data(tensor_alg_size(Width, Depth-1)){
    }

    giles_template_tensor(std::initializer_list<Coeffs> args) : base_type(args) {
        fill_reverse_data();
    }

    explicit giles_template_tensor(size_t index, Coeffs coeff = 1.0) : base_type(index, coeff) {
        fill_reverse_data();
    }

    template<typename... Args>
    explicit giles_template_tensor(Args... args) : base_type(std::forward<Args>(args)...) {
        fill_reverse_data();
    }

    using base_type::range_begin;
    using base_type::range_end;

   private:
    template<typename Op>
    static void
    inplace_multiplication_impl(
        giles_template_tensor &out,
        const giles_template_tensor &lhs,
        const giles_template_tensor &rhs,
        const std::vector<size_type> &start_of_degree_array,
        Op op) noexcept {

        assert(lhs.size() > 0);
        assert(rhs.size() > 0);

#ifdef TENSOR_PLAYGROUND_INLINE_DEF
        Coeffs tile[tile_size];

        for (degree_type out_depth = Depth; out_depth > 2 * tile_letters; --out_depth) {
            const auto stride = power(Width, out_depth - tile_letters);

            for (auto k = 0; k < power(Width, out_depth); ++k) {

                // Handle cases of 0*out_depth and out_depth*0
                {
                    auto lhs_unit = lhs[0], rhs_unit = rhs[0];
                    const auto *lhs_ptr = lhs.range_begin() + k + tensor_start_of_degree(Width, out_depth);
                    const auto *rhs_ptr = rhs.range_begin() + k + tensor_start_of_degree(Width, out_depth);

                    /// First do the zero*Depth and Depth*zero case
                    for (size_type i = 0; i < tile_width; ++i) {
                        for (size_type j = 0; j < tile_width; ++j) {
                            tile[i * tile_width + j] = lhs_unit * rhs_ptr[i * stride + j] + lhs_ptr[i * stride + j] * rhs_unit;
                        }
                    }
                }
                {
                    // Handle the cases where the left hand side degree is too small for tiling
                    const auto *rhs_ptr = rhs.range_begin();
                    /*
                     * Here we have to "transfer" the excess letters to the right-hand word to offset
                     * the fact that on the left-hand side we have at most lhs_deg letters.
                     * For example, for lhs_deg=1 we take the first letter from the tile index i for the
                     * left hand side, and use the remaining letters from index i are taken to the right
                     * hand key
                     */
                    for (degree_type lhs_deg = 1; lhs_deg < tile_letters; ++lhs_deg) {
                        const auto *lhs_ptr = lhs.reverse_data.data() + size_array[lhs_deg];
                        for (size_type i = 0; i < tile_width; ++i) {
                            index_key tmp(i + size_array[lhs_deg]);
                            auto lhs_word = k / power(Width, tile_letters - lhs_deg);
                            auto lhs_middle = k % power(Width, tile_letters - lhs_deg);

//                            auto split = tmp.split(tile_letters - lhs_deg);// split is by number of right-hand letters
//                            auto lhs_word = static_cast<size_type>(split.first) - tensor_start_of_degree(Width, out_depth-lhs_deg);
//                            auto lhs_middle = static_cast<size_type>(split.second) - tensor_start_of_degree(Width, lhs_deg);

                            const auto &lhs_val = lhs_ptr[static_cast<size_type>(lhs_word)];
                            for (size_type j = 0; j < tile_width; ++j) {
                                // We've already accounted for the offset induced by the middle word k.
                                auto t = static_cast<size_type>(lhs_middle) * tile_width + j;
                                tile[i*tile_width + j] += lhs_val * rhs_ptr[t];
                            }
                        }
                    }
                }

                // Middle cases
                for (degree_type split_deg = 1; split_deg < out_depth - 2 * tile_letters - 1; ++split_deg) {
                    index_key tmp(k + size_array[out_depth]);
                    auto split = tmp.split(split_deg);
                    auto lhs_middle = static_cast<size_type>(split.first.reverse()) - size_array[out_depth-split_deg];
                    auto rhs_middle = static_cast<size_type>(split.second) - size_array[split_deg];

                    for (size_type i = 0; i < tile_width; ++i) {
                        for (size_type j = 0; j < tile_width; ++j) {
                            auto lhs_idx = lhs_middle * tile_width + i;
                            auto rhs_idx = rhs_middle * tile_width + j;
                            tile[i * tile_width + j] += lhs.reverse_data[lhs_idx] * rhs[rhs_idx];
                        }
                    }
                }

                // Handle cases where the right hand side degree is too small for tiling
                {
                    const auto *lhs_ptr = lhs.reverse_data.data();
                    /*
                     * This is the reverse of the above. We have to shift the letters from the right-hand
                     * tile index over to the left to account for the fact that there aren't enough letters
                     * on the right to make a full tile.
                     */
                    for (degree_type rhs_deg = 1; rhs_deg < tile_letters; ++rhs_deg) {
                        const auto *rhs_ptr = rhs.range_begin() + size_array[rhs_deg - 1];
                        for (size_type i = 0; i < tile_width; ++i) {
                            for (size_type j = 0; j < tile_width; ++j) {
//                                index_key tmp(j + size_array[rhs_deg]);
//                                auto split = tmp.split(rhs_deg);
//                                auto rhs_middle = split.first.reverse();
//                                auto rhs_word = split.second;
                                auto rhs_middle =
                                tile[static_cast<size_type>(i) * tile_width + static_cast<size_type>(j)] += lhs_ptr[rhs_middle * tile_width + i] * rhs_ptr[static_cast<size_type>(rhs_word)];
                            }
                        }
                    }
                }

                auto sod = tensor_start_of_degree(Width, out_depth);
                // Write the tile out into the result.
                auto *out_ptr = out.range_begin() + static_cast<size_type>(k) + sod;
                for (size_type i = 0; i < tile_width; ++i) {
                    for (size_type j = 0; j < tile_width; ++j) {
                        auto t2 = k*tile_width;
                        auto tmp = (i * stride + j + k);
                        assert(tmp < tensor_alg_size(Width, out_depth));
                        out_ptr[i * stride + j] += tile[i * tile_width + j];
                    }
                }
            }
        }

        /*
         * Now we have to handle all of the lower degree terms, which are supposedly small enough
         * to fit in cache anyway. This is the same as the old FMA code.
         */
        int rhs_deg;

        double *out_ptr;
        const double *lhs_p, *rhs_p;
        for (int out_deg = 2 * tile_letters; out_deg >= 0; --out_deg) {

            for (int lhs_deg = out_deg; lhs_deg >= 0; --lhs_deg) {
                rhs_deg = out_deg - lhs_deg;

                out_ptr = out.data() + start_of_degree_array[out_deg];
                lhs_p = lhs.range_begin() + start_of_degree_array[lhs_deg];
                rhs_p = rhs.range_begin() + start_of_degree_array[rhs_deg];

                for (int i = 0; i < start_of_degree_array[lhs_deg + 1] - start_of_degree_array[lhs_deg]; ++i) {
                    for (int j = 0; j < start_of_degree_array[rhs_deg + 1] - start_of_degree_array[rhs_deg]; ++j) {
                        *(out_ptr++) += op(lhs_p[i] * rhs_p[j]);
                    }
                }
            }
        }

#else
        dtl::multiplication_level_helper<Coeffs, Width, 2, 2>::template do_level<Depth>(out.range_begin(), lhs, rhs, op);
#endif
    }

   public:
    giles_template_tensor &operator+=(const giles_template_tensor &rhs) {
        double *p = range_begin();
        const double *rp = rhs.range_begin();
        for (size_t i = 0; i < tensor_alg_size(Width, Depth); ++i) {
            *(p++) += *(rp++);
        }
        return *this;
    }

    giles_template_tensor &operator-=(const giles_template_tensor &rhs) {
        double *p = range_begin();
        const double *rp = rhs.range_begin();
        for (size_t i = 0; i < tensor_alg_size(Width, Depth); ++i) {
            *(p++) -= *(rp++);
        }
        return *this;
    }

    giles_template_tensor &operator*=(double rhs) {
        double *p = range_begin();
        ;
        for (size_t i = 0; i < tensor_alg_size(Width, Depth); ++i) {
            *(p++) *= rhs;
        }
        return *this;
    }

    giles_template_tensor &operator/=(double rhs) {
        double *p = range_begin();
        ;
        for (size_t i = 0; i < tensor_alg_size(Width, Depth); ++i) {
            *(p++) /= rhs;
        }
        return *this;
    }

    giles_template_tensor operator+(const giles_template_tensor &rhs) const {
        giles_template_tensor result(*this);
        result += rhs;
        return result;
    }

    giles_template_tensor operator-(const giles_template_tensor &rhs) const {
        giles_template_tensor result(*this);
        result -= rhs;
        return result;
    }

    giles_template_tensor operator*(double rhs) const {
        giles_template_tensor result(*this);
        result *= rhs;
        return result;
    }

    giles_template_tensor operator/(double rhs) const {
        giles_template_tensor result(*this);
        result /= rhs;
        return result;
    }

    giles_template_tensor &add_scal_prod(const giles_template_tensor &rhs, double sca) {
        double *p = range_begin();
        const double *rp = rhs.range_begin();
        for (size_t i = 0; i < tensor_alg_size(Width, Depth); ++i) {
            *(p++) += *(rp++) * sca;
        }
        return *this;
    }

    giles_template_tensor &sub_scal_prod(const giles_template_tensor &rhs, double sca) {
        double *p = range_begin();
        const double *rp = rhs.range_begin();
        for (size_t i = 0; i < tensor_alg_size(Width, Depth); ++i) {
            *(p++) -= *(rp++) * sca;
        }
        return *this;
    }

    giles_template_tensor &add_scal_div(const giles_template_tensor &rhs, double sca) {
        double *p = range_begin();
        const double *rp = rhs.range_begin();
        for (size_t i = 0; i < tensor_alg_size(Width, Depth); ++i) {
            *(p++) += *(rp++) / sca;
        }
        return *this;
    }

    giles_template_tensor &sub_scal_div(const giles_template_tensor &rhs, double sca) {
        double *p = range_begin();
        const double *rp = rhs.range_begin();
        for (size_t i = 0; i < tensor_alg_size(Width, Depth); ++i) {
            *(p++) -= *(rp++) / sca;
        }
        return *this;
    }

    giles_template_tensor &operator*=(const giles_template_tensor &rhs) {
        giles_template_tensor tmp;
        pass_through op;
        inplace_multiplication_impl(tmp, *this, rhs, size_array, op);
        base_type::swap(tmp);
        return *this;
    }

    giles_template_tensor operator*(const giles_template_tensor &rhs) const {
        giles_template_tensor result;
        pass_through op;
        inplace_multiplication_impl(result, *this, rhs, size_array, op);
        return result;
    }

    giles_template_tensor &mul_scal_mul(const giles_template_tensor &rhs, double sca) {
        post_mul op(sca);
        giles_template_tensor tmp;
        inplace_multiplication_impl(tmp, *this, rhs, size_array, op);
        base_type::swap(tmp);
        return *this;
    }

    giles_template_tensor &mul_scal_div(const giles_template_tensor &rhs, double sca) {
        post_div op(sca);
        giles_template_tensor tmp;
        inplace_multiplication_impl(tmp, *this, rhs, size_array, op);
        base_type::swap(tmp);
        return *this;
    }

    bool operator==(const giles_template_tensor &other) const {
        for (size_t i = 0; i < tensor_alg_size(Width, Depth); ++i) {
            if ((*this)[i] != other[i])
                return false;
        }
        return true;
    }

    bool operator!=(const giles_template_tensor &other) const {
        return !(*this == other);
    }

    friend std::ostream &operator<<(std::ostream &os, const giles_template_tensor &arg) {
        os << "{ ";
        for (auto &v : arg) {
            os << v << ' ';
        }
        return os << '}';
    }

    friend giles_template_tensor exp(const giles_template_tensor &arg) {
        giles_template_tensor result{1.0};
        for (unsigned d = Depth; d >= 1; --d) {
            result.mul_scal_div(arg, static_cast<double>(d));
            result[0] += 1.0;
        }
        return result;
    }

    friend double L1Norm(const giles_template_tensor &arg) noexcept {
        double result = 0.0;
        for (const auto &v : arg) {
            result = std::abs(v);
        }
        return result;
    }

    friend double LInfNorm(const giles_template_tensor &arg) noexcept {
        double result = 0.0;
        for (const auto &v : arg) {
            auto av = std::abs(v);
            result = (av > result) ? av : result;
        }
        return result;
    }
};

template<unsigned Width, unsigned Depth, typename Coeffs>
const typename giles_template_tensor<Width, Depth, Coeffs>::size_type giles_template_tensor<Width, Depth, Coeffs>::tile_width;

template<unsigned Width, unsigned Depth, typename Coeffs>
const typename giles_template_tensor<Width, Depth, Coeffs>::size_type giles_template_tensor<Width, Depth, Coeffs>::tile_size;

}// namespace playground

#endif//TENSOR_PLAYGROUND_GILES_TEMPLATE_TENSOR_H
