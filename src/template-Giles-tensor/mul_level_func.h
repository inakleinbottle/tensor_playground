//
// Created by sam on 15/02/2022.
//

#ifndef TENSOR_PLAYGROUND_MUL_LEVEL_FUNC_H
#define TENSOR_PLAYGROUND_MUL_LEVEL_FUNC_H

#include "implementation.h"
#include <decreasing_degree_walker.h>
#include <increasing_degree_walker.h>
#include <index_word.h>
#include <reversing_permutation.h>
#include <strided_rw.h>
#include <tile.h>

#include <type_traits>
#include <utility>

namespace playground {
namespace dtl {

template<unsigned... Letters>
struct total_letters;

template<unsigned FirstLetters, unsigned... Letters>
struct total_letters<FirstLetters, Letters...> {
    static constexpr size_type value = FirstLetters + total_letters<Letters...>::value;
};

template<unsigned Letter>
struct total_letters<Letter> {
    static constexpr size_type value = Letter;
};

template<>
struct total_letters<> {
    static constexpr size_type value = 0;
};

template<typename Coeffs, size_type Width, unsigned... TileLetters>
struct multiplication_level_helper;

template<typename Coeffs, size_type Width, unsigned ThisLevelLetters, unsigned... RemainingLetters>
struct multiplication_level_helper<Coeffs, Width, ThisLevelLetters, RemainingLetters...> {
    static_assert(ThisLevelLetters > 0, "tile size must be strictly positive");
    static constexpr size_type tile_width = power(Width, ThisLevelLetters);
    static constexpr size_type tile_size = tile_width * tile_width;
    static constexpr size_type tile_letters = total_letters<ThisLevelLetters, RemainingLetters...>::value;
    static constexpr size_type remaining_letters = total_letters<RemainingLetters...>::value;

    using pointer = Coeffs *;
    using const_pointer = const Coeffs *;
    using index_key = index_word<Width, size_type>;
    using tile_type = tile<Coeffs, tile_size>;
    using next = multiplication_level_helper<Coeffs, Width, RemainingLetters...>;

    template<unsigned Level, typename GilesTensor, typename Op>
    static void
    do_tiled_level(pointer out_ptr,
                   const GilesTensor &lhs,
                   const GilesTensor &rhs,
                   index_key middle,
                   index_key rmiddle,
                   Op &op) noexcept {
        constexpr auto inner_word_length = Level - 2 * tile_letters;
        constexpr size_type this_width = power(Width, ThisLevelLetters);
        constexpr size_type this_size = this_width * this_width;
        constexpr size_type level_shift = power(Width, inner_word_length + ThisLevelLetters);
        using permute = reversing_permutation<Width, ThisLevelLetters>;

        typename next::tile_type this_tile;

        for (size_type i = 0; i < this_width; ++i) {
            for (size_type j = 0; j < this_width; ++j) {
                auto ri = permute::permute_idx(i);
                auto rj = permute::permute_idx(j);
                next::template do_tiled_level<2 * tile_letters>(
                    static_cast<pointer>(this_tile),
                    lhs, rhs,
                    index_key(i * level_shift + middle + j),
                    index_key(rj * level_shift + rmiddle + ri),
                    op);

                auto offset = (i * this_width + j);
                stride_write<double, next::tile_width, tile_width, tile_width>(out_ptr + offset, static_cast<const_pointer>(this_tile));
            }
        }
    }

    template<unsigned Level, typename GilesTensor, typename Op>
    static typename std::enable_if<(Level >= 2 * tile_letters)>::type
    do_level(
        pointer out_ptr,
        const GilesTensor &lhs,
        const GilesTensor &rhs,
        Op &op) noexcept {
        static constexpr auto mid_deg = Level - 2 * tile_letters;
        static constexpr auto stride = power(Width, mid_deg);
        index_key mid(tensor_start_of_degree(Width, mid_deg));
        size_type max_idx(tensor_start_of_degree(Width, mid_deg + 1));

        tile_type this_tile;
        for (; mid < max_idx; ++mid) {
            index_key rmid = mid.reverse();
            do_tiled_level<Level>(static_cast<pointer>(this_tile), lhs, rhs, mid, rmid, op);

            for (size_type i = 0; i < tile_width; ++i) {
                for (size_type j = 0; j < tile_width; ++j) {
                    out_ptr[i * stride + j] += this_tile[i * tile_width + j];
                }
            }
        }

        //        next::template do_level<Level>(out_ptr, lhs, rhs, op);
        do_level<Level - 1>(out_ptr, lhs, rhs, op);
    }

    template<unsigned Level, typename GilesTensor, typename Op>
    static typename std::enable_if<(Level < 2 * tile_letters)>::type
    do_level(
        pointer out_ptr,
        const GilesTensor &lhs,
        const GilesTensor &rhs,
        Op &op) noexcept {
        next::template do_level<Level>(out_ptr, lhs, rhs, op);
    }
};

template<typename Coeffs, size_type Width, unsigned TileLetters>
struct multiplication_level_helper<Coeffs, Width, TileLetters> {
    static_assert(TileLetters > 0, "Tile cannot be zero sized");
    using degree_type = unsigned;
    static constexpr degree_type tile_letters = TileLetters;
    static constexpr size_type tile_width = power(Width, TileLetters);
    static constexpr size_type tile_size = tile_width * tile_width;

    using pointer = Coeffs *;
    using const_pointer = const Coeffs *;
    using index_key = index_word<Width, size_type>;
    using tile_type = tile<Coeffs, tile_size>;

    template<degree_type OutLevel>
    struct level_impls {

        template<degree_type ThisLevel>
        struct left_side_cases {
            static_assert(ThisLevel < tile_letters, "ThisLevel must be less than tile_letters");

            template<typename Op>
            void inner(
                pointer this_tile,
                const_pointer lhs_ptr,
                const_pointer rhs_ptr,
                Op &op,
                size_type i1) const noexcept {
                constexpr size_type loop_bound = power(Width, ThisLevel);
                constexpr size_type lhs_sod = tensor_start_of_degree(Width, ThisLevel);
                constexpr size_type rhs_sod = tensor_start_of_degree(Width, OutLevel - ThisLevel);

                for (size_type i2 = 0; i2 < loop_bound; ++i2) {
                    auto ri1 = dtl::reversing_permutation<Width, ThisLevel>::permute_idx(i1);
                    auto i = ri1 * loop_bound + i2;
                    for (size_type j = 0; j < tile_width; ++j) {
                        this_tile[i * tile_width + j] += op(lhs_ptr[i1 + lhs_sod] * rhs_ptr[i2 * tile_width + j + rhs_sod]);
                    }
                }
            }

            template<typename Op>
            void operator()(
                pointer this_tile,
                const_pointer lhs_ptr,
                const_pointer rhs_ptr,
                Op &op) const noexcept {
                constexpr size_type loop_bound = power(Width, tile_letters - ThisLevel);
                for (size_type i = 0; i < loop_bound; ++i) {
                    inner(this_tile, lhs_ptr, rhs_ptr, op, i);
                }
            }
        };

        template<degree_type ThisLevel>
        struct right_side_cases {
            static_assert(ThisLevel < tile_letters, "ThisLevel must be less than tile_letters");

            template<typename Op>
            void inner(
                pointer this_tile,
                const_pointer lhs_ptr,
                const_pointer rhs_ptr,
                Op &op,
                size_type i,
                size_type j1) const noexcept {
                constexpr size_type loop_bound = power(Width, tile_letters - ThisLevel);
                constexpr size_type shift = power(Width, ThisLevel);
                constexpr size_type lhs_sod = tensor_start_of_degree(Width, OutLevel - ThisLevel);
                constexpr size_type rhs_sod = tensor_start_of_degree(Width, ThisLevel);
                constexpr auto ts = tile_width;

                auto ir = dtl::reversing_permutation<Width, ThisLevel>::permute_idx(i);
                for (size_type j2 = 0; j2 < loop_bound; ++j2) {
                    auto j = j1 * shift + j2;
                    auto jr = dtl::reversing_permutation<Width, ThisLevel>::permute_idx(j2);

                    auto v1 = lhs_ptr[jr * tile_width + ir + lhs_sod];
                    auto v2 = rhs_ptr[j1 + rhs_sod];
                    auto tmp = op(v1 * v2);
                    this_tile[i * tile_width + j] += tmp;
                }
            }

            template<typename Op>
            void operator()(
                pointer this_tile,
                const_pointer lhs_ptr,
                const_pointer rhs_ptr,
                Op &op) const noexcept {
                constexpr size_type loop_bound = power(Width, ThisLevel);
                for (size_type i = 0; i < tile_width; ++i) {
                    for (size_type j = 0; j < loop_bound; ++j) {
                        inner(this_tile, lhs_ptr, rhs_ptr, op, i, j);
                    }
                }
            }
        };

        template<degree_type ThisLevel>
        struct middle_cases {
            template<typename Op>
            void operator()(
                pointer this_tile,
                const_pointer lhs_ptr,
                const_pointer rhs_ptr,
                size_type middle_word,
                Op &op) const noexcept {
                //            auto split = middle_word.template split<ThisLevel>();
                //            auto lhs_middle = split.first.reverse();
                //            auto rhs_middle = split.second;
                constexpr size_type split_n = power(Width, ThisLevel);

                auto lhs_middle = middle_word / split_n;
                auto rhs_middle = middle_word % split_n;

                for (size_type i = 0; i < tile_width; ++i) {
                    for (size_type j = 0; j < tile_width; ++j) {
                        auto lhs_idx = lhs_middle + i;
                        auto rhs_idx = rhs_middle + j;
                        this_tile[i * tile_width + j] += op(lhs_ptr[lhs_idx + tensor_start_of_degree(Width, OutLevel - ThisLevel)] * rhs_ptr[rhs_idx + tensor_start_of_degree(Width, ThisLevel)]);
                    }
                }
            }
        };
    };

    template<unsigned Level, typename GilesTensor, typename Op>
    static void
    do_tiled_level(pointer out_ptr,
                   const GilesTensor &lhs,
                   const GilesTensor &rhs,
                   size_type middle_idx,
                   size_type rmiddle_idx,
                   Op &op) noexcept {

        // Handle cases of 0*out_depth and out_depth*0
        {

            auto lhs_unit = lhs[0], rhs_unit = rhs[0];
            constexpr auto deg_start = tensor_start_of_degree(Width, Level);
            const auto *lhs_ptr = lhs.range_begin() + middle_idx * tile_width + deg_start;
            const auto *rhs_ptr = rhs.range_begin() + middle_idx * tile_width + deg_start;
            //            std::cout << middle_idx*tile_width + deg_start << ' ' << lhs_unit << ' ' << rhs_unit << '\n';
            constexpr auto stride = power(Width, Level - tile_letters);

            /// First do the zero*Depth and Depth*zero case
            for (size_type i = 0; i < tile_width; ++i) {

                for (size_type j = 0; j < tile_width; ++j) {
                    //                    std::cout << rhs_ptr[i*stride+j] << ' ';
                    out_ptr[i * tile_width + j] = op(lhs_unit * rhs_ptr[i * stride + j]);
                }
            }
            //            std::cout << '\n';

            for (size_type i = 0; i < tile_width; ++i) {
                for (size_type j = 0; j < tile_width; ++j) {
                    out_ptr[i * tile_width + j] += op(lhs_ptr[i * stride + j] * rhs_unit);
                }
            }
        }

        constexpr degree_type middle_word_max = (Level > 2 * TileLetters) ? Level - 2 * TileLetters - 1 : 0;

        increasing_degree_walker<1, TileLetters - 1, level_impls<Level>::template left_side_cases> left_walker;
        increasing_degree_walker<1, TileLetters - 1, level_impls<Level>::template right_side_cases> right_walker;
        increasing_degree_walker<1, middle_word_max, level_impls<Level>::template middle_cases> middle_walker;

        left_walker.apply(
            out_ptr,
            lhs.reverse_data.data(),
            rhs.range_begin() + middle_idx * tile_width,
            op);
        middle_walker.apply(
            out_ptr,
            lhs.reverse_data.data(),
            rhs.range_begin(),
            middle_idx,
            op);
        right_walker.apply(
            out_ptr,
            lhs.reverse_data.data() + rmiddle_idx * tile_width,
            rhs.range_begin(),
            op);
    }

    template<unsigned Level, typename GilesTensor, typename Op>
    static typename std::enable_if<(Level >= 2 * tile_letters)>::type
    do_level(
        pointer out_ptr,
        const GilesTensor &lhs,
        const GilesTensor &rhs,
        Op &op) noexcept {
        constexpr auto mid_deg = Level - 2 * tile_letters;
        constexpr auto stride = power(Width, Level - tile_letters);
        //        index_key mid(tensor_start_of_degree(Width, mid_deg));
        size_type max_idx(tensor_start_of_degree(Width, mid_deg + 1));
        constexpr auto out_deg_start = tensor_start_of_degree(Width, Level);
        constexpr auto tile_stride = tile_width;

        static_assert(((tile_width - 1) * power(Width, Level - tile_letters) + (tensor_start_of_degree(Width, mid_deg + 1) - 1) + (tile_width - 1)) < tensor_alg_size(Width, GilesTensor::max_depth),
                      "Index access out of bounds");

        tile_type this_tile;
        for (size_type k = 0; k < power(Width, mid_deg); ++k) {
            auto rmid = reversing_permutation<Width, mid_deg>::permute_idx(k);

            do_tiled_level<Level>(static_cast<pointer>(this_tile), lhs, rhs,k,rmid,op);

            for (size_type i = 0; i < tile_width; ++i) {
                for (size_type j = 0; j < tile_width; ++j) {
                    out_ptr[i * stride + k * tile_width + j + out_deg_start] += this_tile[i * tile_width + j];
                }
            }
        }

        do_level<Level - 1>(out_ptr, lhs, rhs, op);
    }

    template<unsigned Level, typename GilesTensor, typename Op>
    static typename std::enable_if<(Level < 2 * tile_letters)>::type
    do_level(
        pointer out_ptr,
        const GilesTensor &lhs,
        const GilesTensor &rhs,
        Op &op) noexcept {
        multiplication_level_helper<Coeffs, Width>::template do_level<Level>(out_ptr, lhs, rhs, op);
    }
};

template<typename Coeffs, size_type Width>
struct multiplication_level_helper<Coeffs, Width> {

    using degree_type = unsigned;
    using pointer = Coeffs *;
    using const_pointer = const Coeffs *;

    static constexpr size_type tile_size = 0;
    static constexpr degree_type tile_letters = 0;
    static constexpr size_type tile_width = 0;

    template<degree_type ThisLevel>
    struct untiled_outer_helper {

        template<degree_type LhsLevel>
        struct untiled_inner_helper {
            template<typename Op>
            void operator()(pointer out_ptr, const_pointer lhs_ptr, const_pointer rhs_ptr, Op &op) const noexcept {
                constexpr degree_type rhs_level = ThisLevel - LhsLevel;
                constexpr size_type lhs_sod = tensor_start_of_degree(Width, LhsLevel);
                constexpr size_type rhs_sod = tensor_start_of_degree(Width, rhs_level);
                constexpr size_type out_sod = tensor_start_of_degree(Width, ThisLevel);

                auto *lhs_p = lhs_ptr + lhs_sod;
                auto *rhs_p = rhs_ptr + rhs_sod;
                out_ptr += out_sod;

                for (std::ptrdiff_t i = 0; i < power(Width, LhsLevel); ++i) {
                    for (std::ptrdiff_t j = 0; j < power(Width, rhs_level); ++j) {
                        *(out_ptr++) += op(lhs_p[i] * rhs_p[j]);
                    }
                }
            }
        };
    };

    template<unsigned Level, typename GilesTensor, typename Op>
    static typename std::enable_if<(Level > 0)>::type
    do_level(
        pointer out_ptr,
        const GilesTensor &lhs,
        const GilesTensor &rhs,
        Op &op) noexcept {
        using ohelper = untiled_outer_helper<Level>;
        decreasing_degree_walker<Level, 0, ohelper::template untiled_inner_helper> walker;
        walker.template apply(out_ptr, lhs.range_begin(), rhs.range_begin(), op);
        do_level<Level - 1>(out_ptr, lhs, rhs, op);
    }

    template<unsigned Level, typename GilesTensor, typename Op>
    static typename std::enable_if<(Level == 0)>::type
    do_level(
        pointer out_ptr,
        const GilesTensor &lhs,
        const GilesTensor &rhs,
        Op &op) noexcept {
        using ohelper = untiled_outer_helper<Level>;
        decreasing_degree_walker<Level, 0, ohelper::template untiled_inner_helper> walker;
        walker.template apply(out_ptr, lhs.range_begin(), rhs.range_begin(), op);
    }
};

template<typename Coeffs, size_type Width>
constexpr size_type multiplication_level_helper<Coeffs, Width>::tile_size;

template<typename Coeffs, size_type Width>
constexpr size_type multiplication_level_helper<Coeffs, Width>::tile_width;

template<typename Coeffs, size_type Width>
constexpr size_type multiplication_level_helper<Coeffs, Width>::tile_size;

template<typename Coeffs, size_type Width>
constexpr typename multiplication_level_helper<Coeffs, Width>::degree_type multiplication_level_helper<Coeffs, Width>::tile_letters;

}// namespace dtl
}// namespace playground

#endif//TENSOR_PLAYGROUND_MUL_LEVEL_FUNC_H
