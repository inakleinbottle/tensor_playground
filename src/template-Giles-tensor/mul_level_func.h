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
    static typename std::enable_if<(Level >= 2 * tile_letters)>::type
    do_level(pointer out_ptr,
             const GilesTensor &lhs,
             const GilesTensor &rhs,
             index_key middle,
             index_key rmiddle,
             Op &op) noexcept
    {
        constexpr size_type this_size = power(Width, ThisLevelLetters);
        using permute = reversing_permutation<Width, ThisLevelLetters>;

        tile_type this_tile;

        for (size_type i=0; i<this_size; ++i) {
            for (size_type j=0; j<this_size; ++j) {
                auto ri = permute::permute_idx(i);
                auto rj = permute::permute_idx(j);
                next::template do_level<Level>();
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
            next::template do_level<Level>(this_tile, )

            for (size_type i = 0; i < tile_width; ++i) {
                for (size_type j = 0; j < tile_width; ++j) {
                    out_ptr[i * stride + j] += this_tile[i * tile_width];
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
        Op &op) noexcept
    {
        next::template do_level<Level>(out_ptr, lhs, rhs, op);
    }
};

template<typename Coeffs, size_type Width, unsigned TileLetters>
struct multiplication_level_helper<Coeffs, Width, TileLetters> {
    using degree_type = unsigned;
    static constexpr degree_type tile_letters = TileLetters;
    static constexpr size_type tile_width = power(Width, TileLetters);
    static constexpr size_type tile_size = tile_width * tile_width;

    using pointer = Coeffs *;
    using const_pointer = const Coeffs *;
    using index_key = index_word<Width, size_type>;
    using tile_type = tile<Coeffs, tile_size>;

    template<degree_type ThisLevel>
    struct left_side_cases {
        static_assert(ThisLevel < tile_letters, "ThisLevel must be less than tile_letters");

        template<typename Op>
        void inner(
            tile_type &this_tile,
            const_pointer lhs_ptr,
            const_pointer rhs_ptr,
            Op &op,
            size_type i1) const noexcept {
            constexpr size_type loop_bound = power(Width, ThisLevel);

            for (size_type i2 = 0; i2 < loop_bound; ++i2) {
                auto i = i1 * loop_bound + i2;
                for (size_type j = 0; j < tile_width; ++j) {
                    auto ri1 = dtl::reversing_permutation<Width, ThisLevel>::permute_idx(i1);
                    this_tile[i * tile_width + j] += op(lhs_ptr[ri1] * rhs_ptr[i2 * tile_width + j]);
                }
            }
        }

        template<typename Op>
        void operator()(
            tile_type &this_tile,
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
            tile_type &this_tile,
            const_pointer lhs_ptr,
            const_pointer rhs_ptr,
            Op &op,
            size_type i,
            size_type j1) const noexcept {
            constexpr size_type loop_bound = power(Width, tile_letters - ThisLevel);
            constexpr size_type shift = power(Width, ThisLevel);

            for (size_type j2 = 0; j2 < loop_bound; ++j2) {
                auto j = j1 * shift + j2;
                auto jr = dtl::reversing_permutation<Width, ThisLevel>::permute_idx(j2);
                this_tile[i * tile_width + j] += op(lhs_ptr[jr * tile_width + i] * rhs_ptr[j1]);
            }
        }

        template<typename Op>
        void operator()(
            tile_type &this_tile,
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
            tile_type &this_tile,
            const_pointer lhs_ptr,
            const_pointer rhs_ptr,
            index_key middle_word,
            Op &op) const noexcept {
            auto split = middle_word.template split<ThisLevel>();
            auto lhs_middle = split.first.reverse();
            auto rhs_middle = split.second;

            for (size_type i = 0; i < tile_width; ++i) {
                for (size_type j = 0; j < tile_width; ++j) {
                    auto lhs_idx = lhs_middle * tile_width + i;
                    auto rhs_idx = rhs_middle * tile_width + j;
                    this_tile[i * tile_width + j] += op(lhs_ptr[lhs_idx] * rhs_ptr[rhs_idx]);
                }
            }
        }
    };

    template<unsigned Level, typename Op>
    static void
    do_level_tiled(
        pointer out_ptr,
        const_pointer lhs,
        const_pointer lhs_r,
        const_pointer rhs,
        index_key middle_idx,
        index_key rmiddle_idx,
        Op &op) noexcept {

        // Handle cases of 0*out_depth and out_depth*0
        {
            auto lhs_unit = lhs[0], rhs_unit = rhs[0];
            const auto *lhs_ptr = lhs + static_cast<size_type>(middle_idx) * tile_width;
            const auto *rhs_ptr = rhs + static_cast<size_type>(middle_idx) * tile_width;
            static constexpr auto stride = power(Width, Level - tile_letters);

            /// First do the zero*Depth and Depth*zero case
            for (size_type i = 0; i < tile_width; ++i) {
                for (size_type j = 0; j < tile_width; ++j) {
                    out_ptr[i * tile_width + j] = op(lhs_unit * rhs_ptr[i * stride + j]);
                }
            }
            for (size_type i = 0; i < tile_width; ++i) {
                for (size_type j = 0; j < tile_width; ++j) {
                    out_ptr[i * tile_width + j] += op(lhs_ptr[i * stride + j] * rhs_unit);
                }
            }


        }

        increasing_degree_walker<1, TileLetters - 1, left_side_cases> left_walker;
        increasing_degree_walker<1, TileLetters - 1, right_side_cases> right_walker;
        increasing_degree_walker<1, Level - 2 * TileLetters - 1, middle_cases> middle_walker;

        left_walker.apply(
            out_ptr,
            lhs_r,
            rhs + middle_idx * tile_width,
            op);
        middle_walker.apply(
            out_ptr,
            lhs_r,
            rhs,
            middle_idx,
            op);
        right_walker.apply(
            out_ptr,
            lhs_r + rmiddle_idx * tile_width,
            rhs,
            op);
    }

    template <unsigned Level, typename GilesTensor, typename Op>
    static typename std::enable_if<(Level >= 2 * tile_letters)>::type
    do_level(pointer out_ptr,
             const GilesTensor &lhs,
             const GilesTensor &rhs,
             index_key middle,
             index_key rmiddle,
             Op &op) noexcept
    {
        do_level_tiled<Level>(
            out_ptr,
            lhs.range_begin(),
            lhs.reverse_data.data(),
            rhs.range_begin(),
            middle,
            rmiddle,
            op);
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
            do_level<Level>(this_tile, lhs, rhs, mid, rmid, op);

            for (size_type i=0; i<tile_width; ++i) {
                for (size_type j=0; j<tile_width; ++j) {
                    out_ptr[i*stride + j] += this_tile[i*tile_width];
                }
            }
        }

        do_level<Level - 1>(out_ptr, lhs, rhs, op);
    }

    template<degree_type ThisLevel>
    struct untiled_outer_helper {

        template<degree_type LhsLevel>
        struct untiled_inner_helper {
            template<typename Op>
            void operator()(pointer out_ptr, const_pointer lhs_ptr, const_pointer rhs_ptr, Op &op) const noexcept {
                static constexpr degree_type rhs_level = ThisLevel - LhsLevel;
                static constexpr size_type lhs_sod = tensor_start_of_degree(Width, LhsLevel);
                static constexpr size_type rhs_sod = tensor_start_of_degree(Width, rhs_level);

                auto *lhs_p = lhs_ptr + lhs_sod;
                auto *rhs_p = rhs_ptr + rhs_sod;

                for (std::ptrdiff_t i = 0; i < power(Width, LhsLevel); ++i) {
                    for (std::ptrdiff_t j = 0; j < power(Width, rhs_level); ++j) {
                        *(out_ptr++) += op(lhs_p[i] * rhs_p[j]);
                    }
                }
            }
        };
    };

    template<unsigned Level, typename GilesTensor, typename Op>
    static typename std::enable_if<(Level > 0 && Level < 2 * tile_letters)>::type
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
    static typename std::enable_if<(Level == 0 && Level < 2 * tile_letters)>::type
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

}// namespace dtl
}// namespace playground

#endif//TENSOR_PLAYGROUND_MUL_LEVEL_FUNC_H
