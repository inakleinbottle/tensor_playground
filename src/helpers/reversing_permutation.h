//
// Created by sam on 10/12/2021.
//

#ifndef TENSOR_INVERSE_REVERSING_PERMUTATION_H
#define TENSOR_INVERSE_REVERSING_PERMUTATION_H

#include <array>
#include <tensor_size_info.h>

#ifndef TENSOR_INVERSE_MAX_RECURSION_PERMUTE
#define TENSOR_INVERSE_MAX_RECURSION_PERMUTE 2048
#endif

namespace playground {
namespace dtl {

using deg_t = unsigned;

template<deg_t Width>
struct reversing_permutation_rework {

    template<deg_t Level, deg_t ChunkSize = 1>
    struct recurse_tag {
        static constexpr auto factor = power(Width, Level - ChunkSize);
        static constexpr auto truncate = power(Width, ChunkSize);
    };

    template<typename Tag>
    static constexpr size_t prefix(size_t idx, Tag) noexcept {
        return idx / Tag::factor;
    }

    template<typename Tag>
    static constexpr size_t suffix(size_t idx, Tag) noexcept {
        return idx % Tag::truncate;
    }

    template<typename Tag>
    static constexpr size_t middle_word(size_t idx, Tag) noexcept {
        return (idx % Tag::factor) / Tag::truncate;
    }

    template<deg_t Level, deg_t ChunkSize>
    static constexpr size_t
    permute_idx(size_t idx, recurse_tag<Level, ChunkSize> tag) noexcept {
        using tag_t = recurse_tag<Level, ChunkSize>;
        using outer_t = recurse_tag<ChunkSize, 1>;
        using inner_t = recurse_tag<(Level > 2 * ChunkSize) ? Level - 2 * ChunkSize : 0, ChunkSize>;
        return permute_idx(suffix(idx, tag), outer_t()) * tag_t::factor
            + permute_idx(middle_word(idx, tag), inner_t()) * tag_t::truncate
            + permute_idx(prefix(idx, tag), outer_t());
    }

    template<deg_t ChunkSize>
    static constexpr size_t
    permute_idx(size_t idx, recurse_tag<0, ChunkSize> tag) noexcept {
        return 0;
    }

    template<deg_t ChunkSize>
    static constexpr size_t
    permute_idx(size_t idx, recurse_tag<1, ChunkSize> tag) noexcept {
        return idx;
    }
};

template<unsigned Width, unsigned Level>
struct reversing_permutation {
    using size_type = size_t;
    using cycle_type = std::pair<size_type, size_type>;

    //static const std::vector<cycle_type> cycles;
    static constexpr size_type factor = power(Width, Level - 1);

    static constexpr size_type first_letter(size_type idx) {
        return idx / factor;
    }

    static constexpr size_type last_letter(size_type idx) {
        return idx % Width;
    }

    static constexpr size_type middle_word(size_type idx) {
        /*
     * Writing idx = l_2*Width^{Level-1} + index(middle_word)*Width + l1
     * we can rearrange to get index(middle_word) = (idx - l1 - l2*Width^{Level-1})/Width.
     * Although since l1 < Width, we can ignore it since floor division will take care of it.
     * The brackets on the right-hand side can then be realised as idx % Width^{Level-1}
     */
        return (idx % factor) / Width;
    }

    static constexpr size_type permute_idx(size_type idx) {
        static_assert(Level - 2 > 0, "Level must be at least 3 in this specialisation");
        using next = reversing_permutation<Width, Level - 2>;

        constexpr size_type shift = power(Width, Level - 1);
        return last_letter(idx) * shift + next::permute_idx(middle_word(idx)) * Width + first_letter(idx);
    }

   private:
    template<size_type I, size_type N, bool = (permute_idx(I) > I)>
    struct template_tile_permuter {
        template<typename T>
        inline static void eval(T *__restrict tile) noexcept {
            assert(I < power(Width, Level));
            assert(permute_idx(I) < power(Width, Level));
            std::swap(tile[I], tile[permute_idx(I)]);
            template_tile_permuter<(I + 1), N>::eval(tile);
        }

        template<typename T>
        static void eval(std::array<T, N> &tile) noexcept {
            std::swap(tile[I], tile[permute_idx(I)]);
            template_tile_permuter<(I + 1), N>::eval(tile);
        }

        template<typename T, typename Op>
        inline static void eval(const T *__restrict src, T *__restrict dst, Op signer) noexcept {
            assert(I < power(Width, Level));
            assert(permute_idx(I) < power(Width, Level));
            dst[permute_idx(I)] = signer(src[I]);
            template_tile_permuter<(I + 1), N>::eval(src, dst, signer);
        }
    };

    template<size_type I, size_type N>
    struct template_tile_permuter<I, N, false> {
        template<typename T>
        inline static void eval(T *__restrict tile) noexcept {
            template_tile_permuter<(I + 1), N>::eval(tile);
        }

        template<typename T>
        static void eval(std::array<T, N> &tile) noexcept {
            template_tile_permuter<(I + 1), N>::eval(tile);
        }

        template<typename T, typename Op>
        inline static void eval(const T *__restrict src, T *__restrict dst, Op signer) noexcept {
            assert(I < power(Width, Level));
            assert(permute_idx(I) < power(Width, Level));
            dst[permute_idx(I)] = signer(src[I]);
            template_tile_permuter<(I + 1), N>::eval(src, dst, signer);
        }
    };

    template<size_type N>
    struct template_tile_permuter<N, N, false> {
        template<typename T>
        static void eval(T *__restrict tile) noexcept {
        }

        template<typename T>
        static void eval(std::array<T, N> &tile) noexcept {
        }

        template<typename T, typename Op>
        inline static void eval(const T *__restrict src, T *__restrict dst, Op signer) noexcept {
        }
    };

    template<size_type N>
    struct template_tile_permuter<N, N, true> {
        template<typename T>
        static void eval(T *__restrict tile) noexcept {
        }

        template<typename T>
        static void eval(std::array<T, N> &tile) noexcept {
        }

        template<typename T, typename Op>
        inline static void eval(const T *__restrict src, T *__restrict dst, Op signer) noexcept {
        }
    };

    template<size_type S, typename T>
    inline static typename std::enable_if<(S > TENSOR_INVERSE_MAX_RECURSION_PERMUTE)>::type
    impl(T *__restrict tile) noexcept {
        constexpr size_type num = power(Width, Level);

        for (size_type i = 0; i < num; ++i) {
            //const auto j = std::max(i, permute_idx(i));
            //if (j > i) {
            //    std::swap(tile[i], tile[j]);
            //}
            assert(permute_idx(i) < power(Width, Level));
            std::swap(tile[i], tile[std::max(i, permute_idx(i))]);
        }
    }

    template<size_type S, typename T, typename Op>
    inline static typename std::enable_if<(S > TENSOR_INVERSE_MAX_RECURSION_PERMUTE)>::type
    impl(const T *__restrict src, T *__restrict dst, Op signer) noexcept {
        constexpr size_type num = power(Width, Level);

        for (size_type i = 0; i < num; ++i) {
            assert(permute_idx(i) < power(Width, Level));
            dst[permute_idx(i)] = signer(src[i]);
        }
    }

    template<size_type S, typename T>
    inline static typename std::enable_if<(S <= TENSOR_INVERSE_MAX_RECURSION_PERMUTE)>::type
    impl(T *__restrict tile) noexcept {
        template_tile_permuter<1, S>::eval(tile);
    }

    template<typename T, size_type N>
    static typename std::enable_if<(N > TENSOR_INVERSE_MAX_RECURSION_PERMUTE)>::type
    impl(std::array<T, N> &tile) noexcept {
        for (size_type i = 0; i < N; ++i) {
            assert(permute_idx(i) < power(Width, Level));
            std::swap(tile[i], tile[std::max(i, permute_idx(i))]);
        }
    }

    template<size_type N, typename T, typename Op>
    static typename std::enable_if<(N <= TENSOR_INVERSE_MAX_RECURSION_PERMUTE)>::type
    impl(const T *__restrict src, T *__restrict dst, Op signer) noexcept {
        template_tile_permuter<0, N>::eval(src, dst, signer);
    }

   public:
    template<typename T>
    void operator()(T *__restrict tile) const noexcept {
        impl<power(Width, Level)>(tile);
    }

    template<typename T, size_type N>
    void operator()(std::array<T, N> &tile) const noexcept {
        static_assert(N == power(Width, Level), "N must correspond to Width and Level");
        impl(tile);
    }

    template<typename T, typename Op>
    void operator()(const T *__restrict src, T *__restrict dst, Op signer) const noexcept {
        impl<power(Width, Level)>(src, dst, signer);
    }
};

template<unsigned Width>
struct reversing_permutation<Width, 2> {
    using size_type = size_t;
    using cycle_type = std::pair<size_type, size_type>;
    static const unsigned Level = 2;

    static constexpr size_type first_letter(size_type idx) {
        return idx / Width;
    }

    static constexpr size_type last_letter(size_type idx) {
        return idx % Width;
    }

    static constexpr size_type permute_idx(size_type idx) {
        return last_letter(idx) * Width + first_letter(idx);
    }

   private:
    template<size_type I, size_type N, bool = (permute_idx(I) > I)>
    struct template_tile_permuter {
        template<typename T>
        static void eval(T *__restrict tile) noexcept {
            assert(I < power(Width, Level));
            assert(permute_idx(I) < power(Width, Level));
            std::swap(tile[I], tile[permute_idx(I)]);
            //T tmp = std::move(tile[I]);
            //tile[I] = std::move(tile[permute_idx(I)]);
            //tile[permute_idx(I)] = std::move(tmp);
            template_tile_permuter<(I + 1), N>::eval(tile);
        }

        template<typename T>
        static void eval(std::array<T, N> &tile) noexcept {
            assert(I < power(Width, Level));
            assert(permute_idx(I) < power(Width, Level));
            std::swap(tile[I], tile[permute_idx(I)]);
            template_tile_permuter<(I + 1), N>::eval(tile);
        }

        template<typename T, typename Op>
        inline static void eval(const T *__restrict src, T *__restrict dst, Op signer) noexcept {
            assert(I < power(Width, Level));
            assert(permute_idx(I) < power(Width, Level));
            dst[permute_idx(I)] = signer(src[I]);
            template_tile_permuter<(I + 1), N>::eval(src, dst, signer);
        }
    };

    template<size_type I, size_type N>
    struct template_tile_permuter<I, N, false> {
        template<typename T>
        static void eval(T *__restrict tile) noexcept {
            template_tile_permuter<(I + 1), N>::eval(tile);
        }

        template<typename T>
        static void eval(std::array<T, N> &tile) noexcept {
            template_tile_permuter<(I + 1), N>::eval(tile);
        }

        template<typename T, typename Op>
        inline static void eval(const T *__restrict src, T *__restrict dst, Op signer) noexcept {
            assert(I < power(Width, Level));
            assert(permute_idx(I) < power(Width, Level));
            dst[permute_idx(I)] = signer(src[I]);
            template_tile_permuter<(I + 1), N>::eval(src, dst, signer);
        }
    };

    template<size_type N>
    struct template_tile_permuter<N, N, false> {
        template<typename T>
        static void eval(T *__restrict tile) noexcept {
        }

        template<typename T>
        static void eval(std::array<T, N> &tile) noexcept {
        }

        template<typename T, typename Op>
        inline static void eval(const T *__restrict src, T *__restrict dst, Op signer) noexcept {
        }
    };

    template<size_type S, typename T>
    inline static typename std::enable_if<(S > TENSOR_INVERSE_MAX_RECURSION_PERMUTE)>::type
    impl(T *__restrict tile) noexcept {
        constexpr size_type num = power(Width, Level);
        for (size_type i = 0; i < num; ++i) {
            //const auto j = permute_idx(i);
            //if (j > i) {
            //    std::swap(tile[i], tile[j]);
            //}
            assert(permute_idx(i) < power(Width, Level));
            std::swap(tile[i], tile[std::max(i, permute_idx(i))]);
        }
    }

    template<typename T, size_type N>
    static typename std::enable_if<(N > TENSOR_INVERSE_MAX_RECURSION_PERMUTE)>::type
    impl(std::array<T, N> &tile) noexcept {
        for (size_type i = 0; i < N; ++i) {
            //const auto j = permute_idx(i);
            //if (j > i) {
            //    std::swap(tile[i], tile[j]);
            //}
            assert(permute_idx(i) < power(Width, Level));
            std::swap(tile[i], tile[std::max(i, permute_idx(i))]);
        }
    }

    template<size_type S, typename T, typename Op>
    inline static typename std::enable_if<(S > TENSOR_INVERSE_MAX_RECURSION_PERMUTE)>::type
    impl(const T *__restrict src, T *__restrict dst, Op signer) noexcept {
        constexpr size_type num = power(Width, Level);

        for (size_type i = 0; i < num; ++i) {
            assert(permute_idx(i) < power(Width, Level));
            dst[permute_idx(i)] = signer(src[i]);
        }
    }

    template<size_type S, typename T>
    inline static typename std::enable_if<(S <= TENSOR_INVERSE_MAX_RECURSION_PERMUTE)>::type
    impl(T *__restrict tile) noexcept {
        template_tile_permuter<1, S>::eval(tile);
    }

    template<size_type N, typename T>
    inline static typename std::enable_if<(N <= TENSOR_INVERSE_MAX_RECURSION_PERMUTE)>::type
    impl(std::array<T, N> &tile) noexcept {
        template_tile_permuter<1, N>::eval(tile);
    }

    template<size_type N, typename T, typename Op>
    static typename std::enable_if<(N <= TENSOR_INVERSE_MAX_RECURSION_PERMUTE)>::type
    impl(const T *__restrict src, T *__restrict dst, Op signer) noexcept {
        template_tile_permuter<0, N>::eval(src, dst, signer);
    }

   public:
    /// Operate inplace on a single tile
    template<typename T>
    void operator()(T *__restrict tile) const noexcept {
        impl<power(Width, 2)>(tile);
    }

    constexpr size_type operator()(size_type idx) const noexcept {
        return permute_idx(idx);
    }

    template<typename T, typename Op>
    void operator()(const T *__restrict src, T *__restrict dst, Op signer) const noexcept {
        impl<power(Width, Level)>(src, dst, signer);
    }
};

template<unsigned Width>
struct reversing_permutation<Width, 1> {
    using size_type = size_t;
    using cycle_type = std::pair<size_type, size_type>;
    static const unsigned Level = 1;

    static constexpr size_type permute_idx(size_type idx) {
        return idx;
    }

    /// Operate inplace on a single tile
    template<typename T>
    void operator()(T *__restrict tile) const noexcept {
        // Do Nothing!
    }

    template<typename T, size_type N>
    void operator()(std::array<T, N> &tile) const noexcept {}

    constexpr size_type operator()(size_type idx) const noexcept {
        return permute_idx(idx);
    }

    template<typename T, typename Op>
    void operator()(const T *__restrict src, T *__restrict dst, Op signer) const noexcept {
        for (size_type i = 0; i < Width; ++i) {
            dst[i] = signer(src[i]);
        }
    }
};

template<unsigned Width>
struct reversing_permutation<Width, 0> {
    using size_type = size_t;
    using cycle_type = std::pair<size_type, size_type>;
    static const unsigned Level = 0;

    static constexpr size_type permute_idx(size_type idx) {
        return idx;
    }

    /// Operate inplace on a single tile
    template<typename T>
    void operator()(T *__restrict tile) const noexcept {
        // Do nothing!
    }

    template<typename T, size_type N>
    void operator()(std::array<T, N> &tile) const noexcept {}

    constexpr size_type operator()(size_type idx) const noexcept {
        return permute_idx(idx);
    }

    template<typename T, typename Op>
    void operator()(const T *__restrict src, T *__restrict dst, Op signer) const noexcept {
        dst[0] = signer(src[0]);
    }
};

}// namespace dtl
}// namespace playground
#endif//TENSOR_INVERSE_REVERSING_PERMUTATION_H
