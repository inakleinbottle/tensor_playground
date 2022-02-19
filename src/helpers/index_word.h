//
// Created by sam on 17/02/2022.
//

#ifndef TENSOR_PLAYGROUND_INDEX_WORD_H
#define TENSOR_PLAYGROUND_INDEX_WORD_H

#include "tensor_size_info.h"
#include <type_traits>
#include <utility>
#include <iostream>

namespace playground {

enum class index_type {
    TotalOrder,
    WithinDegree
};

template<unsigned Width, typename Integer = std::size_t, index_type Type = index_type::TotalOrder>
class index_word {
    Integer m_data;

   public:
    using integer_type = Integer;

    template<typename I>
    constexpr explicit index_word(I data) : m_data(data) {}

    constexpr index_word(std::initializer_list<integer_type> arg) : m_data() {
        if (Type == index_type::TotalOrder) {
            for (auto a : arg) {
                m_data *= Width;
                m_data += a;
            }
        } else {
            for (auto a : arg) {
                m_data *= Width;
                m_data += a-1;
            }
        }
    }

    template<typename InputIt>
    explicit index_word(InputIt it, InputIt end) : m_data() {
        for (; it != end; ++it) {
            m_data *= Width;
            m_data += *it;
        }
    }

    constexpr index_word &operator=(integer_type data) noexcept {
        m_data = std::move(data);
    }

    template<typename I>
    constexpr index_word &operator=(I &&data) noexcept {
        m_data = std::forward<I>(data);
        return *this;
    }

    constexpr operator integer_type() const noexcept {
        return m_data;
    }

    constexpr operator bool() const noexcept {
        return m_data != 0;
    }

    constexpr std::pair<index_word, index_word> split(unsigned right_length) const noexcept {
        auto split_n = power(Integer(Width), right_length);
        if (Type == index_type::TotalOrder) {
            auto tmp = 1 + (m_data - 1) % split_n;
            return {index_word((m_data - tmp) / split_n), index_word(tmp)};
        } else {
            return {index_word((m_data) / split_n), index_word(m_data % split_n)};
        }
    }

    template<unsigned RightLength>
    constexpr std::pair<index_word, index_word> split() const noexcept {
        static_assert(RightLength > 0, "RightLength must be positive");
        constexpr auto split_n = power(Integer(Width), RightLength);
        if (m_data == 0) {
            return {index_word(0), index_word(0)};
        }
        if (Type == index_type::TotalOrder) {
            auto tmp = 1 + (m_data - 1) % split_n;
            return {index_word((m_data - tmp) / split_n), index_word(tmp)};
        } else {
            return {index_word((m_data) / split_n), index_word(m_data % split_n)};
        }
    }

    index_word reverse() const noexcept {
        index_word tmp(m_data);
        Integer result_inner(0);

        while (tmp) {
            auto tmp2 = tmp.split<1>();
            tmp = tmp2.first;
            result_inner *= Width;
            result_inner += static_cast<Integer>(tmp2.second);
        }
        return index_word{result_inner};
    }

    friend std::ostream &operator<<(std::ostream &os, const index_word &word) noexcept {
        auto l = word.length();

        if (l == 0) {
            return os << "()";
        }
        std::vector<integer_type> letters;
        letters.reserve(l);
        auto tmp = word.m_data;

        while (tmp) {
            auto v = 1 + (tmp - 1) % Width;
            letters.push_back(v);
            tmp -= v;
            tmp /= Width;
        }
        std::reverse(letters.begin(), letters.end());

        os << '(';
        for (auto let : letters) {
            os << let;
            if (--l > 0)
                os << ',';
        }
        os << ')';

        return os;
    }

    constexpr integer_type idx() const noexcept {
        return m_data;
    }

    constexpr unsigned length() const noexcept {
        return (m_data == 0) ? 0 : log<Width>(m_data * (Width - 1) + 1);
    }

    constexpr bool operator<(const integer_type &rhs) const noexcept {
        return m_data < rhs;
    }

    constexpr bool operator<=(const integer_type &rhs) const noexcept {
        return m_data <= rhs;
    }

    constexpr bool operator>(const integer_type &rhs) const noexcept {
        return m_data > rhs;
    }

    constexpr bool operator>=(const integer_type &rhs) const noexcept {
        return m_data >= rhs;
    }

    constexpr index_word concatenate(index_word rhs, unsigned rhs_length) const noexcept {
        return index_word(m_data * power(Width, rhs_length) + rhs.m_data);
    }

    constexpr index_word &operator++() noexcept {
        ++m_data;
        return *this;
    }
};

template<unsigned Width, typename Integer>
Integer operator*(const index_word<Width, Integer> &lhs, const Integer &rhs) noexcept {
    return static_cast<Integer>(lhs) * rhs;
}

template<unsigned Width, typename Integer>
Integer operator*(const Integer &lhs, const index_word<Width, Integer> &rhs) noexcept {
    return lhs * static_cast<Integer>(rhs);
}

template<unsigned Width, typename Integer>
Integer operator+(const index_word<Width, Integer> &lhs, const Integer &rhs) noexcept {
    return static_cast<Integer>(lhs) + rhs;
}

template<unsigned Width, typename Integer>
Integer operator+(const Integer &lhs, const index_word<Width, Integer> &rhs) noexcept {
    return lhs + static_cast<Integer>(rhs);
}

}// namespace playground

#endif//TENSOR_PLAYGROUND_INDEX_WORD_H
