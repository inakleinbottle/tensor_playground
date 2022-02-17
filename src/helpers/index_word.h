//
// Created by sam on 17/02/2022.
//

#ifndef TENSOR_PLAYGROUND_INDEX_WORD_H
#define TENSOR_PLAYGROUND_INDEX_WORD_H

#include <type_traits>
#include <utility>
#include "tensor_size_info.h"

namespace playground {

    enum class index_type {
        TotalOrder,
        WithinDegree
    };


    template <unsigned Width, typename Integer=std::size_t, index_type Type=index_type::TotalOrder>
    class index_word
    {
        Integer m_data;

    public:

        constexpr explicit index_word(Integer data) : m_data(data)
        {}

        constexpr operator Integer() const noexcept
        {
            return m_data;
        }

        template <typename I>
        constexpr index_word(std::initializer_list<I> letters) : m_data(0)
        {
            if (Type == index_type::TotalOrder) {
                for (auto let: letters) {
                    m_data *= Width;
                    m_data += static_cast<Integer>(let);
                }
            } else if (Type == index_type::WithinDegree){
                for (auto let : letters) {
                    m_data *= Width;
                    m_data += static_cast<Integer>(let-1);
                }
            }
        }

        std::pair<index_word, index_word> split(unsigned right_letters) const noexcept
        {
            Integer split_n = power(Width, right_letters);
            if (Type == index_type::TotalOrder) {
                return {index_word(m_data / split_n), index_word(m_data % split_n)};
            } else {
                return {index_word(m_data / split_n), index_word(m_data % split_n)};
            }
        }

        template <unsigned RightLetters>
        std::pair<index_word, index_word> split() const noexcept
        {
            constexpr Integer split_n = power(Width, RightLetters);
            if (Type == index_type::TotalOrder) {
                 return {index_word(m_data / split_n), index_word(m_data % split_n)};
            } else {
                return {index_word(m_data / split_n), index_word(m_data % split_n)};
            }
        }

        index_word reverse() const noexcept {
            index_word tmp(m_data);
            Integer inner = 0;

            while (tmp) {
                auto tmp2 = split<1>();
                tmp = tmp2.first;
                inner *= Width;
                inner += tmp2.second.m_data;
            }

            return index_word(inner);
        }

        index_word& operator++() noexcept
        {
            ++m_data;
            return *this;
        }



    };



}



#endif //TENSOR_PLAYGROUND_INDEX_WORD_H
