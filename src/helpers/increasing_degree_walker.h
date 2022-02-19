//
// Created by sam on 16/02/2022.
//

#ifndef TENSOR_PLAYGROUND_INCREASING_DEGREE_WALKER_H
#define TENSOR_PLAYGROUND_INCREASING_DEGREE_WALKER_H

#include <utility>

namespace playground {

template<unsigned Level, unsigned MaxLevel, template<unsigned> class Fn, bool = (Level <= MaxLevel)>
struct increasing_degree_walker {

    Fn<Level> m_fn;

    template<typename... Args>
    void apply(Args &&...args) const noexcept {
        m_fn(args...);
        increasing_degree_walker<Level + 1, MaxLevel, Fn> next;
        next.apply(std::forward<Args>(args)...);
    }
};

template<unsigned MaxLevel, template<unsigned> class Fn>
struct increasing_degree_walker<MaxLevel, MaxLevel, Fn, true> {

    Fn<MaxLevel> m_fn;

    template<typename... Args>
    void apply(Args &&...args) const noexcept {
        m_fn(args...);
    }
};

template<unsigned Level, unsigned MaxLevel, template<unsigned> class Fn>
struct increasing_degree_walker<Level, MaxLevel, Fn, false> {

    template<typename... Args>
    void apply(Args &&...) const noexcept {}
};

}// namespace playground

#endif//TENSOR_PLAYGROUND_INCREASING_DEGREE_WALKER_H
