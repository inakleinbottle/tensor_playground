//
// Created by sam on 16/02/2022.
//

#ifndef TENSOR_PLAYGROUND_DECREASING_DEGREE_WALKER_H
#define TENSOR_PLAYGROUND_DECREASING_DEGREE_WALKER_H

#include <utility>

namespace playground {

template<unsigned Level, unsigned MinLevel, template<unsigned> class Fn, bool = (Level > MinLevel)>
struct decreasing_degree_walker {

  Fn<Level> m_fn;

  template<typename... Args>
  void apply(Args &&...args) noexcept {
    m_fn(args...);
    increasing_degree_walker<Level - 1, MinLevel, Fn> next;
    next.apply(std::forward<Args>(args)...);
  }
};

template<unsigned MinLevel, template<unsigned> class Fn>
struct decreasing_degree_walker<MinLevel, MinLevel, Fn, false> {

  Fn<MinLevel> m_fn;

  template<typename... Args>
  void apply(Args &&...args) noexcept {
    m_fn(args...);
  }
};

template<unsigned Level, unsigned MinLevel, template<unsigned> class Fn>
struct decreasing_degree_walker<Level, MinLevel, Fn, false> {

  template<typename... Args>
  void apply(Args &&...) noexcept {}
};

}// namespace playground

#endif//TENSOR_PLAYGROUND_DECREASING_DEGREE_WALKER_H
