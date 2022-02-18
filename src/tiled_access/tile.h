//
// Created by sam on 15/02/2022.
//

#ifndef TENSOR_PLAYGROUND_TWO_SIDED_TILE_H
#define TENSOR_PLAYGROUND_TWO_SIDED_TILE_H

#include "implementation.h"
#include <memory>

namespace playground {

template<typename ElementType, size_type N, bool Inline = (N * sizeof(ElementType) <= PLAYGROUND_MAX_INLINE_TILE_SIZE)>
class tile {
  std::unique_ptr<ElementType[]> data;
  static_assert(Inline == false, "this size tile should have been inlined");

 public:
  using value_type = ElementType;
  using pointer = ElementType *;
  using const_pointer = const ElementType *;

  tile() : data(new ElementType[N]) {
  }

  explicit constexpr operator pointer() noexcept {
    return data.get();
  }

  explicit constexpr operator const_pointer() const noexcept {
    return data.get();
  }

  ElementType &operator[](const size_type i) noexcept {
    return data[i];
  }

  const ElementType &operator[](const size_type i) const noexcept {
    return data[i];
  }
};

template<typename ElementType, size_type N>
class tile<ElementType, N, true> {

  static_assert(N * sizeof(ElementType) <= PLAYGROUND_MAX_INLINE_TILE_SIZE,
                "this size tile should not have been inlined");
  ElementType data[N];

 public:
  using value_type = ElementType;
  using pointer = ElementType *;
  using const_pointer = const ElementType *;
  using reference = ElementType &;
  using const_reference = const ElementType &;

  tile() : data() {
  }

  explicit constexpr operator pointer() noexcept {
    return data;
  }

  explicit constexpr operator const_pointer() const noexcept {
    return data.get();
  }

  reference operator[](const size_type i) noexcept {
    return data[i];
  }

  const_reference operator[](const size_type i) const noexcept {
    return data[i];
  }
};

}// namespace playground

#endif//TENSOR_PLAYGROUND_TWO_SIDED_TILE_H
