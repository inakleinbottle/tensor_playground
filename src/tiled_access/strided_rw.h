//
// Created by sam on 15/02/2022.
//

#ifndef TENSOR_PLAYGROUND_STRIDED_RW_H
#define TENSOR_PLAYGROUND_STRIDED_RW_H

#include "implementation.h"

namespace playground {

/**
 * @brief Read data from strided source into contiguous dst
 * @tparam T Element type
 * @tparam Stride Size of the stride
 * @tparam NRows Number of contiguous rows
 * @tparam NCols Number of contiguous columns per row
 * @param dst Pointer to contiguous destination
 * @param src Const pointer to strided source
 */
template<typename T, size_type Stride, size_type NRows, size_type NCols>
void stride_read(T *dst, const T *src) noexcept {
  for (size_type i = 0; i < NRows; ++i) {
    const auto *row_ptr = src + i * Stride;
    for (size_type j = 0; j < NCols; ++j) {
      *(dst++) = *(row_ptr++);
    }
  }
}

/**
 * @brief Write from contigous source to strided destination
 * @tparam T Element type
 * @tparam Stride Size of the stride
 * @tparam NRows Number of contiguous rows
 * @tparam NCols Number of contiguous columns per row
 * @param dst Pointer to strided destination
 * @param src Const pointer to contiguous source
 */
template<typename T, size_type Stride, size_type NRows, size_type NCols>
void stride_write(T *dst, const T *src) noexcept {
  for (size_type i = 0; i < NRows; ++i) {
    auto *row_ptr = dst + i * Stride;
    for (size_type j = 0; j < NCols; ++j) {
      *(row_ptr++) = *(src++);
    }
  }
}

}// namespace playground

#endif//TENSOR_PLAYGROUND_STRIDED_RW_H
