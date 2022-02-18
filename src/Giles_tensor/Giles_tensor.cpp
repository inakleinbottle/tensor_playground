//
// Created by sam on 07/11/2021.
//

#include "Giles_tensor.h"
#include <cassert>

#include <multiplication_helpers.h>
#include <tensor_size_info.h>

using namespace playground;
namespace {

using size_type = typename Giles_tensor::size_type;
using degree_type = typename Giles_tensor::degree_type;

template<typename BinOp, typename UOp>
void do_inplace_binary_operation(
    std::vector<double> &lhs,
    size_type rhs_size,
    const double *prhs,
    BinOp bin_op,
    UOp u_op) {
  auto min_dim = std::min(lhs.size(), rhs_size);

  auto *plhs = lhs.data();
  for (size_type i = 0; i < min_dim; ++i) {
    bin_op(*(plhs++), *(prhs++));
  }

  if (rhs_size > min_dim) {
    lhs.reserve(rhs_size);
    for (size_type i = min_dim; i < rhs_size; ++i) {
      lhs.emplace_back(u_op(*(prhs++)));
    }
  }
}

template<typename BinOp, typename LhsUOp, typename RhsUOp>
void do_binary_operation(
    std::vector<double> &result,
    size_type lhs_size,
    size_type rhs_size,
    const double *plhs,
    const double *prhs,
    BinOp bin_op,
    LhsUOp lhs_uop,
    RhsUOp rhs_uop) noexcept {
  auto min_dim = std::min(lhs_size, rhs_size);

  for (size_type i = 0; i < min_dim; ++i) {
    result.emplace_back(bin_op(*(plhs++), *(prhs++)));
  }

  for (size_type i = min_dim; i < lhs_size; ++i) {
    result.emplace_back(lhs_uop(*(plhs++)));
  }

  for (size_type i = min_dim; i < rhs_size; ++i) {
    result.emplace_back(rhs_uop(*(prhs++)));
  }
}

template<typename Op>
inline void
inplace_multiplication_impl(
    std::vector<double> &out,
    std::vector<double> &result_r,
    unsigned width,
    unsigned lhs_depth,
    unsigned rhs_depth,
    const Giles_tensor &lhs,
    const Giles_tensor &rhs,
    const std::vector<size_type> &start_of_degree_array,
    Op op) noexcept {
  int rhs_deg;

  auto depth = std::max(lhs_depth, rhs_depth);

  double *out_ptr;
  const double *lhs_p, *rhs_p;
  for (int out_deg = static_cast<int>(depth); out_deg >= 0; --out_deg) {
    auto lhs_min_deg = std::max(0, out_deg - static_cast<int>(rhs_depth));
    auto lhs_max_deg = std::min(static_cast<int>(lhs_depth), out_deg);
    for (int lhs_deg = lhs_max_deg; lhs_deg >= lhs_min_deg; --lhs_deg) {
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
}

size_type reverse_index(size_type index, size_type width, degree_type d) noexcept {
  if (d <= 1) return index;

  auto lfact = power(width, d - 1);
  auto rfact = width;

  auto a = index / lfact;
  auto b = (index % lfact) / rfact;
  auto c = index % rfact;

  return c * lfact + reverse_index(b, width, d - 1) * rfact + a;
}

void reverse_data(
    std::vector<double> &reverse_data,
    const std::vector<double> &fwd_data,
    const std::vector<size_type> &start_of_degree_array,
    size_type width,
    degree_type deg) {
  assert(deg < start_of_degree_array.size());

  if (deg == 0)
    return;

  reverse_data.resize(start_of_degree_array[deg]);
  reverse_data[0] = fwd_data[0];

  if (deg == 1)
    return;

  for (size_type i = 1; i <= start_of_degree_array[2]; ++i) {
    reverse_data[i] = fwd_data[i];
  }

  for (degree_type d = 2; d < deg; ++d) {
    for (size_type i = start_of_degree_array[d]; i < start_of_degree_array[d + 1]; ++i) {
      auto idx = i - start_of_degree_array[d];
      reverse_data[reverse_index(idx, width, d) + start_of_degree_array[d]] = fwd_data[i];
    }
  }
}

}// namespace

Giles_tensor::Giles_tensor(degree_type width)
    : m_width(width), m_degree(0), base_type(), size_array(), transposed_data() {
}

void Giles_tensor::populate_size_array(degree_type maxd) {
  size_array.reserve(maxd);
  if (size_array.empty())
    size_array.push_back(0);

  auto sz = size_array.size();

  for (auto i = sz - 1; i < maxd; ++i) {
    size_array.emplace_back(size_array[i] + power(m_width, i));
  }
}

Giles_tensor::Giles_tensor(degree_type width, degree_type depth)
    : m_width(width),
      m_degree(depth),
      base_type(tensor_alg_size(width, depth)),
      size_array(),
      transposed_data(width, (depth > 0) ? depth - 1 : 0) {
  populate_size_array(depth + 2);
}

Giles_tensor::Giles_tensor(degree_type width, degree_type depth, std::initializer_list<double> args)
    : m_width(width), m_degree(depth), base_type(args), size_array(), transposed_data() {
  resize(tensor_alg_size(width, depth));
  transposed_data.resize(tensor_alg_size(width, (depth > 0) ? depth - 1 : 0));
  populate_size_array(depth + 2);
  reverse_data(transposed_data, *this, size_array, width, depth);
}

double *Giles_tensor::range_begin() {
  return data();
}
double *Giles_tensor::range_end() {
  return data() + size();
}
const double *Giles_tensor::range_begin() const {
  return data();
}
const double *Giles_tensor::range_end() const {
  return data() + size();
}

const double *Giles_tensor::trange_begin() const {
  return transposed_data.data();
}
const double *Giles_tensor::trange_end() const {
  return transposed_data.data() + transposed_data.size();
}

degree_type Giles_tensor::degree() const {
  return m_degree;
}
degree_type Giles_tensor::width() const {
  return m_width;
}
Giles_tensor &Giles_tensor::operator+=(const Giles_tensor &rhs) {
  return *this;
}
Giles_tensor &Giles_tensor::operator-=(const Giles_tensor &rhs) {
  return *this;
}
Giles_tensor &Giles_tensor::operator*=(double rhs) {
  return *this;
}
Giles_tensor &Giles_tensor::operator/=(double rhs) {
  return *this;
}
Giles_tensor Giles_tensor::operator+(const Giles_tensor &rhs) const {
  return Giles_tensor(m_width);
}
Giles_tensor Giles_tensor::operator-(const Giles_tensor &rhs) const {
  return Giles_tensor(m_width);
}
Giles_tensor Giles_tensor::operator*(double rhs) const {
  return Giles_tensor(m_width);
}
Giles_tensor Giles_tensor::operator/(double rhs) const {
  return Giles_tensor(m_width);
}
Giles_tensor &Giles_tensor::add_scal_prod(const Giles_tensor &rhs, double sca) {
  return *this;
}
Giles_tensor &Giles_tensor::sub_scal_prod(const Giles_tensor &rhs, double sca) {
  return *this;
}
Giles_tensor &Giles_tensor::add_scal_div(const Giles_tensor &rhs, double sca) {
  return *this;
}
Giles_tensor &Giles_tensor::sub_scal_div(const Giles_tensor &rhs, double sca) {
  return *this;
}
Giles_tensor &Giles_tensor::operator*=(const Giles_tensor &rhs) {
  return *this;
}
Giles_tensor Giles_tensor::operator*(const Giles_tensor &rhs) const {
  Giles_tensor result(m_width, std::max(m_degree, rhs.m_degree));
  pass_through op;
  inplace_multiplication_impl(result, result.transposed_data,
                              m_width, m_degree, rhs.m_degree,
                              *this, rhs, size_array, op);
  return result;
}
Giles_tensor &Giles_tensor::mul_scal_mul(const Giles_tensor &rhs, double sca) {
  return *this;
}
Giles_tensor &Giles_tensor::mul_scal_div(const Giles_tensor &rhs, double sca) {
  return *this;
}
bool Giles_tensor::operator==(const Giles_tensor &other) const {
  return false;
}
bool Giles_tensor::operator!=(const Giles_tensor &other) const {
  return false;
}

std::ostream &operator<<(std::ostream &os, const Giles_tensor &arg) {
  return os;
}
