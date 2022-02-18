//
// Created by sam on 06/11/2021.
//

#include "free_tensor.h"
#include <cassert>
#include <functional>
#include <iostream>
#include <multiplication_helpers.h>
#include <tensor_size_info.h>
#include <utility>

using namespace playground;
namespace {

using size_type = typename free_tensor::size_type;
using degree_type = typename free_tensor::degree_type;

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
    unsigned width,
    unsigned lhs_depth,
    unsigned rhs_depth,
    const free_tensor &lhs,
    const free_tensor &rhs,
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

}// namespace

free_tensor::free_tensor(free_tensor::degree_type width)
    : m_width(width), m_degree(0), base_type(), size_array() {
}

free_tensor::free_tensor(free_tensor::degree_type width, free_tensor::degree_type depth)
    : m_width(width), m_degree(depth), base_type(tensor_alg_size(width, depth)), size_array() {
    populate_size_array(depth + 2);
}

free_tensor::free_tensor(free_tensor::degree_type width, free_tensor::degree_type depth,
                         std::initializer_list<double> args)
    : m_width(width), m_degree(depth), base_type(args), size_array() {
    resize(tensor_alg_size(width, depth));
    populate_size_array(depth + 2);
}

double *free_tensor::range_begin() {
    return data();
}
double *free_tensor::range_end() {
    return data() + size();
}
const double *free_tensor::range_begin() const {
    return data();
}
const double *free_tensor::range_end() const {
    return data() + size();
}
free_tensor &free_tensor::operator+=(const free_tensor &rhs) {
    assert(m_width == rhs.m_width);
    bin_add_inplace bo;
    u_add uo;
    do_inplace_binary_operation(*this, rhs.size(), rhs.data(), bo, uo);
    return *this;
}
free_tensor &free_tensor::operator-=(const free_tensor &rhs) {
    assert(m_width == rhs.m_width);
    bin_sub_inplace bo;
    u_sub uo;
    do_inplace_binary_operation(*this, rhs.size(), rhs.data(), bo, uo);
    return *this;
}
free_tensor &free_tensor::operator*=(double rhs) {
    for (auto &v : *this) {
        v *= rhs;
    }

    return *this;
}
free_tensor &free_tensor::operator/=(double rhs) {
    for (auto &v : *this) {
        v /= rhs;
    }

    return *this;
}

free_tensor free_tensor::operator+(const free_tensor &rhs) const {
    free_tensor result(m_width);
    assert(m_width == rhs.m_width);
    auto rhs_size = rhs.size();
    auto lhs_size = size();

    auto max_dim = std::max(lhs_size, rhs_size);
    result.reserve(max_dim);

    do_binary_operation(
        result,
        lhs_size,
        rhs_size,
        data(),
        rhs.data(),
        std::plus<>(),
        u_add(),
        u_add());
    return result;
}
free_tensor free_tensor::operator-(const free_tensor &rhs) const {
    free_tensor result(m_width);
    assert(m_width == rhs.m_width);
    auto rhs_size = rhs.size();
    auto lhs_size = size();
    auto max_dim = std::max(lhs_size, rhs_size);

    result.reserve(max_dim);

    do_binary_operation(
        result,
        lhs_size,
        rhs_size,
        data(),
        rhs.data(),
        std::minus<>(),
        u_add(),
        u_sub());

    return result;
}

free_tensor free_tensor::operator*(double rhs) const {
    free_tensor result(m_width);
    result.reserve(size());
    for (const auto &v : *this) {
        result.emplace_back(v * rhs);
    }
    return result;
}

free_tensor free_tensor::operator/(double rhs) const {
    free_tensor result(m_width);
    result.reserve(size());
    for (const auto &v : *this) {
        result.emplace_back(v / rhs);
    }
    return result;
}
free_tensor &free_tensor::add_scal_prod(const free_tensor &rhs, double sca) {
    assert(m_width == rhs.m_width);
    bin_fmadd_inplace bo(sca);
    u_fmadd uo(sca);
    do_inplace_binary_operation(*this, rhs.size(), rhs.data(), bo, uo);
    return *this;
}
free_tensor &free_tensor::sub_scal_prod(const free_tensor &rhs, double sca) {
    assert(m_width == rhs.m_width);
    bin_fmsub_inplace bo(sca);
    u_fmsub uo(sca);
    do_inplace_binary_operation(*this, rhs.size(), rhs.data(), bo, uo);
    return *this;
}
free_tensor &free_tensor::add_scal_div(const free_tensor &rhs, double sca) {
    assert(m_width == rhs.m_width);
    bin_fmadd_inplace bo(1.0 / sca);
    u_fmadd uo(1.0 / sca);
    do_inplace_binary_operation(*this, rhs.size(), rhs.data(), bo, uo);
    return *this;
}
free_tensor &free_tensor::sub_scal_div(const free_tensor &rhs, double sca) {
    assert(m_width == rhs.m_width);
    bin_fmsub_inplace bo(1.0 / sca);
    u_fmsub uo(1.0 / sca);
    do_inplace_binary_operation(*this, rhs.size(), rhs.data(), bo, uo);
    return *this;
}
free_tensor &free_tensor::operator*=(const free_tensor &rhs) {
    free_tensor tmp(m_width, std::max(m_degree, rhs.m_degree));
    pass_through op;
    inplace_multiplication_impl(tmp,
                                m_width, m_degree, rhs.m_degree,
                                *this, rhs, size_array, op);
    base_type::swap(tmp);
    return *this;
}

free_tensor free_tensor::operator*(const free_tensor &rhs) const {
    free_tensor result(m_width, std::max(m_degree, rhs.m_degree));
    pass_through op;
    inplace_multiplication_impl(result,
                                m_width, m_degree, rhs.m_degree,
                                *this, rhs, size_array, op);
    return result;
}
free_tensor &free_tensor::mul_scal_mul(const free_tensor &rhs, double sca) {
    post_mul op(sca);
    free_tensor tmp(m_width, std::max(m_degree, rhs.m_degree));
    inplace_multiplication_impl(tmp,
                                m_width, m_degree, rhs.m_degree,
                                *this, rhs, size_array, op);
    base_type::swap(tmp);
    return *this;
}
free_tensor &free_tensor::mul_scal_div(const free_tensor &rhs, double sca) {
    post_div op(sca);
    free_tensor tmp(m_width, std::max(m_degree, rhs.m_degree));
    inplace_multiplication_impl(tmp,
                                m_width, m_degree, rhs.m_degree,
                                *this, rhs, size_array, op);
    base_type::swap(tmp);
    return *this;
}

bool free_tensor::operator==(const free_tensor &other) const {
    if (m_width != other.m_width) {
        return false;
    }

    auto min_size = std::min(size(), other.size());

    for (size_type i = 0; i < min_size; ++i) {
        if (operator[](i) != other[i]) {
            return false;
        }
    }

    for (size_type i = min_size; i < size(); ++i) {
        if (operator[](i) != 0.0) {
            return false;
        }
    }

    for (size_type i = min_size; i < other.size(); ++i) {
        if (other[i] != 0.0) {
            return false;
        }
    }

    return true;
}

bool free_tensor::operator!=(const free_tensor &other) const {
    return !(*this == other);
}

free_tensor::degree_type free_tensor::degree() const {
    return m_degree;
}

free_tensor::degree_type free_tensor::width() const {
    return m_width;
}

void free_tensor::populate_size_array(degree_type maxd) {
    size_array.reserve(maxd);
    if (size_array.empty())
        size_array.push_back(0);

    auto sz = size_array.size();

    for (auto i = sz - 1; i < maxd; ++i) {
        size_array.emplace_back(size_array[i] + power(m_width, i));
    }
}

std::ostream &operator<<(std::ostream &os, const free_tensor &arg) {
    os << "{ ";
    for (auto &v : arg) {
        os << v << ' ';
    }
    return os << '}';
}

double L1Norm(const free_tensor &arg) noexcept {
    double result = 0.0;
    for (const auto &v : arg) {
        result = std::abs(v);
    }
    return result;
}

double LInfNorm(const free_tensor &arg) noexcept {
    double result = 0.0;
    for (const auto &v : arg) {
        auto av = std::abs(v);
        result = (av > result) ? av : result;
    }
    return result;
}
