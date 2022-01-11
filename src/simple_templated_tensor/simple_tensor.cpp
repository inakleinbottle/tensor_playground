//
// Created by sam on 29/09/2021.
//

#include "simple_template_tensor.h"

#include <cmath>

simple_tensor& simple_tensor::operator+=(const simple_tensor& rhs)
{
    double* p = range_begin();
    const double* rp = rhs.range_begin();
    for (size_t i = 0; i < tensor_alg_size(DEPTH); ++i) {
        *(p++) += *(rp++);
    }
    return *this;
}
simple_tensor& simple_tensor::operator-=(const simple_tensor& rhs)
{
    double* p = range_begin();
    const double* rp = rhs.range_begin();
    for (size_t i = 0; i<tensor_alg_size(DEPTH); ++i) {
        *(p++) -= *(rp++);
    }
    return *this;
}
simple_tensor& simple_tensor::operator*=(double rhs)
{
    double* p = range_begin();;
    for (size_t i = 0; i<tensor_alg_size(DEPTH); ++i) {
        *(p++) *= rhs;
    }
    return *this;
}
simple_tensor& simple_tensor::operator/=(double rhs)
{
    double* p = range_begin();;
    for (size_t i = 0; i<tensor_alg_size(DEPTH); ++i) {
        *(p++) /= rhs;
    }
    return *this;
}
simple_tensor simple_tensor::operator+(const simple_tensor& rhs) const
{
    simple_tensor result(*this);
    result += rhs;
    return result;
}
simple_tensor simple_tensor::operator-(const simple_tensor& rhs) const
{
    simple_tensor result(*this);
    result -= rhs;
    return result;
}
simple_tensor simple_tensor::operator*(double rhs) const
{
    simple_tensor result(*this);
    result *= rhs;
    return result;
}
simple_tensor simple_tensor::operator/(double rhs) const
{
    simple_tensor result(*this);
    result /= rhs;
    return result;
}
simple_tensor& simple_tensor::add_scal_prod(const simple_tensor& rhs, double sca)
{
    double* p = range_begin();
    const double* rp = rhs.range_begin();
    for (size_t i = 0; i<tensor_alg_size(DEPTH); ++i) {
        *(p++) += *(rp++) * sca;
    }
    return *this;
}
simple_tensor& simple_tensor::sub_scal_prod(const simple_tensor& rhs, double sca)
{
    double* p = range_begin();
    const double* rp = rhs.range_begin();
    for (size_t i = 0; i<tensor_alg_size(DEPTH); ++i) {
        *(p++) -= *(rp++) * sca;
    }
    return *this;
}
simple_tensor& simple_tensor::add_scal_div(const simple_tensor& rhs, double sca)
{
    double* p = range_begin();
    const double* rp = rhs.range_begin();
    for (size_t i = 0; i<tensor_alg_size(DEPTH); ++i) {
        *(p++) += *(rp++) / sca;
    }
    return *this;
}
simple_tensor& simple_tensor::sub_scal_div(const simple_tensor& rhs, double sca)
{
    double* p = range_begin();
    const double* rp = rhs.range_begin();
    for (size_t i = 0; i<tensor_alg_size(DEPTH); ++i) {
        *(p++) -= *(rp++) / sca;
    }
    return *this;
}

template<typename Op>
inline void
inplace_multiplication_impl(
        simple_tensor& out,
        const simple_tensor& lhs,
        const simple_tensor& rhs,
        const std::vector<size_type>& start_of_degree_array,
        Op op
        ) noexcept
{
    int rhs_deg;

    double* out_ptr;
    const double* lhs_p, *rhs_p;
    for (int out_deg = DEPTH; out_deg >= 0; --out_deg) {
        for (int lhs_deg = out_deg; lhs_deg >= 0; --lhs_deg) {
            rhs_deg = out_deg - lhs_deg;

            out_ptr = out.data()+start_of_degree_array[out_deg];
            lhs_p = lhs.range_begin()+start_of_degree_array[lhs_deg];
            rhs_p = rhs.range_begin()+start_of_degree_array[rhs_deg];

            for (int i = 0; i<start_of_degree_array[lhs_deg+1]-start_of_degree_array[lhs_deg]; ++i) {
                for (int j = 0; j<start_of_degree_array[rhs_deg+1]-start_of_degree_array[rhs_deg]; ++j) {
                    *(out_ptr++) += op(lhs_p[i]*rhs_p[j]);
                }
            }

        }
    }

}




struct pass_through
{
    constexpr pass_through() = default;
    constexpr double operator()(double arg) const { return arg; }
};

struct scalar_minus
{
    constexpr scalar_minus() = default;
    constexpr double operator()(double arg) const { return -arg; }
};

struct post_mul
{
    explicit constexpr post_mul(double s)
            :sca(s) { }

    constexpr double operator()(double arg) const { return arg * sca; }
private:
    double sca;
};

struct post_div
{
    explicit constexpr post_div(double s) : sca(s) {}

    constexpr double operator()(double arg) const { return arg / sca; }
private:
    double sca;
};





simple_tensor& simple_tensor::operator*=(const simple_tensor& rhs)
{
    simple_tensor tmp;
    pass_through op;
    inplace_multiplication_impl(tmp, *this, rhs, op);
    base_type::swap(tmp);
    return *this;
}
simple_tensor simple_tensor::operator*(const simple_tensor& rhs) const
{
    simple_tensor result;
    pass_through op;
    inplace_multiplication_impl(result, *this, rhs, op);
    return result;
}
simple_tensor& simple_tensor::mul_scal_mul(const simple_tensor& rhs, double sca)
{
    post_mul op(sca);
    simple_tensor tmp;
    inplace_multiplication_impl(tmp, *this, rhs, op);
    base_type::swap(tmp);
    return *this;
}
simple_tensor& simple_tensor::mul_scal_div(const simple_tensor& rhs, double sca)
{
    post_div op(sca);
    simple_tensor tmp;
    inplace_multiplication_impl(tmp, *this, rhs, op);
    base_type::swap(tmp);
    return *this;
}

bool simple_tensor::operator==(const simple_tensor& other) const
{
    for (size_t i=0; i <tensor_alg_size(DEPTH); ++i) {
        if ((*this)[i] != other[i])
            return false;
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, const simple_tensor& arg)
{
    os << "{ ";
    for (auto& v : arg) {
        os << v << ' ';
    }
    return os << '}';
}


simple_tensor exp(const simple_tensor& arg)
{
    simple_tensor result {1.0};
    for (deg_t d = DEPTH; d >= 1; --d) {
        result.mul_scal_div(arg, static_cast<double>(d));
        result[0] += 1.0;
    }
    return result;
}

double L1Norm(const simple_tensor& arg) noexcept
{
    double result = 0.0;
    for (const auto& v : arg) {
        result = std::abs(v);
    }
    return result;
}

double LInfNorm(const simple_tensor& arg) noexcept
{
    double result = 0.0;
    for (const auto& v: arg) {
        auto av = std::abs(v);
        result = (av > result) ? av : result;
    }
    return result;

}
