//
// Created by sam on 29/09/2021.
//

#ifndef TENSOR_INVERSE_SIMPLE_TENSOR_H
#define TENSOR_INVERSE_SIMPLE_TENSOR_H

#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
#include <multiplication_helpers.h>
#include <tensor_size_info.h>


template <unsigned Width, unsigned Depth, typename Coeffs=double>
class simple_template_tensor : std::vector<Coeffs>
{
protected:
    using base_type = std::vector<Coeffs>;
public:

    using base_type::operator[];
    using base_type::begin;
    using base_type::end;
    using base_type::emplace_back;
    using base_type::assign;
    using base_type::data;
    using base_type::size;
    using typename base_type::size_type;

    using base_type::operator=;

    static const std::vector<size_type> size_array;

    simple_template_tensor() : base_type(tensor_alg_size(Width, Depth))
    {}

    simple_template_tensor(std::initializer_list<Coeffs> args) : base_type(args)
    {
        base_type::resize(tensor_alg_size(Width, Depth));
    }

    explicit simple_template_tensor(size_t index, Coeffs coeff=1.0) : base_type(tensor_alg_size(Width, Depth))
    {
        base_type::operator[](index) = coeff;
    }

    template <typename... Args>
    explicit simple_template_tensor(Args... args) : base_type(std::forward<Args>(args)...)
    {
        base_type::resize(tensor_alg_size(Width, Depth));
    }

    double* range_begin()
    {
        return &*base_type::begin();
    }

    double* range_end()
    {
        return &*base_type::end();
    }

    const double* range_begin() const
    {
        return &*base_type::begin();
    }

    const double* range_end() const
    {
        return &*base_type::end();
    }

private:

    template<typename Op>
    inline void
    inplace_multiplication_impl(
            simple_template_tensor& out,
            const simple_template_tensor& lhs,
            const simple_template_tensor& rhs,
            const std::vector<size_type>& start_of_degree_array,
            Op op) const noexcept
    {
        int rhs_deg;

        double* out_ptr;
        const double* lhs_p, * rhs_p;
        for (int out_deg = Depth; out_deg>=0; --out_deg) {


            for (int lhs_deg = out_deg; lhs_deg>=0; --lhs_deg) {
                rhs_deg = out_deg-lhs_deg;

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

public:

    simple_template_tensor& operator+=(const simple_template_tensor& rhs)
    {
        double* p = range_begin();
        const double* rp = rhs.range_begin();
        for (size_t i = 0; i<tensor_alg_size(Width, Depth); ++i) {
            *(p++) += *(rp++);
        }
        return *this;
    }

    simple_template_tensor& operator-=(const simple_template_tensor& rhs)
    {
        double* p = range_begin();
        const double* rp = rhs.range_begin();
        for (size_t i = 0; i<tensor_alg_size(Width, Depth); ++i) {
            *(p++) -= *(rp++);
        }
        return *this;
    }

    simple_template_tensor& operator*=(double rhs)
    {
        double* p = range_begin();;
        for (size_t i = 0; i<tensor_alg_size(Width, Depth); ++i) {
            *(p++) *= rhs;
        }
        return *this;
    }

    simple_template_tensor& operator/=(double rhs)
    {
        double* p = range_begin();;
        for (size_t i = 0; i<tensor_alg_size(Width, Depth); ++i) {
            *(p++) /= rhs;
        }
        return *this;
    }

    simple_template_tensor operator+(const simple_template_tensor& rhs) const
    {
        simple_template_tensor result(*this);
        result += rhs;
        return result;
    }

    simple_template_tensor operator-(const simple_template_tensor& rhs) const
    {
        simple_template_tensor result(*this);
        result -= rhs;
        return result;
    }

    simple_template_tensor operator*(double rhs) const
    {
        simple_template_tensor result(*this);
        result *= rhs;
        return result;
    }

    simple_template_tensor operator/(double rhs) const
    {
        simple_template_tensor result(*this);
        result /= rhs;
        return result;
    }

    simple_template_tensor& add_scal_prod(const simple_template_tensor& rhs, double sca)
    {
        double* p = range_begin();
        const double* rp = rhs.range_begin();
        for (size_t i = 0; i<tensor_alg_size(Width, Depth); ++i) {
            *(p++) += *(rp++)*sca;
        }
        return *this;
    }

    simple_template_tensor& sub_scal_prod(const simple_template_tensor& rhs, double sca)
    {
        double* p = range_begin();
        const double* rp = rhs.range_begin();
        for (size_t i = 0; i<tensor_alg_size(Width, Depth); ++i) {
            *(p++) -= *(rp++)*sca;
        }
        return *this;
    }

    simple_template_tensor& add_scal_div(const simple_template_tensor& rhs, double sca)
    {
        double* p = range_begin();
        const double* rp = rhs.range_begin();
        for (size_t i = 0; i<tensor_alg_size(Width, Depth); ++i) {
            *(p++) += *(rp++)/sca;
        }
        return *this;
    }

    simple_template_tensor& sub_scal_div(const simple_template_tensor& rhs, double sca)
    {
        double* p = range_begin();
        const double* rp = rhs.range_begin();
        for (size_t i = 0; i<tensor_alg_size(Width, Depth); ++i) {
            *(p++) -= *(rp++)/sca;
        }
        return *this;
    }


    simple_template_tensor& operator*=(const simple_template_tensor& rhs)
    {
        simple_template_tensor tmp;
        pass_through op;
        inplace_multiplication_impl(tmp, *this, rhs, size_array, op);
        base_type::swap(tmp);
        return *this;
    }

    simple_template_tensor operator*(const simple_template_tensor& rhs) const
    {
        simple_template_tensor result;
        pass_through op;
        inplace_multiplication_impl(result, *this, rhs, size_array, op);
        return result;
    }

    simple_template_tensor& mul_scal_mul(const simple_template_tensor& rhs, double sca)
    {
        post_mul op(sca);
        simple_template_tensor tmp;
        inplace_multiplication_impl(tmp, *this, rhs, size_array, op);
        base_type::swap(tmp);
        return *this;
    }

    simple_template_tensor& mul_scal_div(const simple_template_tensor& rhs, double sca)
    {
        post_div op(sca);
        simple_template_tensor tmp;
        inplace_multiplication_impl(tmp, *this, rhs, size_array, op);
        base_type::swap(tmp);
        return *this;
    }

    bool operator==(const simple_template_tensor& other) const
    {
        for (size_t i = 0; i<tensor_alg_size(Width, Depth); ++i) {
            if ((*this)[i]!=other[i])
                return false;
        }
        return true;
    }

    bool operator!=(const simple_template_tensor& other) const
    {
        return !(*this == other);
    }

    friend std::ostream& operator<<(std::ostream& os, const simple_template_tensor& arg)
    {
        os << "{ ";
        for (auto& v: arg) {
            os << v << ' ';
        }
        return os << '}';
    }

    friend simple_template_tensor exp(const simple_template_tensor& arg)
    {
        simple_template_tensor result{1.0};
        for (unsigned d = Depth; d>=1; --d) {
            result.mul_scal_div(arg, static_cast<double>(d));
            result[0] += 1.0;
        }
        return result;
    }

    friend double L1Norm(const simple_template_tensor& arg) noexcept
    {
        double result = 0.0;
        for (const auto& v: arg) {
            result = std::abs(v);
        }
        return result;
    }

    friend double LInfNorm(const simple_template_tensor& arg) noexcept
    {
        double result = 0.0;
        for (const auto& v: arg) {
            auto av = std::abs(v);
            result = (av>result) ? av : result;
        }
        return result;
    }

};

template <unsigned Width, unsigned Depth, typename Coeffs>
const std::vector<typename simple_template_tensor<Width, Depth, Coeffs>::size_type>
simple_template_tensor<Width, Depth, Coeffs>::size_array = []() {
    std::vector<typename simple_template_tensor<Width, Depth, Coeffs>::size_type> result{0};
    result.reserve(Depth+2);
    for (int i=0; i <= Depth+1; ++i)
        result.emplace_back(result.back() + power(Width, i));
    return result;
}();



#endif //TENSOR_INVERSE_SIMPLE_TENSOR_H
