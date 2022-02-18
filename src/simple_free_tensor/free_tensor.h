//
// Created by sam on 06/11/2021.
//

#ifndef TENSOR_PLAYGROUND_FREE_TENSOR_H
#define TENSOR_PLAYGROUND_FREE_TENSOR_H

#include <iosfwd>
#include <vector>

class free_tensor : std::vector<double> {
    using base_type = std::vector<double>;

   public:
    using typename base_type::size_type;
    using base_type::operator[];
    using base_type::assign;
    using base_type::begin;
    using base_type::data;
    using base_type::emplace_back;
    using base_type::end;
    using base_type::size;
    using base_type::operator=;

    using degree_type = unsigned;

   private:
    degree_type m_degree;
    degree_type m_width;
    std::vector<size_type> size_array;

    free_tensor(degree_type width);

    void populate_size_array(degree_type maxd);

   public:
    free_tensor() = delete;

    free_tensor(degree_type width, degree_type depth);

    free_tensor(degree_type width, degree_type depth, std::initializer_list<double> args);

    double *range_begin();

    double *range_end();

    const double *range_begin() const;

    const double *range_end() const;

    degree_type degree() const;
    degree_type width() const;

    free_tensor &operator+=(const free_tensor &rhs);
    free_tensor &operator-=(const free_tensor &rhs);
    free_tensor &operator*=(double rhs);
    free_tensor &operator/=(double rhs);

    free_tensor operator+(const free_tensor &rhs) const;
    free_tensor operator-(const free_tensor &rhs) const;
    free_tensor operator*(double rhs) const;
    free_tensor operator/(double rhs) const;

    free_tensor &add_scal_prod(const free_tensor &rhs, double sca);
    free_tensor &sub_scal_prod(const free_tensor &rhs, double sca);
    free_tensor &add_scal_div(const free_tensor &rhs, double sca);
    free_tensor &sub_scal_div(const free_tensor &rhs, double sca);

    free_tensor &operator*=(const free_tensor &rhs);
    free_tensor operator*(const free_tensor &rhs) const;

    free_tensor &mul_scal_mul(const free_tensor &rhs, double sca);
    free_tensor &mul_scal_div(const free_tensor &rhs, double sca);

    bool operator==(const free_tensor &other) const;
    bool operator!=(const free_tensor &other) const;

    friend std::ostream &operator<<(std::ostream &os, const free_tensor &arg);

    friend free_tensor exp(const free_tensor &arg);

    friend double L1Norm(const free_tensor &arg) noexcept;
    friend double LInfNorm(const free_tensor &arg) noexcept;
};

std::ostream &operator<<(std::ostream &os, const free_tensor &arg);
double L1Norm(const free_tensor &arg) noexcept;
double LInfNorm(const free_tensor &arg) noexcept;

#endif//TENSOR_PLAYGROUND_FREE_TENSOR_H
