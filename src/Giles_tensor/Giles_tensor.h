//
// Created by sam on 07/11/2021.
//

#ifndef TENSOR_PLAYGROUND_GILES_TENSOR_H
#define TENSOR_PLAYGROUND_GILES_TENSOR_H

#include <vector>
#include <iosfwd>

class Giles_tensor : std::vector<double>
{
    using base_type = std::vector<double>;

public:
    using typename base_type::size_type;
    using base_type::operator[];
    using base_type::begin;
    using base_type::end;
    using base_type::emplace_back;
    using base_type::assign;
    using base_type::data;
    using base_type::size;
    using base_type::operator=;

    using degree_type = unsigned;

private:

    degree_type m_degree;
    degree_type m_width;
    std::vector<size_type> size_array;
    std::vector<double> transposed_data;

    Giles_tensor(degree_type width);

    void populate_size_array(degree_type maxd);

public:
    Giles_tensor() = delete;

    Giles_tensor(degree_type
    width,
    degree_type depth
    );

    Giles_tensor(degree_type
    width,
    degree_type depth, std::initializer_list<double>
    args);

    double* range_begin();

    double* range_end();

    const double* range_begin() const;

    const double* range_end() const;

    const double* trange_begin() const;
    const double* trange_end() const;

    degree_type degree() const;
    degree_type width() const;

    Giles_tensor& operator+=(const Giles_tensor& rhs);
    Giles_tensor& operator-=(const Giles_tensor& rhs);
    Giles_tensor& operator*=(double rhs);
    Giles_tensor& operator/=(double rhs);

    Giles_tensor operator+(const Giles_tensor& rhs) const;
    Giles_tensor operator-(const Giles_tensor& rhs) const;
    Giles_tensor operator*(double rhs) const;
    Giles_tensor operator/(double rhs) const;

    Giles_tensor& add_scal_prod(const Giles_tensor& rhs, double sca);
    Giles_tensor& sub_scal_prod(const Giles_tensor& rhs, double sca);
    Giles_tensor& add_scal_div(const Giles_tensor& rhs, double sca);
    Giles_tensor& sub_scal_div(const Giles_tensor& rhs, double sca);

    Giles_tensor& operator*=(const Giles_tensor& rhs);
    Giles_tensor operator*(const Giles_tensor& rhs) const;

    Giles_tensor& mul_scal_mul(const Giles_tensor& rhs, double sca);
    Giles_tensor& mul_scal_div(const Giles_tensor& rhs, double sca);

    bool operator==(const Giles_tensor& other) const;
    bool operator!=(const Giles_tensor& other) const;

    friend std::ostream& operator<<(std::ostream& os, const Giles_tensor& arg);

    friend Giles_tensor exp(const Giles_tensor& arg);

    friend double L1Norm(const Giles_tensor& arg) noexcept;
    friend double LInfNorm(const Giles_tensor& arg) noexcept;

};


std::ostream& operator<<(std::ostream& os, const Giles_tensor& arg);
double L1Norm(const Giles_tensor& arg) noexcept;
double LInfNorm(const Giles_tensor& arg) noexcept;


#endif //TENSOR_PLAYGROUND_GILES_TENSOR_H
