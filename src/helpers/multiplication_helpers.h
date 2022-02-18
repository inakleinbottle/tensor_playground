//
// Created by sam on 06/11/2021.
//

#ifndef TENSOR_PLAYGROUND_MULTIPLICATION_HELPERS_H
#define TENSOR_PLAYGROUND_MULTIPLICATION_HELPERS_H

struct pass_through {
  constexpr pass_through() = default;
  constexpr double operator()(double arg) const { return arg; }
};

struct scalar_minus {
  constexpr scalar_minus() = default;
  constexpr double operator()(double arg) const { return -arg; }
};

struct post_mul {
  explicit constexpr post_mul(double s)
      : sca(s) {}

  constexpr double operator()(double arg) const { return arg * sca; }

 private:
  double sca;
};

struct post_div {
  explicit constexpr post_div(double s)
      : sca(s) {}

  constexpr double operator()(double arg) const { return arg / sca; }

 private:
  double sca;
};

struct bin_add_inplace {
  inline void operator()(double &lhs, const double &rhs) noexcept {
    lhs += rhs;
  }
};

struct bin_sub_inplace {
  inline void operator()(double &lhs, const double &rhs) noexcept {
    lhs -= rhs;
  }
};

struct bin_fmadd_inplace {
  double factor;

  explicit constexpr bin_fmadd_inplace(double f)
      : factor(f) {}

  inline void operator()(double &lhs, const double &rhs) noexcept {
    lhs += rhs * factor;
  }
};

struct bin_fmsub_inplace {
  double factor;

  explicit constexpr bin_fmsub_inplace(double f)
      : factor(f) {}

  inline void operator()(double &lhs, const double &rhs) const noexcept {
    lhs -= rhs * factor;
  }
};

struct bin_fdadd_inplace {
  double factor;

  explicit constexpr bin_fdadd_inplace(double f)
      : factor(f) {}

  inline void operator()(double &lhs, const double &rhs) const noexcept {
    lhs += rhs / factor;
  }
};

struct bin_fdsub_inplace {
  double factor;

  explicit constexpr bin_fdsub_inplace(double f)
      : factor(f) {}

  inline void operator()(double &lhs, const double &rhs) const noexcept {
    lhs -= rhs / factor;
  }
};

struct u_add {
  constexpr double operator()(const double &arg) const noexcept {
    return arg;
  }
};

struct u_sub {
  constexpr double operator()(const double &arg) const noexcept {
    return -arg;
  }
};

struct u_fmadd {
  double factor;

  explicit constexpr u_fmadd(double f)
      : factor(f) {}

  constexpr double operator()(const double &arg) const noexcept {
    return arg * factor;
  }
};

struct u_fmsub {
  double factor;

  explicit constexpr u_fmsub(double f)
      : factor(f) {}

  constexpr double operator()(const double &arg) const noexcept {
    return -arg * factor;
  }
};

struct u_fdadd {
  double factor;

  explicit constexpr u_fdadd(double f)
      : factor(f) {}

  constexpr double operator()(const double &arg) const noexcept {
    return arg / factor;
  }
};

struct u_fdsub {
  double factor;

  explicit constexpr u_fdsub(double f)
      : factor(f) {}

  constexpr double operator()(const double &arg) const noexcept {
    return -arg / factor;
  }
};

#endif//TENSOR_PLAYGROUND_MULTIPLICATION_HELPERS_H
