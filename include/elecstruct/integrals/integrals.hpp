#pragma once

#include <cmath>
#include <cstdint>
#include <tuple>

#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/geometry.hpp"
#include "elecstruct/orbitals.hpp"

namespace impl_elec
{

inline auto gaussian_product_coefficient(
    const coord::Cartesian3D& centre0,
    const coord::Cartesian3D& centre1,
    const elec::GaussianInfo& info0,
    const elec::GaussianInfo& info1
) -> double
{
    const auto diff = centre1 - centre0;
    const auto norm_sq = coord::dot_product(diff, diff);
    const auto expon_scaling = -(info0.exponent * info1.exponent) / (info0.exponent + info1.exponent);
    const auto coefficients = info0.coefficient * info1.coefficient;

    return coefficients * std::exp(expon_scaling * norm_sq);
}

inline auto double_factorial(std::uint64_t n) -> std::uint64_t
{
    auto result = std::uint64_t {1};

    for (std::uint64_t value {n}; value >= 2; value -= 2) {
        result *= value;
    }

    return result;
}

inline auto angular_momentum_double_factorial_term(std::uint64_t angular_momentum) -> double
{
    if (angular_momentum == 0) {
        return 1.0;
    }
    else {
        return static_cast<double>(double_factorial(2 * angular_momentum - 1));
    }
}

}  // namespace impl_elec

namespace elec
{

inline auto gaussian_product(
    const coord::Cartesian3D& centre0,
    const coord::Cartesian3D& centre1,
    const GaussianInfo& info0,
    const GaussianInfo& info1
) -> std::tuple<coord::Cartesian3D, GaussianInfo>
{
    const auto new_centre = (centre0 * info0.exponent + centre1 * info1.exponent) / (info0.exponent + info1.exponent);
    const auto new_coefficient = impl_elec::gaussian_product_coefficient(centre0, centre1, info0, info1);
    const auto new_exponent = info0.exponent + info1.exponent;

    return {
        new_centre, {new_coefficient, new_exponent}
    };
}

inline auto gaussian_norm(const elec::AngularMomentumNumbers& angular_momenta, double gauss_exponent) -> double
{
    const auto ang_mom_sum = static_cast<double>(total_angular_momentum(angular_momenta));

    const auto gauss1d_component = std::pow(2.0 * gauss_exponent / M_PI, 3.0 / 4.0);
    const auto ang_mom_numerator = std::pow(4.0 * gauss_exponent, ang_mom_sum / 2.0);

    // NOTE: maybe the behaviour of scipy changed since the reference was created?
    // - the author used `scipy.misc.factorial2()`, which returns 0 for negative numbers (according to docs I've seen)
    // - if the angular momentum component is 0, then the argument `2l - 1` to the double factorial is negative
    //   - this would make the denominator 0
    // - this probably isn't what is supposed to happen
    const auto ang_mom_denom_x = impl_elec::angular_momentum_double_factorial_term(angular_momenta.x);
    const auto ang_mom_denom_y = impl_elec::angular_momentum_double_factorial_term(angular_momenta.y);
    const auto ang_mom_denom_z = impl_elec::angular_momentum_double_factorial_term(angular_momenta.z);
    const auto ang_mom_denominator = std::sqrt(ang_mom_denom_x * ang_mom_denom_y * ang_mom_denom_z);

    return gauss1d_component * ang_mom_numerator / ang_mom_denominator;
}

}  // namespace elec
