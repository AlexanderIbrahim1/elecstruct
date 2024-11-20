#pragma once

#include <cmath>
#include <cstdint>
#include <tuple>

#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/geometry.hpp"
#include "elecstruct/mathtools/factorial.hpp"
#include "elecstruct/orbitals.hpp"

namespace elec::math
{

namespace impl_elec
{

/*
    NOTE: this function calculates product of two gaussians with coefficients of 1 each
      - the actual coefficients of the two input gaussians are handled by the functions that
        calculate the matrices
*/
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

    return std::exp(expon_scaling * norm_sq);
}

}  // namespace elec::impl_elec


/*
    NOTE: this function calculates product of two gaussians with coefficients of 1 each
      - the actual coefficients of the two input gaussians are handled by the functions that
        calculate the matrices
*/
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
        new_centre, GaussianInfo {new_coefficient, new_exponent}
    };
}

inline auto gaussian_norm(const elec::AngularMomentumNumbers& angular_momenta, double gauss_exponent) -> double
{
    const auto angmom_sum = static_cast<double>(total_angular_momentum(angular_momenta));

    const auto gauss1d_component = std::pow(2.0 * gauss_exponent / M_PI, 3.0 / 4.0);
    const auto angmom_numerator = std::pow(4.0 * gauss_exponent, angmom_sum / 2.0);

    // NOTE: maybe the behaviour of scipy changed since the reference was created?
    // - the author used `scipy.misc.factorial2()`, which returns 0 for negative numbers (according to docs I've seen)
    // - if the angular momentum component is 0, then the argument `2l - 1` to the double factorial is negative
    //   - this would make the denominator 0
    // - this probably isn't what is supposed to happen
    const auto angmom_denom_x = elec::math::double_factorial(2 * angular_momenta.x - 1);
    const auto angmom_denom_y = elec::math::double_factorial(2 * angular_momenta.y - 1);
    const auto angmom_denom_z = elec::math::double_factorial(2 * angular_momenta.z - 1);
    const auto angmom_denominator = std::sqrt(angmom_denom_x * angmom_denom_y * angmom_denom_z);

    return gauss1d_component * angmom_numerator / angmom_denominator;
}

}  // namespace elec::math
