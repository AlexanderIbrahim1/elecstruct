#pragma once

#include <cstdint>
#include <cmath>
#include <tuple>

#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/geometry.hpp"
#include "elecstruct/orbitals.hpp"
#include "elecstruct/basis/gaussian_info.hpp"

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
    const auto expon_scaling = - (info0.exponent * info1.exponent) / (info0.exponent + info1.exponent);
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

}

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

    return {new_centre, {new_coefficient, new_exponent}};
}

// inline auto gaussian_norm(const elec::AngularMomentumNumbers& angular_momenta, const elec::GaussianInfo& info) -> double
// {
// 
// }

}  // namespace elec
