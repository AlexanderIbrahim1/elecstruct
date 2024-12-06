#include <cmath>
#include <tuple>

#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/geometry.hpp"
#include "elecstruct/mathtools/misc.hpp"
#include "elecstruct/orbitals.hpp"

#include "elecstruct/mathtools/gaussian.hpp"

namespace
{

/*
    NOTE: this function calculates product of two gaussians with coefficients of 1 each
      - the actual coefficients of the two input gaussians are handled by the functions that
        calculate the matrices
*/
auto gaussian_product_coefficient(
    const coord::Cartesian3D& centre0,
    const coord::Cartesian3D& centre1,
    double exponent0,
    double exponent1
) -> double
{
    const auto diff = centre1 - centre0;
    const auto norm_sq = coord::dot_product(diff, diff);
    const auto expon_scaling = -(exponent0 * exponent1) / (exponent0 + exponent1);

    return std::exp(expon_scaling * norm_sq);
}

}  // anonymous namespace

namespace elec::math
{

/*
    NOTE: this function calculates product of two gaussians with coefficients of 1 each
      - the actual coefficients of the two input gaussians are handled by the functions that
        calculate the matrices
*/
auto gaussian_product(
    const coord::Cartesian3D& centre0,
    const coord::Cartesian3D& centre1,
    double exponent0,
    double exponent1
) -> std::tuple<coord::Cartesian3D, double>
{
    const auto new_centre = (centre0 * exponent0 + centre1 * exponent1) / (exponent0 + exponent1);
    const auto new_coefficient = gaussian_product_coefficient(centre0, centre1, exponent0, exponent1);

    return {new_centre, new_coefficient};
}

auto gaussian_norm(const elec::AngularMomentumNumbers& angular_momenta, double gauss_exponent) -> double
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
