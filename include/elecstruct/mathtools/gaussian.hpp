#pragma once

#include <tuple>

#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/orbitals.hpp"

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
) -> std::tuple<coord::Cartesian3D, double>;

auto gaussian_norm(const elec::AngularMomentumNumbers& angular_momenta, double gauss_exponent) -> double;

}  // namespace elec::math
