#pragma once

#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/orbitals.hpp"

namespace elec
{

auto electron_electron_integral_contraction(
    const AngularMomentumNumbers& angmom_0,
    const AngularMomentumNumbers& angmom_1,
    const AngularMomentumNumbers& angmom_2,
    const AngularMomentumNumbers& angmom_3,
    const coord::Cartesian3D& pos_gauss0,
    const coord::Cartesian3D& pos_gauss1,
    const coord::Cartesian3D& pos_gauss2,
    const coord::Cartesian3D& pos_gauss3,
    double exponent0,
    double exponent1,
    double exponent2,
    double exponent3
) -> double;

}  // namespace elec
