#pragma once

#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/orbitals.hpp"

namespace elec
{

auto nuclear_electron_integral_contraction(
    const AngularMomentumNumbers& angmom_0,
    const AngularMomentumNumbers& angmom_1,
    const coord::Cartesian3D& pos_gauss0,
    const coord::Cartesian3D& pos_gauss1,
    const coord::Cartesian3D& pos_nuclear,
    double exponent0,
    double exponent1,
    double nuclear_charge
) -> double;

}  // namespace elec
