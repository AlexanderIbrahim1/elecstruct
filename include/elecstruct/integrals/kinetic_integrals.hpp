#pragma once

#include <array>
#include <cstdint>

#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/integrals/integrals.hpp"
#include "elecstruct/integrals/overlap_integrals.hpp"
#include "elecstruct/orbitals.hpp"

namespace elec
{

/*
    The calculation for the component of the kinetic integral along each axis is complicated,
    and involves mixing the different x-, y-, and z-axis compoments of the angular momenta and
    the cartesian positions.

    Instead of just having an (x, y, z) split, we have a split into:
      - the direction of interest (main)
      - the "second" direction (other0)
      - the "third" direction (other1)
*/

struct AngularMomentumKineticIntegral
{
    std::uint64_t main;
    std::uint64_t other0;
    std::uint64_t other1;

};

struct CartesianKineticIntegral
{
    double main;
    double other0;
    double other1;
};

inline auto unnormalized_kinetic_integral_1d(
    const AngularMomentumKineticIntegral& angmom_a,
    const AngularMomentumKineticIntegral& angmom_b,
    const CartesianKineticIntegral& position_a,
    const CartesianKineticIntegral& position_b,
    const CartesianKineticIntegral& position_centre
) -> double
{

}

inline auto kinetic_integral(
    const AngularMomentumNumbers& angmom0,
    const AngularMomentumNumbers& angmom1,
    const coord::Cartesian3D& position0,
    const coord::Cartesian3D& position1,
    const GaussianInfo& info0,
    const GaussianInfo& info1
) -> double
{
}

}  // namespace elec
