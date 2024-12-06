#pragma once

#include <cstdint>

#include "elecstruct/cartesian3d.hpp"
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
struct DirectedAngularMomentumNumbers
{
    std::int64_t main;
    std::int64_t other0;
    std::int64_t other1;

    DirectedAngularMomentumNumbers(std::int64_t main_, std::int64_t other0_, std::int64_t other1_)
        : main {main_}
        , other0 {other0_}
        , other1 {other1_}
    {}

    explicit DirectedAngularMomentumNumbers(const AngularMomentumNumbers& angmom)
        : main {angmom.x}
        , other0 {angmom.y}
        , other1 {angmom.z}
    {}
};

struct DirectedCartesian3D
{
    double main;
    double other0;
    double other1;

    DirectedCartesian3D(double main_, double other0_, double other1_)
        : main {main_}
        , other0 {other0_}
        , other1 {other1_}
    {}

    explicit DirectedCartesian3D(const coord::Cartesian3D& position)
        : main {position.x}
        , other0 {position.y}
        , other1 {position.z}
    {}
};

template <typename T>
auto left_cyclic_shift(const T& coordinates) -> T
{
    return T {coordinates.other0, coordinates.other1, coordinates.main};
}

template <typename T>
auto right_cyclic_shift(const T& coordinates) -> T
{
    return T {coordinates.other1, coordinates.main, coordinates.other0};
}

auto unnormalized_kinetic_integral_1d(
    const DirectedAngularMomentumNumbers& angmom_a,
    const DirectedAngularMomentumNumbers& angmom_b,
    const DirectedCartesian3D& position_a,
    const DirectedCartesian3D& position_b,
    const DirectedCartesian3D& position_centre,
    double exponent_a,
    double exponent_b,
    double centre_coefficient
) -> double;

auto kinetic_integral_contraction(
    const AngularMomentumNumbers& angmom0,
    const AngularMomentumNumbers& angmom1,
    const coord::Cartesian3D& position0,
    const coord::Cartesian3D& position1,
    double exponent0,
    double exponent1
) -> double;

}  // namespace elec
