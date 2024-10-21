#pragma once

#include <array>
#include <cstdint>

#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/integrals/integrals.hpp"
#include "elecstruct/integrals/overlap_integrals.hpp"
#include "elecstruct/orbitals.hpp"

namespace impl_elec::unorm_kinetic_integral
{

inline auto term_minusa_minusb(
    const elec::AngularMomentumKineticIntegral& angmom_a,
    const elec::AngularMomentumKineticIntegral& angmom_b,
    const elec::CartesianKineticIntegral& position_a,
    const elec::CartesianKineticIntegral& position_b,
    const elec::CartesianKineticIntegral& position_centre,
    double exponent_a,
    double exponent_b
) -> double
{
    if (angmom_a.main == 0 || angmom_b.main == 0) {
        return 0.0;
    }
    else {
        const auto coeff = 0.5 * static_cast<double>(angmom_a.main * angmom_b.main);
        const auto unorm_overlap = elec::unnormalized_overlap_integral_1d(
            {angmom_a.main - 1, exponent_a, position_a.main},
            {angmom_b.main - 1, exponent_b, position_b.main},
            position_centre.main
        );

        return coeff * unorm_overlap;
    }
}

inline auto term_plusa_minusb(
    const elec::AngularMomentumKineticIntegral& angmom_a,
    const elec::AngularMomentumKineticIntegral& angmom_b,
    const elec::CartesianKineticIntegral& position_a,
    const elec::CartesianKineticIntegral& position_b,
    const elec::CartesianKineticIntegral& position_centre,
    double exponent_a,
    double exponent_b
) -> double
{
    if (angmom_b.main == 0) {
        return 0.0;
    }
    else {
        const auto coeff = -1.0 * exponent_a * static_cast<double>(angmom_b.main);
        const auto unorm_overlap = elec::unnormalized_overlap_integral_1d(
            {angmom_a.main + 1, exponent_a, position_a.main},
            {angmom_b.main - 1, exponent_b, position_b.main},
            position_centre.main
        );

        return coeff * unorm_overlap;
    }
}

inline auto term_minusa_plusb(
    const elec::AngularMomentumKineticIntegral& angmom_a,
    const elec::AngularMomentumKineticIntegral& angmom_b,
    const elec::CartesianKineticIntegral& position_a,
    const elec::CartesianKineticIntegral& position_b,
    const elec::CartesianKineticIntegral& position_centre,
    double exponent_a,
    double exponent_b
) -> double
{
    if (angmom_a.main == 0) {
        return 0.0;
    }
    else {
        const auto coeff = -1.0 * static_cast<double>(angmom_a.main) * exponent_b;
        const auto unorm_overlap = elec::unnormalized_overlap_integral_1d(
            {angmom_a.main - 1, exponent_a, position_a.main},
            {angmom_b.main + 1, exponent_b, position_b.main},
            position_centre.main
        );

        return coeff * unorm_overlap;
    }
}

inline auto term_plusa_plusb(
    const elec::AngularMomentumKineticIntegral& angmom_a,
    const elec::AngularMomentumKineticIntegral& angmom_b,
    const elec::CartesianKineticIntegral& position_a,
    const elec::CartesianKineticIntegral& position_b,
    const elec::CartesianKineticIntegral& position_centre,
    double exponent_a,
    double exponent_b
) -> double
{
    const auto coeff = 2.0 * exponent_a * exponent_b;
    const auto unorm_overlap = elec::unnormalized_overlap_integral_1d(
        {angmom_a.main + 1, exponent_a, position_a.main},
        {angmom_b.main + 1, exponent_b, position_b.main},
        position_centre.main
    );

    return coeff * unorm_overlap;
}

inline auto term_second_direction(
    const elec::AngularMomentumKineticIntegral& angmom_a,
    const elec::AngularMomentumKineticIntegral& angmom_b,
    const elec::CartesianKineticIntegral& position_a,
    const elec::CartesianKineticIntegral& position_b,
    const elec::CartesianKineticIntegral& position_centre,
    double exponent_a,
    double exponent_b
) -> double
{
    return elec::unnormalized_overlap_integral_1d(
        {angmom_a.other0, exponent_a, position_a.other0},
        {angmom_b.other0, exponent_b, position_b.other0},
        position_centre.other0
    );
}

inline auto term_third_direction(
    const elec::AngularMomentumKineticIntegral& angmom_a,
    const elec::AngularMomentumKineticIntegral& angmom_b,
    const elec::CartesianKineticIntegral& position_a,
    const elec::CartesianKineticIntegral& position_b,
    const elec::CartesianKineticIntegral& position_centre,
    double exponent_a,
    double exponent_b
) -> double
{
    return elec::unnormalized_overlap_integral_1d(
        {angmom_a.other1, exponent_a, position_a.other1},
        {angmom_b.other1, exponent_b, position_b.other1},
        position_centre.other1
    );
}

}  // namespace impl_elec

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
    const CartesianKineticIntegral& position_centre,
    double exponent_a,
    double exponent_b,
    double centre_coefficient
) -> double
{
    namespace uki = impl_elec::unorm_kinetic_integral;

    // clang-format off
    const auto term_ma_mb = uki::term_minusa_minusb(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_pa_mb = uki::term_plusa_minusb(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_ma_pb = uki::term_minusa_plusb(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_pa_pb = uki::term_plusa_plusb(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_dir2 = uki::term_second_direction(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_dir3 = uki::term_third_direction(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    // clang-format on

    double kinetic_coefficient = centre_coefficient * overlap_integral_3d_norm(exponent_a, exponent_b);

    return kinetic_coefficient * term_ma_mb * term_pa_mb * term_ma_pb * term_pa_pb * term_dir2 * term_dir3;
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
