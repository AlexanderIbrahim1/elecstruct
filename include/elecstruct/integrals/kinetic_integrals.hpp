#pragma once

#include <array>
#include <cstdint>

#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/mathtools/gaussian.hpp"
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
auto left_cyclic_shift(const T& coordinates) -> T {
    return T {coordinates.other0, coordinates.other1, coordinates.main};
}

template <typename T>
auto right_cyclic_shift(const T& coordinates) -> T {
    return T {coordinates.other1, coordinates.main, coordinates.other0};
}

namespace impl_elec::unorm_kinetic_integral
{

inline auto term_minusa_minusb(
    const elec::DirectedAngularMomentumNumbers& angmom_a,
    const elec::DirectedAngularMomentumNumbers& angmom_b,
    const elec::DirectedCartesian3D& position_a,
    const elec::DirectedCartesian3D& position_b,
    const elec::DirectedCartesian3D& position_centre,
    double exponent_a,
    double exponent_b
) -> double
{
    const auto coeff = 0.5 * static_cast<double>(angmom_a.main * angmom_b.main);
    const auto unorm_overlap = elec::unnormalized_overlap_integral_1d(
        {angmom_a.main - 1, exponent_a, position_a.main},
        {angmom_b.main - 1, exponent_b, position_b.main},
        position_centre.main
    );

    return coeff * unorm_overlap;
}

inline auto term_plusa_minusb(
    const elec::DirectedAngularMomentumNumbers& angmom_a,
    const elec::DirectedAngularMomentumNumbers& angmom_b,
    const elec::DirectedCartesian3D& position_a,
    const elec::DirectedCartesian3D& position_b,
    const elec::DirectedCartesian3D& position_centre,
    double exponent_a,
    double exponent_b
) -> double
{
    const auto coeff = -1.0 * exponent_a * static_cast<double>(angmom_b.main);
    const auto unorm_overlap = elec::unnormalized_overlap_integral_1d(
        {angmom_a.main + 1, exponent_a, position_a.main},
        {angmom_b.main - 1, exponent_b, position_b.main},
        position_centre.main
    );

    return coeff * unorm_overlap;
}

inline auto term_minusa_plusb(
    const elec::DirectedAngularMomentumNumbers& angmom_a,
    const elec::DirectedAngularMomentumNumbers& angmom_b,
    const elec::DirectedCartesian3D& position_a,
    const elec::DirectedCartesian3D& position_b,
    const elec::DirectedCartesian3D& position_centre,
    double exponent_a,
    double exponent_b
) -> double
{
    const auto coeff = -1.0 * static_cast<double>(angmom_a.main) * exponent_b;
    const auto unorm_overlap = elec::unnormalized_overlap_integral_1d(
        {angmom_a.main - 1, exponent_a, position_a.main},
        {angmom_b.main + 1, exponent_b, position_b.main},
        position_centre.main
    );

    return coeff * unorm_overlap;
}

inline auto term_plusa_plusb(
    const elec::DirectedAngularMomentumNumbers& angmom_a,
    const elec::DirectedAngularMomentumNumbers& angmom_b,
    const elec::DirectedCartesian3D& position_a,
    const elec::DirectedCartesian3D& position_b,
    const elec::DirectedCartesian3D& position_centre,
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

inline auto term_direction_other0(
    const elec::DirectedAngularMomentumNumbers& angmom_a,
    const elec::DirectedAngularMomentumNumbers& angmom_b,
    const elec::DirectedCartesian3D& position_a,
    const elec::DirectedCartesian3D& position_b,
    const elec::DirectedCartesian3D& position_centre,
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

inline auto term_direction_other1(
    const elec::DirectedAngularMomentumNumbers& angmom_a,
    const elec::DirectedAngularMomentumNumbers& angmom_b,
    const elec::DirectedCartesian3D& position_a,
    const elec::DirectedCartesian3D& position_b,
    const elec::DirectedCartesian3D& position_centre,
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

}  // namespace elec::impl_elec::unorm_kinetic_integral


inline auto unnormalized_kinetic_integral_1d(
    const DirectedAngularMomentumNumbers& angmom_a,
    const DirectedAngularMomentumNumbers& angmom_b,
    const DirectedCartesian3D& position_a,
    const DirectedCartesian3D& position_b,
    const DirectedCartesian3D& position_centre,
    double exponent_a,
    double exponent_b,
    double centre_coefficient
) -> double
{
    namespace uki = impl_elec::unorm_kinetic_integral;

    // clang-format off
    const auto kinetic_coefficient = centre_coefficient * overlap_integral_3d_norm(exponent_a, exponent_b);
    const auto term_other0 = uki::term_direction_other0(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_other1 = uki::term_direction_other1(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_ma_mb = uki::term_minusa_minusb(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_pa_mb = uki::term_plusa_minusb(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_ma_pb = uki::term_minusa_plusb(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_pa_pb = uki::term_plusa_plusb(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    // clang-format on

    return kinetic_coefficient * term_other0 * term_other1 * (term_ma_mb + term_pa_mb + term_ma_pb + term_pa_pb);
}

inline auto kinetic_integral_contraction(
    const AngularMomentumNumbers& angmom0,
    const AngularMomentumNumbers& angmom1,
    const coord::Cartesian3D& position0,
    const coord::Cartesian3D& position1,
    double exponent0,
    double exponent1
) -> double
{
    const auto [pos_product, coeff_product] = elec::math::gaussian_product(position0, position1, exponent0, exponent1);

    const auto angmom0_integ_x = DirectedAngularMomentumNumbers {angmom0};
    const auto angmom1_integ_x = DirectedAngularMomentumNumbers {angmom1};
    const auto original0_integ_x = DirectedCartesian3D {position0};
    const auto original1_integ_x = DirectedCartesian3D {position1};
    const auto product_integ_x = DirectedCartesian3D {pos_product};

    const auto kinetic_x = unnormalized_kinetic_integral_1d(
        angmom0_integ_x,
        angmom1_integ_x,
        original0_integ_x,
        original1_integ_x,
        product_integ_x,
        exponent0,
        exponent1,
        coeff_product
    );

    const auto kinetic_y = unnormalized_kinetic_integral_1d(
        left_cyclic_shift(angmom0_integ_x),
        left_cyclic_shift(angmom1_integ_x),
        left_cyclic_shift(original0_integ_x),
        left_cyclic_shift(original1_integ_x),
        left_cyclic_shift(product_integ_x),
        exponent0,
        exponent1,
        coeff_product
    );

    const auto kinetic_z = unnormalized_kinetic_integral_1d(
        right_cyclic_shift(angmom0_integ_x),
        right_cyclic_shift(angmom1_integ_x),
        right_cyclic_shift(original0_integ_x),
        right_cyclic_shift(original1_integ_x),
        right_cyclic_shift(product_integ_x),
        exponent0,
        exponent1,
        coeff_product
    );

    const auto norm0 = elec::math::gaussian_norm(angmom0, exponent0);
    const auto norm1 = elec::math::gaussian_norm(angmom1, exponent1);

    return norm0 * norm1 * (kinetic_x + kinetic_y + kinetic_z);
}

}  // namespace elec
