#include <array>
#include <cstdint>

#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/mathtools/gaussian.hpp"
#include "elecstruct/integrals/overlap_integrals.hpp"
#include "elecstruct/orbitals.hpp"

#include "elecstruct/integrals/kinetic_integrals.hpp"

namespace
{

auto term_minusa_minusb(
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

auto term_plusa_minusb(
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

auto term_minusa_plusb(
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

auto term_plusa_plusb(
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

auto term_direction_other0(
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

auto term_direction_other1(
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

}  // anonymous namespace

namespace elec
{

auto unnormalized_kinetic_integral_1d(
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
    // clang-format off
    const auto kinetic_coefficient = centre_coefficient * overlap_integral_3d_norm(exponent_a, exponent_b);
    const auto term_other0 = term_direction_other0(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_other1 = term_direction_other1(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_ma_mb = term_minusa_minusb(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_pa_mb = term_plusa_minusb(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_ma_pb = term_minusa_plusb(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    const auto term_pa_pb = term_plusa_plusb(angmom_a, angmom_b, position_a, position_b, position_centre, exponent_a, exponent_b);
    // clang-format on

    return kinetic_coefficient * term_other0 * term_other1 * (term_ma_mb + term_pa_mb + term_ma_pb + term_pa_pb);
}

auto kinetic_integral_contraction(
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
