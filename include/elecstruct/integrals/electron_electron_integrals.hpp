#pragma once

#include <cmath>

#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/integrals/boys.hpp"
#include "elecstruct/integrals/electron_electron_index_iterator.hpp"
#include "elecstruct/integrals/expansion_coefficient.hpp"
#include "elecstruct/integrals/integrals.hpp"
#include "elecstruct/mathtools/factorial.hpp"
#include "elecstruct/mathtools/misc.hpp"
#include "elecstruct/orbitals.hpp"

namespace elec
{

namespace impl_elec::electron_electron_integrals
{

struct Positions1D
{
    double position_0;
    double position_1;
    double position_2;
    double position_3;
    double position_prod01;
    double position_prod23;
};

struct GaussianExponentInfo
{
    double g_value_01;
    double g_value_23;
    double delta;
};

/*
    This function calculates the "theta-factor" that repeatedly appears in the calculation of
    the electron-electron repulsion integrals.

    title: Handbook of Computational Quantum Chemistry
    author: David B. Cook
    year: 2005
    published: Dover Publications

    on page 229.

    TODO: there might be some numerical instability risks here with the products of factorials?
*/
inline auto electron_electron_theta_factor(
    std::int64_t idx_ltot,
    std::int64_t idx_l0,
    std::int64_t idx_l1,
    std::int64_t idx_r,
    double separation0,
    double separation1,
    double gauss_exponent
) -> double
{
    const auto f_factor = elec::expansion_coefficient(idx_ltot, idx_l0, idx_l1, separation0, separation1);
    const auto ltot_fact = static_cast<double>(elec::math::factorial(idx_ltot));
    const auto r_fact = static_cast<double>(elec::math::factorial(idx_r));
    const auto ltot_r_fact = static_cast<double>(elec::math::factorial(idx_ltot - 2 * idx_r));
    const auto expon = std::pow(gauss_exponent, idx_r - idx_ltot);

    return f_factor * ltot_fact * expon / (r_fact * ltot_r_fact);
}

// NOTE: according to the reference code, the book has the wrong formula
inline auto electron_electron_b_factor(
    const elec::ElectronElectronIntegralIndices& indices,
    const elec::AngularMomenta1D& angmoms,
    const Positions1D& positions,
    const GaussianExponentInfo& info
) -> double
{
    const auto diff_0 = positions.position_prod01 - positions.position_0;
    const auto diff_1 = positions.position_prod01 - positions.position_1;
    const auto diff_2 = positions.position_prod23 - positions.position_2;
    const auto diff_3 = positions.position_prod23 - positions.position_3;
    const auto diff_prod = positions.position_prod01 - positions.position_prod23;

    const auto idx_k = indices.idx_l_01 + indices.idx_l_23 - 2 * (indices.idx_r_01 + indices.idx_r_23);

    // NOTE: code and book have a different variable here:
    // code: indices.idx_l_01 (i.e. l in the reference,  l in the book)
    // book: indices.idx_l_23 (i.e. ll in the reference, l' in the book)
    const auto sign = elec::math::neg_1_power(indices.idx_l_01 + indices.idx_i);

    // numerator
    const auto theta01 = electron_electron_theta_factor(
        indices.idx_l_01, angmoms.angmom_0, angmoms.angmom_1, indices.idx_r_01, diff_0, diff_1, info.g_value_01
    );
    const auto theta23 = electron_electron_theta_factor(
        indices.idx_l_23, angmoms.angmom_2, angmoms.angmom_3, indices.idx_r_23, diff_2, diff_3, info.g_value_23
    );
    const auto k_factorial = static_cast<double>(elec::math::factorial(idx_k));
    const auto expon_position = std::pow(diff_prod, idx_k - 2 * indices.idx_i);

    // denominator
    const auto pow_2_coeff = std::pow(2.0, -(indices.idx_l_01 + indices.idx_l_23 + indices.idx_i));
    const auto delta_factor = std::pow(2.0 * info.delta, idx_k - indices.idx_i);
    const auto i_factorial = static_cast<double>(elec::math::factorial(indices.idx_i));
    const auto k2i_factorial = static_cast<double>(elec::math::factorial(idx_k - 2 * indices.idx_i));

    const auto numerator = theta01 * theta23 * k_factorial * expon_position;
    const auto denominator = pow_2_coeff * delta_factor * i_factorial * k2i_factorial;

    return sign * numerator / denominator;
}

inline auto boys_index(
    const elec::ElectronElectronIntegralIndices& indices_x,
    const elec::ElectronElectronIntegralIndices& indices_y,
    const elec::ElectronElectronIntegralIndices& indices_z
) noexcept -> std::int64_t
{
    // clang-format off
    const auto idx_l_sum = indices_x.idx_l_01 + indices_x.idx_l_23
                         + indices_y.idx_l_01 + indices_y.idx_l_23
                         + indices_z.idx_l_01 + indices_z.idx_l_23;

    const auto idx_r_sum = indices_x.idx_r_01 + indices_x.idx_r_23
                         + indices_y.idx_r_01 + indices_y.idx_r_23
                         + indices_z.idx_r_01 + indices_z.idx_r_23;

    const auto idx_i_sum = indices_x.idx_i + indices_y.idx_i + indices_z.idx_i;
    // clang-format on

    return idx_l_sum - 2 * idx_r_sum - idx_i_sum;
}

}  // namespace elec::impl_elec::electron_electron_integrals


inline auto electron_electron_integral(
    const AngularMomentumNumbers& angmom_0,
    const AngularMomentumNumbers& angmom_1,
    const AngularMomentumNumbers& angmom_2,
    const AngularMomentumNumbers& angmom_3,
    const coord::Cartesian3D& pos_gauss0,
    const coord::Cartesian3D& pos_gauss1,
    const coord::Cartesian3D& pos_gauss2,
    const coord::Cartesian3D& pos_gauss3,
    const GaussianInfo& info0,
    const GaussianInfo& info1,
    const GaussianInfo& info2,
    const GaussianInfo& info3
) -> double
{
    namespace eli = impl_elec::electron_electron_integrals;

    const auto [pos_product_01, info_product_01] = gaussian_product(pos_gauss0, pos_gauss1, info0, info1);
    const auto [pos_product_23, info_product_23] = gaussian_product(pos_gauss2, pos_gauss3, info2, info3);

    const auto norm0 = gaussian_norm(angmom_0, info0.exponent);
    const auto norm1 = gaussian_norm(angmom_1, info1.exponent);
    const auto norm2 = gaussian_norm(angmom_2, info2.exponent);
    const auto norm3 = gaussian_norm(angmom_3, info3.exponent);

    const auto g_value_01 = info0.exponent + info1.exponent;
    const auto g_value_23 = info2.exponent + info3.exponent;
    const auto delta = 0.25 * (1.0 / g_value_01 + 1.0 / g_value_23);
    const auto expon_info = eli::GaussianExponentInfo {g_value_01, g_value_23, delta};

    // clang-format off
    const auto angmoms_x = AngularMomenta1D {angmom_0.x, angmom_1.x, angmom_2.x, angmom_3.x};
    const auto angmoms_y = AngularMomenta1D {angmom_0.y, angmom_1.y, angmom_2.y, angmom_3.y};
    const auto angmoms_z = AngularMomenta1D {angmom_0.z, angmom_1.z, angmom_2.z, angmom_3.z};
    const auto positions_x = eli::Positions1D {pos_gauss0.x, pos_gauss1.x, pos_gauss2.x, pos_gauss3.x, pos_product_01.x, pos_product_23.x};
    const auto positions_y = eli::Positions1D {pos_gauss0.y, pos_gauss1.y, pos_gauss2.y, pos_gauss3.y, pos_product_01.y, pos_product_23.y};
    const auto positions_z = eli::Positions1D {pos_gauss0.z, pos_gauss1.z, pos_gauss2.z, pos_gauss3.z, pos_product_01.z, pos_product_23.z};
    // clang-format on

    // clang-format off
    auto integral = double {0.0};

    for (const auto indices_x : ElectronElectronIndexGenerator {angmoms_x}) {
        const auto b_factor_x = eli::electron_electron_b_factor(indices_x, angmoms_x, positions_x, expon_info);

        for (const auto indices_y : ElectronElectronIndexGenerator {angmoms_y}) {
            const auto b_factor_y = eli::electron_electron_b_factor(indices_y, angmoms_y, positions_y, expon_info);

            for (const auto indices_z : ElectronElectronIndexGenerator {angmoms_z}) {
                const auto b_factor_z = eli::electron_electron_b_factor(indices_z, angmoms_z, positions_z, expon_info);

                const auto idx_boys = eli::boys_index(indices_x, indices_y, indices_z);
                const auto boys_arg = 0.25 * coord::norm_squared(pos_product_01 - pos_product_23) / delta;
                const auto boys_factor = boys_function_via_series_expansion(boys_arg, idx_boys);
                
                const auto contribution = b_factor_x * b_factor_y * b_factor_z * boys_factor;
                integral += contribution;
            }
        }
    }
    // clang-format on

    const auto norm_tot = norm0 * norm1 * norm2 * norm3;
    const auto coeff_tot = info_product_01.coefficient * info_product_23.coefficient;
    const auto expon_tot = std::pow(M_PI, 5.0/2.0) / std::pow(g_value_01 + g_value_23, 3.0/2.0);

    return 2.0 * integral * coeff_tot * norm_tot * expon_tot;
}

}  // namespace elec
