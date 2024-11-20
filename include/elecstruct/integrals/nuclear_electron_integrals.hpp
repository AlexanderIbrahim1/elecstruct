#pragma once

// TODO: remove
#include <iostream>

#include <cmath>
#include <cstdint>

#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/integrals/expansion_coefficient.hpp"
#include "elecstruct/integrals/boys.hpp"
#include "elecstruct/mathtools/gaussian.hpp"
#include "elecstruct/integrals/nuclear_electron_index_iterator.hpp"
#include "elecstruct/mathtools/factorial.hpp"
#include "elecstruct/mathtools/misc.hpp"
#include "elecstruct/orbitals.hpp"

namespace elec
{

namespace impl_elec::nuclear_electron_integrals
{

struct AngularMomenta1D
{
    std::int64_t angmom_0;
    std::int64_t angmom_1;
};

struct Positions1D
{
    double position_a;
    double position_b;
    double position_c;
    double position_relative;
};

/*
    This function calculates the "A-factor" that repeatedly appears in the calculation of
    the nuclear-electron attraction integrals.

    title: Handbook of Computational Quantum Chemistry
    author: David B. Cook
    year: 2005
    published: Dover Publications

    on page 227.

    TODO: there might be some numerical instability risks here with the products of factorials?
*/
inline auto nuclear_a_factor(
    std::int64_t idx_l,
    std::int64_t idx_r,
    std::int64_t idx_i,
    const AngularMomenta1D& angmoms,
    const Positions1D& positions,
    double epsilon
) -> double
{
    const auto diff_a = positions.position_relative - positions.position_a;
    const auto diff_b = positions.position_relative - positions.position_b;
    const auto diff_c = positions.position_relative - positions.position_c;
    const auto idx_c = idx_l - 2 * (idx_r + idx_i);

    const auto sign = elec::math::neg_1_power(idx_l + idx_i);
    const auto expansion = elec::expansion_coefficient(idx_l, angmoms.angmom_0, angmoms.angmom_1, diff_a, diff_b);
    const auto epsilon_exponent = std::pow(epsilon, idx_r + idx_i);
    const auto diff_c_exponent = std::pow(diff_c, idx_c);
    const auto fact_l = elec::math::factorial(idx_l);

    const auto fact_r = elec::math::factorial(idx_r);
    const auto fact_i = elec::math::factorial(idx_i);
    const auto fact_diff_c = elec::math::factorial(idx_c);

    const auto numerator = static_cast<double>(sign * fact_l) * expansion * epsilon_exponent * diff_c_exponent;
    const auto denominator = static_cast<double>(fact_r * fact_i * fact_diff_c);

    return numerator / denominator;
}

}  // namespace elec::impl_elec::nuclear_electron_integrals


// TODO: maybe turn this into a class with a callable member function?
// - the interface is very similar to the overlap and kinetic integrals, except for the nuclear position and charge
// - but those two variables stay the same when calculating a nuclear-electron integral for a given atom
//
// - if we fix those variables, then we can get three different matrices created with the same interface
//   - and we can pass it as a template argument and reduce a ton of code duplication 
inline auto nuclear_electron_integral(
    const AngularMomentumNumbers& angmom_0,
    const AngularMomentumNumbers& angmom_1,
    const coord::Cartesian3D& pos_gauss0,
    const coord::Cartesian3D& pos_gauss1,
    const coord::Cartesian3D& pos_nuclear,
    const GaussianContractionInfo& info0,
    const GaussianContractionInfo& info1,
    double nuclear_charge
) -> double
{
    namespace nui = impl_elec::nuclear_electron_integrals;

    const auto [pos_product, info_product] = elec::math::gaussian_product(pos_gauss0, pos_gauss1, info0, info1);
    const auto norm0 = elec::math::gaussian_norm(angmom_0, info0.exponent_coeff);
    const auto norm1 = elec::math::gaussian_norm(angmom_1, info1.exponent_coeff);

    const auto g_value = info0.exponent_coeff + info1.exponent_coeff;
    const auto epsilon = 0.25 / g_value;
    const auto boys_arg = g_value * coord::norm_squared(pos_product - pos_nuclear);

    const auto angmoms_x = nui::AngularMomenta1D {angmom_0.x, angmom_1.x};
    const auto angmoms_y = nui::AngularMomenta1D {angmom_0.y, angmom_1.y};
    const auto angmoms_z = nui::AngularMomenta1D {angmom_0.z, angmom_1.z};
    const auto positions_x = nui::Positions1D {pos_gauss0.x, pos_gauss1.x, pos_nuclear.x, pos_product.x};
    const auto positions_y = nui::Positions1D {pos_gauss0.y, pos_gauss1.y, pos_nuclear.y, pos_product.y};
    const auto positions_z = nui::Positions1D {pos_gauss0.z, pos_gauss1.z, pos_nuclear.z, pos_product.z};

    auto integral = double {0.0};

    // clang-format off
    for (const auto [idx_l, idx_r, idx_i] : NuclearElectronIndices(angmom_0.x, angmom_1.x)) {
        const auto a_factor_x = nui::nuclear_a_factor(idx_l, idx_r, idx_i, angmoms_x, positions_x, epsilon);

        for (const auto [idx_m, idx_s, idx_j] : NuclearElectronIndices(angmom_0.y, angmom_1.y)) {
            const auto a_factor_y = nui::nuclear_a_factor(idx_m, idx_s, idx_j, angmoms_y, positions_y, epsilon);

            for (const auto [idx_n, idx_t, idx_k] : NuclearElectronIndices(angmom_0.z, angmom_1.z)) {
                const auto a_factor_z = nui::nuclear_a_factor(idx_n, idx_t, idx_k, angmoms_z, positions_z, epsilon);

                const auto idx_boys = idx_l + idx_m + idx_n - 2 * (idx_r + idx_s + idx_t) - (idx_i + idx_j + idx_k);
                const auto boys_factor = boys_function_via_series_expansion(boys_arg, idx_boys);

                std::cout << "a_factor_x = " << a_factor_x << '\n';
                std::cout << "a_factor_y = " << a_factor_y << '\n';
                std::cout << "a_factor_z = " << a_factor_z << '\n';
                std::cout << "boys_factor = " << boys_factor << '\n';
                
                const auto contribution = a_factor_x * a_factor_y * a_factor_z * boys_factor;
                integral += contribution;
            }
        }
    }
    // clang-format on

    std::cout << "norm0 = " << norm0 << '\n';
    std::cout << "norm1 = " << norm1 << '\n';
    std::cout << "nuclear_charge = " << nuclear_charge << '\n';
    std::cout << "integral = " << integral << '\n';
    std::cout << "g_value = " << g_value << '\n';

    return - 2.0 * M_PI * norm0 * norm1 * nuclear_charge * integral / g_value;

    // Vn *= - Zn * Na * Nb * c * 8 * np.pi * epsilon
}

}  // namespace elec
