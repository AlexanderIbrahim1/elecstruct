#include <cmath>
#include <cstdint>

#include "elecstruct/basis/basis.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/geometry.hpp"
#include "elecstruct/integrals/boys.hpp"
#include "elecstruct/integrals/f_coefficient.hpp"
#include "elecstruct/integrals/nuclear_electron_index_iterator.hpp"
#include "elecstruct/mathtools/gaussian.hpp"
#include "elecstruct/mathtools/misc.hpp"
#include "elecstruct/orbitals.hpp"

#include "elecstruct/integrals/nuclear_electron_integrals.hpp"

namespace
{

struct AngularMomenta1D
{
    std::int64_t angmom_0;
    std::int64_t angmom_1;
};

struct Positions1D
{
    double position_0;
    double position_1;
    double position_nuclear;
    double position_product;
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
auto nuclear_a_factor(
    std::int64_t idx_l,
    std::int64_t idx_r,
    std::int64_t idx_i,
    const AngularMomenta1D& angmoms,
    const Positions1D& positions,
    double epsilon
) -> double
{
    const auto diff_0 = positions.position_product - positions.position_0;
    const auto diff_1 = positions.position_product - positions.position_1;
    const auto diff_n = positions.position_product - positions.position_nuclear;
    const auto idx_n = idx_l - 2 * (idx_r + idx_i);

    const auto sign = elec::math::neg_1_power(idx_l + idx_i);
    const auto expansion = elec::f_coefficient(idx_l, angmoms.angmom_0, angmoms.angmom_1, diff_0, diff_1);
    const auto epsilon_exponent = std::pow(epsilon, idx_r + idx_i);
    const auto diff_n_exponent = std::pow(diff_n, idx_n);
    const auto fact_l = elec::math::factorial(idx_l);

    const auto fact_r = elec::math::factorial(idx_r);
    const auto fact_i = elec::math::factorial(idx_i);
    const auto fact_diff_n = elec::math::factorial(idx_n);

    const auto numerator = static_cast<double>(sign * fact_l) * expansion * epsilon_exponent * diff_n_exponent;
    const auto denominator = static_cast<double>(fact_r * fact_i * fact_diff_n);

    return numerator / denominator;
}

}  // anonymous namespace

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
) -> double
{
    const auto [pos_product, coeff_product] =
        elec::math::gaussian_product(pos_gauss0, pos_gauss1, exponent0, exponent1);

    const auto norm0 = elec::math::gaussian_norm(angmom_0, exponent0);
    const auto norm1 = elec::math::gaussian_norm(angmom_1, exponent1);

    const auto g_value = exponent0 + exponent1;
    const auto epsilon = 0.25 / g_value;
    const auto boys_arg = g_value * coord::norm_squared(pos_product - pos_nuclear);

    const auto angmoms_x = AngularMomenta1D {angmom_0.x, angmom_1.x};
    const auto angmoms_y = AngularMomenta1D {angmom_0.y, angmom_1.y};
    const auto angmoms_z = AngularMomenta1D {angmom_0.z, angmom_1.z};
    const auto positions_x = Positions1D {pos_gauss0.x, pos_gauss1.x, pos_nuclear.x, pos_product.x};
    const auto positions_y = Positions1D {pos_gauss0.y, pos_gauss1.y, pos_nuclear.y, pos_product.y};
    const auto positions_z = Positions1D {pos_gauss0.z, pos_gauss1.z, pos_nuclear.z, pos_product.z};

    auto integral = double {0.0};

    // clang-format off
    for (const auto [idx_l, idx_r, idx_i] : NuclearElectronIndexGenerator(angmom_0.x, angmom_1.x)) {
        const auto a_factor_x = nuclear_a_factor(idx_l, idx_r, idx_i, angmoms_x, positions_x, epsilon);

        for (const auto [idx_m, idx_s, idx_j] : NuclearElectronIndexGenerator(angmom_0.y, angmom_1.y)) {
            const auto a_factor_y = nuclear_a_factor(idx_m, idx_s, idx_j, angmoms_y, positions_y, epsilon);

            for (const auto [idx_n, idx_t, idx_k] : NuclearElectronIndexGenerator(angmom_0.z, angmom_1.z)) {
                const auto a_factor_z = nuclear_a_factor(idx_n, idx_t, idx_k, angmoms_z, positions_z, epsilon);

                const auto idx_boys = idx_l + idx_m + idx_n - 2 * (idx_r + idx_s + idx_t) - (idx_i + idx_j + idx_k);
                const auto boys_factor = boys_beylkin_sharma(boys_arg, static_cast<std::size_t>(idx_boys));

                const auto contribution = a_factor_x * a_factor_y * a_factor_z * boys_factor;
                integral += contribution;
            }
        }
    }
    // clang-format on

    return -(2.0 * M_PI / g_value) * coeff_product * nuclear_charge * norm0 * norm1 * integral;
}

}  // namespace elec
