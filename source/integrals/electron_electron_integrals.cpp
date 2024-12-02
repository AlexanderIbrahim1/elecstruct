#include <cmath>

#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/geometry.hpp"
#include "elecstruct/integrals/boys.hpp"
#include "elecstruct/integrals/electron_electron_index_iterator.hpp"
#include "elecstruct/integrals/f_coefficient.hpp"
#include "elecstruct/mathtools/gaussian.hpp"
#include "elecstruct/mathtools/factorial.hpp"
#include "elecstruct/mathtools/misc.hpp"
#include "elecstruct/orbitals.hpp"

namespace
{

struct PositionDifferences1D
{
    double diff_0;
    double diff_1;
    double diff_2;
    double diff_3;
    double diff_prod;

};

auto make_position_differences(
    double pos_gauss0,
    double pos_gauss1,
    double pos_gauss2,
    double pos_gauss3,
    double pos_product_01,
    double pos_product_23
) -> PositionDifferences1D
{
    const auto diff_0 = pos_product_01 - pos_gauss0;
    const auto diff_1 = pos_product_01 - pos_gauss1;
    const auto diff_2 = pos_product_23 - pos_gauss2;
    const auto diff_3 = pos_product_23 - pos_gauss3;
    const auto diff_prod = pos_product_01 - pos_product_23;

    return PositionDifferences1D {diff_0, diff_1, diff_2, diff_3, diff_prod};
}

struct GaussianExponentInfo
{
    double g_value_01;
    double g_value_23;
    double delta;
};

struct BFactorCalculationInfo
{
    double value;
    bool continue_flag;
};

constexpr auto B_FACTOR_SMALL_TOLERANCE = double {1.0e-8};

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
auto electron_electron_theta_factor(
    std::int64_t idx_ltot,
    std::int64_t idx_l0,
    std::int64_t idx_l1,
    std::int64_t idx_r,
    double separation0,
    double separation1,
    double gauss_exponent
) -> double
{
    const auto f_factor = elec::f_coefficient(idx_ltot, idx_l0, idx_l1, separation0, separation1);
    const auto ltot_fact = static_cast<double>(elec::math::factorial(idx_ltot));
    const auto r_fact = static_cast<double>(elec::math::factorial(idx_r));
    const auto ltot_r_fact = static_cast<double>(elec::math::factorial(idx_ltot - 2 * idx_r));
    const auto expon = std::pow(gauss_exponent, idx_r - idx_ltot);

    return f_factor * ltot_fact * expon / (r_fact * ltot_r_fact);
}

// NOTE: according to the reference code, the book has the wrong formula
auto electron_electron_b_factor(
    const elec::ElectronElectronIntegralIndices& indices,
    const elec::AngularMomenta1D& angmoms,
    const PositionDifferences1D& differences,
    const GaussianExponentInfo& info
) -> BFactorCalculationInfo
{
    const auto idx_k = indices.idx_l_01 + indices.idx_l_23 - 2 * (indices.idx_r_01 + indices.idx_r_23);

    // NOTE: code and book have a different variable here:
    // code: indices.idx_l_01 (i.e. l in the reference,  l in the book)
    // book: indices.idx_l_23 (i.e. ll in the reference, l' in the book)

    // numerator
    const auto theta01 = electron_electron_theta_factor(
        indices.idx_l_01, angmoms.angmom_0, angmoms.angmom_1, indices.idx_r_01,
        differences.diff_0, differences.diff_1,
        info.g_value_01
    );

    if (std::fabs(theta01) < B_FACTOR_SMALL_TOLERANCE) {
        return {0.0, true};
    }

    const auto theta23 = electron_electron_theta_factor(
        indices.idx_l_23, angmoms.angmom_2, angmoms.angmom_3, indices.idx_r_23,
        differences.diff_2, differences.diff_3,
        info.g_value_23
    );

    if (std::fabs(theta23) < B_FACTOR_SMALL_TOLERANCE) {
        return {0.0, true};
    }

    const auto expon_position = std::pow(differences.diff_prod, idx_k - 2 * indices.idx_i);

    if (std::fabs(expon_position) < B_FACTOR_SMALL_TOLERANCE) {
        return {0.0, true};
    }

    const auto sign = static_cast<double>(elec::math::neg_1_power(indices.idx_l_01 + indices.idx_i));
    const auto k_factorial = static_cast<double>(elec::math::factorial(idx_k));

    // denominator
    const auto i_factorial = static_cast<double>(elec::math::factorial(indices.idx_i));
    const auto k2i_factorial = static_cast<double>(elec::math::factorial(idx_k - 2 * indices.idx_i));
    const auto delta_factor = std::pow(info.delta, idx_k - indices.idx_i);
    const auto pow2_factor = std::pow(2.0, idx_k + indices.idx_l_01 + indices.idx_l_23);

    const auto numerator = sign * theta01 * theta23 * k_factorial * expon_position;
    const auto denominator = pow2_factor * delta_factor * i_factorial * k2i_factorial;
    const auto b_factor = numerator / denominator;

    return {b_factor, false};
}


auto boys_index(
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

}  // anonymous namespace


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
) -> double
{
    const auto [pos_product_01, coeff_product_01] = elec::math::gaussian_product(pos_gauss0, pos_gauss1, exponent0, exponent1);
    const auto [pos_product_23, coeff_product_23] = elec::math::gaussian_product(pos_gauss2, pos_gauss3, exponent2, exponent3);

    // clang-format off
    const auto differences_x = make_position_differences(pos_gauss0.x, pos_gauss1.x, pos_gauss2.x, pos_gauss3.x, pos_product_01.x, pos_product_23.x);
    const auto differences_y = make_position_differences(pos_gauss0.y, pos_gauss1.y, pos_gauss2.y, pos_gauss3.y, pos_product_01.y, pos_product_23.y);
    const auto differences_z = make_position_differences(pos_gauss0.z, pos_gauss1.z, pos_gauss2.z, pos_gauss3.z, pos_product_01.z, pos_product_23.z);

    const auto norm0 = elec::math::gaussian_norm(angmom_0, exponent0);
    const auto norm1 = elec::math::gaussian_norm(angmom_1, exponent1);
    const auto norm2 = elec::math::gaussian_norm(angmom_2, exponent2);
    const auto norm3 = elec::math::gaussian_norm(angmom_3, exponent3);

    const auto g_value_01 = exponent0 + exponent1;
    const auto g_value_23 = exponent2 + exponent3;
    const auto delta = 0.25 * (1.0 / g_value_01 + 1.0 / g_value_23);
    const auto expon_info = GaussianExponentInfo {g_value_01, g_value_23, delta};

    const auto angmoms_x = AngularMomenta1D {angmom_0.x, angmom_1.x, angmom_2.x, angmom_3.x};
    const auto angmoms_y = AngularMomenta1D {angmom_0.y, angmom_1.y, angmom_2.y, angmom_3.y};
    const auto angmoms_z = AngularMomenta1D {angmom_0.z, angmom_1.z, angmom_2.z, angmom_3.z};

    auto integral = double {0.0};

    for (const auto indices_x : ElectronElectronIndexGenerator {angmoms_x}) {
        const auto [b_factor_x, continue_x_flag] = electron_electron_b_factor(indices_x, angmoms_x, differences_x, expon_info);
        if (continue_x_flag) {
            continue;
        }

        for (const auto indices_y : ElectronElectronIndexGenerator {angmoms_y}) {
            const auto [b_factor_y, continue_y_flag] = electron_electron_b_factor(indices_y, angmoms_y, differences_y, expon_info);
            if (continue_y_flag) {
                continue;
            }

            for (const auto indices_z : ElectronElectronIndexGenerator {angmoms_z}) {
                const auto [b_factor_z, continue_z_flag] = electron_electron_b_factor(indices_z, angmoms_z, differences_z, expon_info);
                if (continue_z_flag) {
                    continue;
                }

                const auto idx_boys = boys_index(indices_x, indices_y, indices_z);
                const auto boys_arg = 0.25 * coord::norm_squared(pos_product_01 - pos_product_23) / delta;
                const auto boys_factor = boys_beylkin_sharma(boys_arg, static_cast<std::size_t>(idx_boys));
                
                const auto contribution = b_factor_x * b_factor_y * b_factor_z * boys_factor;
                integral += contribution;
            }
        }
    }
    // clang-format on

    const auto norm_tot = norm0 * norm1 * norm2 * norm3;
    const auto coeff_tot = coeff_product_01 * coeff_product_23;
    const auto expon_tot = 2.0 * M_PI * M_PI / (g_value_01 * g_value_23) * std::sqrt(M_PI / (g_value_01 + g_value_23));

    return integral * coeff_tot * norm_tot * expon_tot;
}

}  // namespace elec
