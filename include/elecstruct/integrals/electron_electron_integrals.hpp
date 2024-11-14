#pragma once

#include <cmath>

#include "elecstruct/integrals/expansion_coefficient.hpp"
#include "elecstruct/mathtools/factorial.hpp"
#include "elecstruct/mathtools/misc.hpp"

namespace impl_elec::electron_electron_integrals
{

struct AngularMomenta1D
{
    std::int64_t angmom_0;
    std::int64_t angmom_1;
    std::int64_t angmom_2;
    std::int64_t angmom_3;
};

struct Positions1D
{
    double position_0;
    double position_1;
    double position_2;
    double position_3;
    double position_prod01;
    double position_prod23;
};

struct ElectronElectronIntegralIndices
{
    std::int64_t idx_l_01;
    std::int64_t idx_r_01;
    std::int64_t idx_l_23;
    std::int64_t idx_r_23;
    std::int64_t idx_i;
};

struct GaussianExponentInfo
{
    double exponent_01;
    double exponent_23;
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

inline auto electron_electron_b_factor(
    const ElectronElectronIntegralIndices& indices,
    const AngularMomenta1D& angmoms,
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
    // book: indices.idx_l_01 (i.e. ll in the reference, l' in the book)
    const auto sign = elec::math::neg_1_power(indices.idx_l_01 + indices.idx_i);

    // numerator
    const auto theta01 = electron_electron_theta_factor(
        indices.idx_l_01, angmoms.angmom_0, angmoms.angmom_1, indices.idx_r_01, diff_0, diff_1, info.exponent_01
    );
    const auto theta23 = electron_electron_theta_factor(
        indices.idx_l_23, angmoms.angmom_2, angmoms.angmom_3, indices.idx_r_23, diff_2, diff_3, info.exponent_23
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

}  // namespace elec
