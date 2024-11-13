#pragma once

#include <cstdint>

#include "elecstruct/integrals/expansion_coefficient.hpp"
#include "elecstruct/mathtools/factorial.hpp"

namespace impl_elec::nuclear_electron_integrals
{

inline auto neg_1_power(std::int64_t arg) -> std::int64_t
{
    if (arg % 2 == 0) {
        return 1;
    } else {
        return -1;
    }
}

}  // namespace impl_elec::nuclear_electron_integrals

namespace elec
{

/*
    The expansion coefficient for the calculation of Gaussian-Type Functions, given in:
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
    std::int64_t angmom0,
    std::int64_t angmom1,
    double position_a,
    double position_b,
    double position_c,
    double position_relative,
    double epsilon
) -> double
{
    namespace nei = impl_elec::nuclear_electron_integrals;

    const auto diff_a = position_relative - position_a;
    const auto diff_b = position_relative - position_b;
    const auto diff_c = position_relative - position_c;
    const auto idx_c = idx_l - 2 * (idx_r + idx_i);

    const auto sign = nei::neg_1_power(idx_l + idx_i);
    const auto expansion = expansion_coefficient(idx_l, angmom0, angmom1, diff_a, diff_b);
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

}  // namespace elec
