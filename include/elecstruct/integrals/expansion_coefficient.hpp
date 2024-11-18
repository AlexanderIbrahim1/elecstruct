#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <tuple>

#include "elecstruct/mathtools/n_choose_k.hpp"


namespace elec
{

namespace impl_elec::expansion_coefficient
{

inline auto range_limits_(
    std::int64_t angmom_j,
    std::int64_t angmom_l,
    std::int64_t angmom_m
) -> std::tuple<std::int64_t, std::int64_t>
{
    const auto minimum = [&]() {
        if (angmom_j < angmom_m) {
            return std::int64_t {0};
        } else {
            return angmom_j - angmom_m;
        }
    }();

    const auto maximum = std::max(angmom_j, angmom_l) + 1;

    return {minimum, maximum};
}

}  // namespace elec::impl_elec::expansion_coefficient


/*
    The expansion coefficient for the calculation of Gaussian-Type Functions, given in:

    title: Handbook of Computational Quantum Chemistry
    author: David B. Cook
    year: 2005
    published: Dover Publications

    at the bottom on page 219.
*/
inline auto expansion_coefficient(
    std::int64_t angmom_j,
    std::int64_t angmom_l,
    std::int64_t angmom_m,
    double separation0,
    double separation1
) -> double
{
    namespace ieec = impl_elec::expansion_coefficient;

    auto result = double {0.0};

    const auto [minimum, maximum] = ieec::range_limits_(angmom_j, angmom_l, angmom_m);
    for (std::int64_t angmom_k {minimum}; angmom_k < maximum; ++angmom_k) {
        const auto binom0 = elec::math::N_CHOOSE_K_GRID.at(angmom_l, angmom_k);
        const auto coeff0 = std::pow(separation0, angmom_l - angmom_k);

        const auto binom1 = elec::math::N_CHOOSE_K_GRID.at(angmom_m, angmom_j - angmom_k);
        const auto coeff1 = std::pow(separation1, angmom_m - angmom_j + angmom_k);

        const auto contribution = static_cast<double>(binom0 * binom1) * coeff0 * coeff1;

        result += contribution;
    }

    return result;
}

}  // namespace elec
