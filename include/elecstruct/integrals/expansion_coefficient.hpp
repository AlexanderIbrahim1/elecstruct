#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <tuple>

#include "elecstruct/mathtools/n_choose_k.hpp"


namespace impl_elec::expansion_coefficient
{

inline auto range_limits_(
    std::uint64_t angmom_j,
    std::uint64_t angmom_l,
    std::uint64_t angmom_m
) -> std::tuple<std::uint64_t, std::uint64_t>
{
    const auto minimum = [&]() {
        if (angmom_j < angmom_m) {
            return std::uint64_t {0};
        } else {
            return angmom_j - angmom_m;
        }
    }();

    const auto maximum = std::max(angmom_j, angmom_l) + 1;

    return {minimum, maximum};
}

}  // namespace impl_elec::expansion_coefficient


namespace elec
{

/*
    The expansion coefficient for the calculation of Gaussian-Type Functions, given in:

    title: Handbook of Computational Quantum Chemistry
    author: David B. Cook
    year: 2005
    published: Dover Publications

    at the bottom on page 219.
*/
inline auto expansion_coefficient(
    std::uint64_t angmom_j,
    std::uint64_t angmom_l,
    std::uint64_t angmom_m,
    double separation0,
    double separation1
) -> double
{
    const auto [minimum, maximum] = impl_elec::expansion_coefficient::range_limits_(angmom_j, angmom_l, angmom_m);


}

}  // namespace elec
