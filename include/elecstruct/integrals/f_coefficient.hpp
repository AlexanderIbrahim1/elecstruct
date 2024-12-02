#pragma once

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
auto f_coefficient(
    std::int64_t angmom_j,
    std::int64_t angmom_l,
    std::int64_t angmom_m,
    double separation0,
    double separation1
) -> double;

}  // namespace elec
