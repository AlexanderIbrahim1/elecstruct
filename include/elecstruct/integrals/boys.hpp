#pragma once

#include <cmath>
#include <cstdint>

#include "elecstruct/mathtools/factorial.hpp"

namespace elec
{

constexpr auto N_TERMS_BOYS_SERIES_EXPANSION = std::int64_t {50};

/*
    NOTE: there are many much faster ways of calculating the Boys function, and its evaluation is
    a topic of ongoing research. For the time being, we can settle for a (relatively) slow series
    expansion.
*/
inline auto boys_function_via_series_expansion(double x, std::int64_t order) -> double
{
    const auto coeff = std::exp(-x);

    auto summ = double {0.0};
    for (std::int64_t i {0}; i < N_TERMS_BOYS_SERIES_EXPANSION; ++i) {
        const auto numer_binom = static_cast<double>(elec::math::double_factorial(2 * order - 1));
        const auto denom_binom = static_cast<double>(elec::math::double_factorial(2 * order + 2 * i + 1));
        const auto numer_expon = std::pow(2.0 * x, i);
        const auto contribution = numer_binom * numer_expon / denom_binom;

        summ += contribution;
    }

    return summ;
}

}  // namespace elec
