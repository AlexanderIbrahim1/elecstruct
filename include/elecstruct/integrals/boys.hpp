#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>

#include "elecstruct/mathtools/factorial.hpp"

namespace elec
{

namespace impl_elec::boys
{

inline auto boys_small(double x, std::int64_t order, std::int64_t n_max_terms) -> double
{
    auto result = double {0.0};

    for (std::int64_t k {0}; k < n_max_terms; ++k) {
        const auto numerator = std::pow(-x, k);
        const auto denominator = (2 * order + 2 * k + 1) * elec::math::factorial(k);

        result += numerator / static_cast<double>(denominator);
    }

    return result;
}

inline auto boys_large(double x, std::int64_t order) -> double
{
    const auto term_fact = static_cast<double>(elec::math::double_factorial(2 * order - 1));
    const auto term_expon = std::pow(2.0, order + 1);
    const auto term_sqrt = std::sqrt(M_PI / std::pow(x, 2 * order + 1));

    return term_fact * term_sqrt / term_expon;
}

}  // namespace elec::impl_elec::boys

// NOTE: too large a number, and it won't fit in a 64-bit signed integer
constexpr auto N_TERMS_BOYS_SERIES_EXPANSION = std::int64_t {10};
constexpr auto N_MAX_TERMS_BOYS_SMALL = std::int64_t {11};

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

    return coeff * summ;
}

inline auto boys_mix_small_large(double x, std::int64_t order, std::int64_t n_small_terms) -> double
{
    if (n_small_terms < 0 || n_small_terms % 2 == 0) {
        throw std::runtime_error {
            "The number of terms in the small boys approximation must be a positive odd integer."
        };
    }

    const auto small = elec::impl_elec::boys::boys_small(x, order, n_small_terms);
    const auto large = elec::impl_elec::boys::boys_large(x, order);

    return std::min(small, large);
}

}  // namespace elec
