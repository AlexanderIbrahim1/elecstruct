#pragma once

#include <array>
#include <cstdint>

#include "elecstruct/integrals/integrals.hpp"
#include "elecstruct/mathtools/factorial.hpp"
#include "elecstruct/mathtools/n_choose_k.hpp"

namespace elec
{

struct OverlapIntegralGaussianInfo1D
{
    std::uint64_t angular_momentum;
    double exponent;
    double centre;
};

inline auto overlap_integral_1d(
    const OverlapIntegralGaussianInfo1D& gaussian0,
    const OverlapIntegralGaussianInfo1D& gaussian1,
    double total_centre
) -> double
{
    if (gaussian0.angular_momentum >= elec::math::N_CHOOSE_K_GRID.size0()) {
        throw std::runtime_error {"Encountered an angular momentum beyond the N-choose-K grid size bounds."};
    }

    if (gaussian1.angular_momentum >= elec::math::N_CHOOSE_K_GRID.size1()) {
        throw std::runtime_error {"Encountered an angular momentum beyond the N-choose-K grid size bounds."};
    }

    auto overlap = double {0.0};
    const auto ang_mom0 = gaussian0.angular_momentum;
    const auto ang_mom1 = gaussian1.angular_momentum;

    for (std::uint64_t i0 {0}; i0 < ang_mom0 + 1; ++i0) {
        for (std::uint64_t i1 {0}; i1 < ang_mom1 + 1; ++i1) {
            if ((i0 + i1) % 2 != 0) {
                continue;
            }

            const auto choose_term0 = elec::math::N_CHOOSE_K_GRID.at(ang_mom0, i0);
            const auto choose_term1 = elec::math::N_CHOOSE_K_GRID.at(ang_mom1, i1);
            const auto factorial_term = elec::math::double_factorial(i0 + i1 - 1);

            const auto gauss_1d_coeff = [&]()
            {
                const auto argument = 2.0 * (gaussian0.exponent + gaussian1.exponent);
                const auto power = 0.5 * static_cast<double>(i0 + i1);
                return std::pow(argument, power);
            }();

            const auto gauss0_contrib = std::pow(total_centre - gaussian0.centre, ang_mom0 - i0);
            const auto gauss1_contrib = std::pow(total_centre - gaussian1.centre, ang_mom1 - i1);

            const auto combinatoric_part = static_cast<double>(choose_term0 * choose_term1 * factorial_term);
            const auto gaussian_part = gauss0_contrib * gauss1_contrib / gauss_1d_coeff;

            overlap += (combinatoric_part * gaussian_part);
        }
    }

    return overlap;
}

}  // namespace elec
