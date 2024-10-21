#pragma once

#include <array>
#include <cstdint>

#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/integrals/integrals.hpp"
#include "elecstruct/mathtools/factorial.hpp"
#include "elecstruct/mathtools/n_choose_k.hpp"
#include "elecstruct/orbitals.hpp"

namespace elec
{

struct OverlapIntegralGaussianInfo1D
{
    std::uint64_t angular_momentum;
    double exponent;
    double centre;
};

inline auto unnormalized_overlap_integral_1d(
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

    auto unnormalized_overlap = double {0.0};
    const auto ang_mom0 = gaussian0.angular_momentum;
    const auto ang_mom1 = gaussian1.angular_momentum;

    for (std::uint64_t i0 {0}; i0 < ang_mom0 + 1; ++i0) {
        for (std::uint64_t i1 {0}; i1 < ang_mom1 + 1; ++i1) {
            if ((i0 + i1) % 2 != 0) {
                continue;
            }

            const auto choose_term0 = elec::math::N_CHOOSE_K_GRID.at(ang_mom0, i0);
            const auto choose_term1 = elec::math::N_CHOOSE_K_GRID.at(ang_mom1, i1);
            const auto factorial_term = elec::math::double_factorial_minus_1(i0 + i1);

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

            unnormalized_overlap += (combinatoric_part * gaussian_part);
        }
    }

    return unnormalized_overlap;
}

inline auto overlap_integral_3d_norm(const GaussianInfo& info0, const GaussianInfo& info1) -> double
{
    const auto argument = M_PI / (info0.exponent + info1.exponent);
    return std::sqrt(argument * argument * argument);
}

inline auto overlap_integral(
    const AngularMomentumNumbers& angmom0,
    const AngularMomentumNumbers& angmom1,
    const coord::Cartesian3D& position0,
    const coord::Cartesian3D& position1,
    const GaussianInfo& info0,
    const GaussianInfo& info1
) -> double
{
    [[maybe_unused]] const auto [new_position, new_info] = gaussian_product(position0, position1, info0, info1);

    const auto norm0 = gaussian_norm(angmom0, info0.exponent);
    const auto norm1 = gaussian_norm(angmom1, info1.exponent);
    const auto overlap_norm = overlap_integral_3d_norm(info0, info1);
    const auto total_norm = norm0 * norm1 * overlap_norm;

    // clang-format off
    const auto unorm_overlap_x = unnormalized_overlap_integral_1d(
        {angmom0.x, info0.exponent, position0.x},
        {angmom1.x, info1.exponent, position1.x},
        new_position.x
    );

    const auto unorm_overlap_y = unnormalized_overlap_integral_1d(
        {angmom0.y, info0.exponent, position0.y},
        {angmom1.y, info1.exponent, position1.y},
        new_position.y
    );

    const auto unorm_overlap_z = unnormalized_overlap_integral_1d(
        {angmom0.z, info0.exponent, position0.z},
        {angmom1.z, info1.exponent, position1.z},
        new_position.z
    );
    // clang-format on

    return new_info.coefficient * total_norm * unorm_overlap_x * unorm_overlap_y * unorm_overlap_z;
}

}  // namespace elec
