#pragma once

#include <array>
#include <cstdint>

#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/mathtools/gaussian.hpp"
#include "elecstruct/mathtools/factorial.hpp"
#include "elecstruct/mathtools/n_choose_k.hpp"
#include "elecstruct/orbitals.hpp"

namespace elec
{

struct OverlapIntegralGaussianInfo1D
{
    std::int64_t angular_momentum;
    double exponent;
    double centre;
};

inline auto gauss_product_coeff(double exponent0, double exponent1, std::int64_t angmom0, std::int64_t angmom1) -> double
{
    const auto argument = 2.0 * (exponent0 + exponent1);
    const auto power = 0.5 * static_cast<double>(angmom0 + angmom1);
    return std::pow(argument, power);
}

inline auto unnormalized_overlap_integral_1d(
    const OverlapIntegralGaussianInfo1D& gaussian0,
    const OverlapIntegralGaussianInfo1D& gaussian1,
    double total_centre
) -> double
{
    // neither of the nested loops to calculate the overlap integral will run, thus giving 0
    if (gaussian0.angular_momentum < 0 || gaussian1.angular_momentum < 0) {
        return 0.0;
    }

    if (static_cast<std::size_t>(gaussian0.angular_momentum) >= elec::math::N_CHOOSE_K_GRID.size0()) {
        throw std::runtime_error {"Encountered an angular momentum beyond the N-choose-K grid size bounds."};
    }

    if (static_cast<std::size_t>(gaussian1.angular_momentum) >= elec::math::N_CHOOSE_K_GRID.size1()) {
        throw std::runtime_error {"Encountered an angular momentum beyond the N-choose-K grid size bounds."};
    }

    auto unnormalized_overlap = double {0.0};
    const auto angmom0 = gaussian0.angular_momentum;
    const auto angmom1 = gaussian1.angular_momentum;

    for (std::int64_t i0 {0}; i0 < angmom0 + 1; ++i0) {
        for (std::int64_t i1 {0}; i1 < angmom1 + 1; ++i1) {
            if ((i0 + i1) % 2 != 0) {
                continue;
            }

            const auto choose_term0 = elec::math::N_CHOOSE_K_GRID.at(angmom0, i0);
            const auto choose_term1 = elec::math::N_CHOOSE_K_GRID.at(angmom1, i1);
            const auto factorial_term = elec::math::double_factorial(i0 + i1 - 1);

            const auto gauss0_contrib = std::pow(total_centre - gaussian0.centre, angmom0 - i0);
            const auto gauss1_contrib = std::pow(total_centre - gaussian1.centre, angmom1 - i1);
            const auto gauss_1d_coeff = gauss_product_coeff(gaussian0.exponent, gaussian1.exponent, i0, i1);

            const auto combinatoric_part = static_cast<double>(choose_term0 * choose_term1 * factorial_term);
            const auto gaussian_part = gauss0_contrib * gauss1_contrib / gauss_1d_coeff;

            unnormalized_overlap += (combinatoric_part * gaussian_part);
        }
    }

    return unnormalized_overlap;
}

inline auto overlap_integral_3d_norm(double exponent0, double exponent1) -> double
{
    const auto argument = M_PI / (exponent0 + exponent1);
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
    [[maybe_unused]] const auto [new_position, new_info] = elec::math::gaussian_product(position0, position1, info0, info1);

    const auto norm0 = elec::math::gaussian_norm(angmom0, info0.exponent);
    const auto norm1 = elec::math::gaussian_norm(angmom1, info1.exponent);
    const auto overlap_norm = overlap_integral_3d_norm(info0.exponent, info1.exponent);
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

    // return new_info.coefficient * total_norm * unorm_overlap_x * unorm_overlap_y * unorm_overlap_z;
    return new_info.coefficient * unorm_overlap_x * unorm_overlap_y * unorm_overlap_z * total_norm;
}

}  // namespace elec
