#include <array>
#include <cmath>
#include <cstdint>

#include "elecstruct/basis/basis.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/mathtools/gaussian.hpp"
#include "elecstruct/mathtools/misc.hpp"
#include "elecstruct/mathtools/n_choose_k.hpp"
#include "elecstruct/orbitals.hpp"

#include "elecstruct/integrals/overlap_integrals.hpp"

namespace
{

auto gauss_product_coeff(double exponent0, double exponent1, std::int64_t angmom0, std::int64_t angmom1) -> double
{
    const auto argument = 2.0 * (exponent0 + exponent1);
    const auto power = 0.5 * static_cast<double>(angmom0 + angmom1);
    return std::pow(argument, power);
}

}  // anonymous namespace

namespace elec
{

auto overlap_integral_3d_norm(double exponent0, double exponent1) -> double
{
    const auto argument = M_PI / (exponent0 + exponent1);
    return std::sqrt(argument * argument * argument);
}

auto unnormalized_overlap_integral_1d(
    const OverlapIntegralGaussianContractionInfo1D& gaussian0,
    const OverlapIntegralGaussianContractionInfo1D& gaussian1,
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
            const auto gauss_1d_coeff = gauss_product_coeff(gaussian0.exponent_coeff, gaussian1.exponent_coeff, i0, i1);

            const auto combinatoric_part = static_cast<double>(choose_term0 * choose_term1 * factorial_term);
            const auto gaussian_part = gauss0_contrib * gauss1_contrib / gauss_1d_coeff;

            unnormalized_overlap += (combinatoric_part * gaussian_part);
        }
    }

    return unnormalized_overlap;
}

auto overlap_integral_contraction(
    const AngularMomentumNumbers& angmom0,
    const AngularMomentumNumbers& angmom1,
    const coord::Cartesian3D& position0,
    const coord::Cartesian3D& position1,
    double exponent0,
    double exponent1
) -> double
{
    const auto [pos_product, coeff_product] = elec::math::gaussian_product(position0, position1, exponent0, exponent1);

    const auto norm0 = elec::math::gaussian_norm(angmom0, exponent0);
    const auto norm1 = elec::math::gaussian_norm(angmom1, exponent1);
    const auto overlap_norm = overlap_integral_3d_norm(exponent0, exponent1);
    const auto total_norm = norm0 * norm1 * overlap_norm;

    // clang-format off
    const auto unorm_overlap_x = unnormalized_overlap_integral_1d(
        {angmom0.x, exponent0, position0.x},
        {angmom1.x, exponent1, position1.x},
        pos_product.x
    );

    const auto unorm_overlap_y = unnormalized_overlap_integral_1d(
        {angmom0.y, exponent0, position0.y},
        {angmom1.y, exponent1, position1.y},
        pos_product.y
    );

    const auto unorm_overlap_z = unnormalized_overlap_integral_1d(
        {angmom0.z, exponent0, position0.z},
        {angmom1.z, exponent1, position1.z},
        pos_product.z
    );
    // clang-format on

    return coeff_product * unorm_overlap_x * unorm_overlap_y * unorm_overlap_z * total_norm;
}

}  // namespace elec
