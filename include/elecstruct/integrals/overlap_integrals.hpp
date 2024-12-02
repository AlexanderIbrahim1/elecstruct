#pragma once

#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/orbitals.hpp"

namespace elec
{

struct OverlapIntegralGaussianContractionInfo1D
{
    std::int64_t angular_momentum;
    double exponent_coeff;
    double centre;
};

auto overlap_integral_3d_norm(double exponent0, double exponent1) -> double;

auto unnormalized_overlap_integral_1d(
    const OverlapIntegralGaussianContractionInfo1D& gaussian0,
    const OverlapIntegralGaussianContractionInfo1D& gaussian1,
    double total_centre
) -> double;

auto overlap_integral_contraction(
    const AngularMomentumNumbers& angmom0,
    const AngularMomentumNumbers& angmom1,
    const coord::Cartesian3D& position0,
    const coord::Cartesian3D& position1,
    double exponent0,
    double exponent1
) -> double;

}  // namespace elec
