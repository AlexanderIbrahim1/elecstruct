#pragma once

#include <vector>

#include <Eigen/Dense>

#include "elecstruct/basis/basis_sets/sto3g.hpp"
#include "elecstruct/integrals/overlap_integrals.hpp"


namespace elec
{

using Basis = std::vector<AtomicOrbitalInfoSTO3G>;

// TODO: change the type of the gaussian coefficients to complex, because technically we take
// the product between the complex conjugate of `info0.coefficient1` and `info1.coefficient`
auto overlap_matrix_s(const Basis& basis) -> Eigen::MatrixXd
{
    const auto size = basis.size();

    auto output = Eigen::MatrixXd {size, size};

    for (std::size_t i0 {0}; i0 < size; ++i0) {
        const auto& basis0 = basis[i0];
        const auto pos0 = basis0.position;
        const auto angmom_0 = basis0.angular_momentum;
        for (std::size_t i1 {0}; i1 < size; ++i1) {
            const auto& basis1 = basis[i1];
            const auto pos1 = basis1.position;
            const auto angmom_1 = basis1.angular_momentum;

            for (const auto info0 : basis0.gaussians) {
                for (const auto info1 : basis1.gaussians) {
                    const auto coeff = info0.coefficient * info1.coefficient;
                    const auto overlap = overlap_integral(angmom_0, angmom_1, pos0, pos1, info0, info1);

                    output(static_cast<int>(i0), static_cast<int>(i1)) = coeff * overlap;
                }
            }
        }
    }

    return output;
}

// inline auto overlap_integral(
//     const AngularMomentumNumbers& angmom0,
//     const AngularMomentumNumbers& angmom1,
//     const coord::Cartesian3D& position0,
//     const coord::Cartesian3D& position1,
//     const GaussianInfo& info0,
//     const GaussianInfo& info1
// ) -> double

}  // namespace elec
