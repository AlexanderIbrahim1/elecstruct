#pragma once

#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include "elecstruct/atoms/atoms.hpp"
#include "elecstruct/basis/basis_sets/sto3g.hpp"
#include "elecstruct/integrals/overlap_integrals.hpp"
#include "elecstruct/integrals/kinetic_integrals.hpp"
#include "elecstruct/integrals/nuclear_electron_integrals.hpp"


// --- TODO ----------
// 1. change the type of the gaussian coefficients to complex, because technically we take
//   the product between the complex conjugate of `info0.coefficient1` and `info1.coefficient`
//
// 2. these matrices are Hermitian, so we can save calculating half of the entries!
//
// 3. many of these matrices have nearly identical bodies
//  - it looks like a bit of `if constexpr` can cut down the repetition by a lot
// 
// 4. maybe turn `nuclear_electron_integral()` into a class with a callable member function?
// - the interface is very similar to the overlap and kinetic integrals, except for the nuclear position and charge
// - but those two variables stay the same when calculating a nuclear-electron integral for a given atom
//
// - if we fix those variables, then we can get three different matrices created with the same interface
//   - and we can pass it as a template argument and reduce a ton of code duplication 

namespace elec
{

using Basis = std::vector<AtomicOrbitalInfoSTO3G>;

inline auto overlap_matrix(const Basis& basis) -> Eigen::MatrixXd
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

inline auto transformation_matrix(const Eigen::MatrixXd& s_overlap) -> Eigen::MatrixXd
{
    const auto eigensolver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> {s_overlap};
    if (eigensolver.info() != Eigen::Success) {
        throw std::runtime_error {"Failed to perform eigenvalue decomposition on overlap matrix"};
    }

    const auto eigenvectors = eigensolver.eigenvectors();

    const auto eigenvalues = eigensolver.eigenvalues();
    const auto inv_sqrt_eigenvalues = eigenvalues.array().pow(-0.5).matrix();
    const auto diagonal_inv_sqrt = inv_sqrt_eigenvalues.asDiagonal();

    return eigenvectors * diagonal_inv_sqrt;
}

inline auto kinetic_matrix(const Basis& basis) -> Eigen::MatrixXd
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
                    const auto overlap = kinetic_integral(angmom_0, angmom_1, pos0, pos1, info0, info1);

                    output(static_cast<int>(i0), static_cast<int>(i1)) = coeff * overlap;
                }
            }
        }
    }

    return output;
}

inline auto nuclear_electron_matrix(const Basis& basis, const AtomInfo& atom) -> Eigen::MatrixXd
{
    const auto size = basis.size();
    const auto charge = atom_charge_map.at(atom.label);

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
                    const auto overlap = nuclear_electron_integral(angmom_0, angmom_1, pos0, pos1, atom.position, info0, info1, charge);

                    output(static_cast<int>(i0), static_cast<int>(i1)) = coeff * overlap;
                }
            }
        }
    }

    return output;
}


}  // namespace elec
