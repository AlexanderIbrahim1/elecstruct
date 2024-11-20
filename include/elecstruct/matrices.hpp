#pragma once

// TODO: remove
#include <iostream>

#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include "elecstruct/atoms/atoms.hpp"
#include "elecstruct/grids/grid4d.hpp"
#include "elecstruct/basis/basis_sets/sto3g.hpp"
#include "elecstruct/integrals/overlap_integrals.hpp"
#include "elecstruct/integrals/kinetic_integrals.hpp"
#include "elecstruct/integrals/nuclear_electron_integrals.hpp"
#include "elecstruct/integrals/electron_electron_integrals.hpp"


// --- TODO ----------
// 1. change the type of the gaussian coefficients to complex, because technically we take
//   the product between the complex conjugate of `info0.contraction_coeff1` and `info1.contraction_coeff`
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

namespace impl_elec::matrices
{

// TODO: change coefficients of the integral so that they might be complex, since technically this is required
inline auto two_electron_integral_grid(
    const elec::AtomicOrbitalInfoSTO3G& orbital0,
    const elec::AtomicOrbitalInfoSTO3G& orbital1,
    const elec::AtomicOrbitalInfoSTO3G& orbital2,
    const elec::AtomicOrbitalInfoSTO3G& orbital3
) -> double
{
    const auto pos_gauss0 = orbital0.position;
    const auto pos_gauss1 = orbital1.position;
    const auto pos_gauss2 = orbital2.position;
    const auto pos_gauss3 = orbital3.position;
    const auto angmom_0 = orbital0.angular_momentum;
    const auto angmom_1 = orbital1.angular_momentum;
    const auto angmom_2 = orbital2.angular_momentum;
    const auto angmom_3 = orbital3.angular_momentum;

    // clang-format off
    auto output = double {0.0};
    for (const auto& info0 : orbital0.gaussians)
    for (const auto& info1 : orbital1.gaussians)
    for (const auto& info2 : orbital2.gaussians)
    for (const auto& info3 : orbital3.gaussians) {
        const auto coeff = info0.contraction_coeff * info1.contraction_coeff * info2.contraction_coeff * info3.contraction_coeff;
        const auto integral = elec::electron_electron_integral(
            angmom_0, angmom_1, angmom_2, angmom_3,
            pos_gauss0, pos_gauss1, pos_gauss2, pos_gauss3,
            info0.exponent_coeff, info1.exponent_coeff, info2.exponent_coeff, info3.exponent_coeff
        );

        output += coeff * integral;
    }
    // clang-format on

    return output;
}

}  // namespace elec::impl_elec::matrices


inline auto overlap_matrix(const std::vector<AtomicOrbitalInfoSTO3G>& basis) -> Eigen::MatrixXd
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

            auto element = double {0.0};
            for (const auto info0 : basis0.gaussians) {
                for (const auto info1 : basis1.gaussians) {
                    const auto coeff = info0.contraction_coeff * info1.contraction_coeff;
                    const auto overlap = overlap_integral(angmom_0, angmom_1, pos0, pos1, info0.exponent_coeff, info1.exponent_coeff);

                    element += coeff * overlap;
                }
            }
            output(static_cast<Eigen::Index>(i0), static_cast<Eigen::Index>(i1)) = element;
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

inline auto kinetic_matrix(const std::vector<AtomicOrbitalInfoSTO3G>& basis) -> Eigen::MatrixXd
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

            auto element = double {0.0};
            for (const auto info0 : basis0.gaussians) {
                for (const auto info1 : basis1.gaussians) {
                    const auto coeff = info0.contraction_coeff * info1.contraction_coeff;
                    const auto overlap = kinetic_integral(angmom_0, angmom_1, pos0, pos1, info0.exponent_coeff, info1.exponent_coeff);
                    element += coeff * overlap;
                }
            }
            output(static_cast<Eigen::Index>(i0), static_cast<Eigen::Index>(i1)) = element;
        }
    }

    return output;
}

inline auto nuclear_electron_matrix(const std::vector<AtomicOrbitalInfoSTO3G>& basis, const AtomInfo& atom) -> Eigen::MatrixXd
{
    const auto size = basis.size();
    const auto charge = atom_charge_map.at(atom.label);

    auto output = Eigen::MatrixXd {size, size};

    // clang-format off
    for (std::size_t i0 {0}; i0 < size; ++i0) {
        const auto& basis0 = basis[i0];
        const auto pos0 = basis0.position;
        const auto angmom_0 = basis0.angular_momentum;
        for (std::size_t i1 {0}; i1 < size; ++i1) {
            const auto& basis1 = basis[i1];
            const auto pos1 = basis1.position;
            const auto angmom_1 = basis1.angular_momentum;

            auto element = double {0.0};
            for (const auto info0 : basis0.gaussians) {
                for (const auto info1 : basis1.gaussians) {
                    const auto coeff = info0.contraction_coeff * info1.contraction_coeff;
                    const auto overlap = nuclear_electron_integral(
                        angmom_0, angmom_1,
                        pos0, pos1, atom.position,
                        info0.exponent_coeff, info1.exponent_coeff,
                        charge
                    );

                    element += coeff * overlap;
                }
            }
            std::cout << "NE(i0, i1) = (" << i0 << ", " << i1 << ") = " << element << '\n';

            output(static_cast<Eigen::Index>(i0), static_cast<Eigen::Index>(i1)) = element;
        }
    }
    // clang-format on

    return output;
}


/*
    Calculates the core Hamiltonian matrix, which is the sum of the kinetic energy matrix
    and all the nuclear-electron interaction matrices.
*/
inline auto core_hamiltonian_matrix(
    const std::vector<AtomicOrbitalInfoSTO3G>& basis,
    const std::vector<AtomInfo>& atoms
) -> Eigen::MatrixXd
{
    auto output = kinetic_matrix(basis);

    std::cout << "kinetic_mtx\n";
    std::cout << output << "\n\n";

    for (const auto& atom : atoms) {
        const auto nuclear_mtx = nuclear_electron_matrix(basis, atom);

        std::cout << "nuclear_mtx\n";
        std::cout << nuclear_mtx << "\n\n";

        output += nuclear_mtx;
    }

    return output;
}


/*
    NOTE: I'm pretty sure exchanging i0 <-> i1, and i2 <-> i3, only change the result by a complex
    conjugate; once I get the code running, I should take advantage of that symmetry
*/
inline auto two_electron_integral_grid(const std::vector<AtomicOrbitalInfoSTO3G>& basis) -> grid::Grid4D
{
    namespace emat = impl_elec::matrices;

    const auto size = basis.size();
    auto integral_grid = grid::Grid4D {size, size, size, size};

    for (std::size_t i0 {0}; i0 < size; ++i0)
    for (std::size_t i1 {0}; i1 < size; ++i1)
    for (std::size_t i2 {0}; i2 < size; ++i2)
    for (std::size_t i3 {0}; i3 < size; ++i3) {
        integral_grid.set(i0, i1, i2, i3, emat::two_electron_integral_grid(basis[i0], basis[i1], basis[i2], basis[i3]));
    }

    return integral_grid;
}


/*
    NOTE: the density elements are actually given by the sum of:
        `coefficient_mtx(i0, j)` * CONJUGATE(coefficient_mtx(i1, j))`;

    I think the examples used here all have real elements, so it won't matter
      - but it might change if the orbitals start off using complex coefficients
*/
inline auto density_matrix_restricted_hartree_fock(const Eigen::MatrixXd& coefficient_mtx, std::size_t n_electrons) -> Eigen::MatrixXd
{
    const auto size = coefficient_mtx.cols();
    const auto half = static_cast<Eigen::Index>(n_electrons / 2);

    auto output = Eigen::MatrixXd {size, size};

    for (Eigen::Index i0 {0}; i0 < size; ++i0) {
        for (Eigen::Index i1 {0}; i1 < size; ++i1) {
            auto density_element = double {0.0};
            for (Eigen::Index j {0}; j < half; ++j) {
                density_element += coefficient_mtx(i0, j) * coefficient_mtx(i1, j);
            }
            output(i0, i1) = 2.0 * density_element;
        }
    }

    return output;
}

inline auto electron_electron_matrix(
    const std::vector<AtomicOrbitalInfoSTO3G>& basis,
    const Eigen::MatrixXd& density_matrix,
    const grid::Grid4D& two_electron_integrals
) -> Eigen::MatrixXd
{
    const auto size = basis.size();
    auto output = Eigen::MatrixXd {size, size};

    for (std::size_t i0 {0}; i0 < size; ++i0) {
        for (std::size_t i1 {0}; i1 < size; ++i1) {
            auto element = double {0.0};
            for (std::size_t i2 {0}; i2 < size; ++i2) {
                for (std::size_t i3 {0}; i3 < size; ++i3) {
                    const auto i2_eig = static_cast<Eigen::Index>(i2);
                    const auto i3_eig = static_cast<Eigen::Index>(i3);
                    const auto density_part = density_matrix(i2_eig, i3_eig);

                    const auto coulombic_term = two_electron_integrals.get(i0, i1, i2, i3);
                    const auto exchange_term = two_electron_integrals.get(i0, i3, i2, i1);
                    const auto electron_part = coulombic_term - 0.5 * exchange_term;

                    element += density_part * electron_part;
                }
            }
            output(static_cast<Eigen::Index>(i0), static_cast<Eigen::Index>(i1)) = element;
        }
    }

    return output;
}

}  // namespace elec
