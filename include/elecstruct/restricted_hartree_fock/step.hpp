#pragma once

#include <algorithm>
#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include "elecstruct/atoms/atoms.hpp"
#include "elecstruct/basis/basis_sets/sto3g.hpp"
#include "elecstruct/grids/grid4d.hpp"
#include "elecstruct/matrices.hpp"


namespace impl_elec::step
{

/*
    Create a new matrix with the sorted columns
*/
inline auto matrix_with_sorted_columns(const Eigen::MatrixXd& matrix, const std::vector<Eigen::Index>& indices) -> Eigen::MatrixXd
{
    const auto size = matrix.cols();

    if (size != static_cast<Eigen::Index>(indices.size())) {
        throw std::runtime_error {"Mismatch between number of columns in matrix, and number of indices."};
    }

    auto output = Eigen::MatrixXd {size, size};

    for (Eigen::Index i {0}; i < size; ++i) {
        output.col(i) = matrix.col(indices[i]);
    }

    return output;
}

inline auto sorted_indices(const Eigen::VectorXd& elements) -> std::vector<Eigen::Index>
{
    auto indices = [&]() {
        const auto size = static_cast<std::size_t>(elements.size());

        auto indices_ = std::vector<Eigen::Index> {};
        indices_.reserve(size);
        for (std::size_t i {0}; i < size; ++i) {
            indices_[i] = static_cast<Eigen::Index>(i);
        }

        return indices_;
    }();

    const auto compare = [&](Eigen::Index i0, Eigen::Index i1) { return elements[i0] < elements[i1]; };

    std::sort(indices.begin(), indices.end(), compare);

    return indices;
}

}  // namespace impl_elec::step


namespace elec
{

inline auto fock_matrix(
    const Eigen::MatrixXd& old_density_mtx,
    const std::vector<AtomicOrbitalInfoSTO3G>& basis,
    const std::vector<AtomInfo>& atoms,
    const grid::Grid4D& two_electron_integrals,
    const Eigen::MatrixXd& core_hamiltonian_mtx
) -> Eigen::MatrixXd
{
    const auto electron_electron_mtx = electron_electron_matrix(basis, atoms, old_density_mtx, two_electron_integrals);
    const auto fock_mtx = core_hamiltonian_mtx + electron_electron_mtx;

    return fock_mtx;
}

/*
    Performs an iteration of the Restricted Hartree-Fock method.

    Returns:
      - the new density matrix
      - the Fock matrix in the original basis (not the orthonormalized basis where S -> I)
      - the eigenvalues of the Fock matrix
*/
inline auto rhf_step(
    const Eigen::MatrixXd& fock_mtx,
    const Eigen::MatrixXd& basis_transformation_mtx,
    std::size_t n_electrons
) -> std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd>
{
    namespace ies = impl_elec::step;

    const auto fock_mtx_trans = basis_transformation_mtx.transpose() * (fock_mtx * basis_transformation_mtx);

    const auto eigensolver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> {fock_mtx_trans};
    if (eigensolver.info() != Eigen::Success) {
        throw std::runtime_error {"Failed to perform eigenvalue decomposition of Fock matrix."};
    }

    // get the sorted eigenvalues and eigenvectors
    auto eigenvalues = eigensolver.eigenvalues();
    const auto sorted_indices = ies::sorted_indices(eigenvalues);

    std::sort(eigenvalues.begin(), eigenvalues.end());

    const auto eigenvalue_mtx = eigenvalues.asDiagonal();
    const auto coefficient_mtx_trans = ies::matrix_with_sorted_columns(eigensolver.eigenvectors(), sorted_indices);
    const auto coefficient_mtx = basis_transformation_mtx * coefficient_mtx_trans;

    const auto new_density_mtx = density_matrix_restricted_hartree_fock(coefficient_mtx, n_electrons);
}


inline auto density_matrix_difference(
    const Eigen::MatrixXd& old_density_mtx,
    const Eigen::MatrixXd& new_density_mtx
) -> double
{
    const auto size = old_density_mtx.cols();

    auto difference = double {0.0};
    for (Eigen::Index i0 {0}; i0 < size; ++i0) {
        for (Eigen::Index i1 {0}; i1 < size; ++i1) {
            const auto element_diff = new_density_mtx(i0, i1) - old_density_mtx(i0, i1);
            difference += element_diff * element_diff;
        }
    }

    return std::sqrt(difference / 4.0);
}


inline auto electronic_energy(
    const Eigen::MatrixXd& density_mtx,
    const Eigen::MatrixXd& fock_mtx,
    const Eigen::MatrixXd& core_hamiltonain_mtx
) -> double
{
    const auto size = density_mtx.cols();

    auto energy = double {0.0};

    for (Eigen::Index i0 {0}; i0 < size; ++i0) {
        for (Eigen::Index i1 {0}; i1 < size; ++i1) {
            energy += 0.5 * density_mtx(i0, i1) * (fock_mtx(i0, i1) + core_hamiltonain_mtx(i0, i1));
        }
    }

    return energy;
}


}  // namespace elec
