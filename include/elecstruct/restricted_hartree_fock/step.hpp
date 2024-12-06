#pragma once

#include <vector>

#include <Eigen/Dense>

#include "elecstruct/atoms.hpp"

namespace elec
{

/*
    Create a new matrix with the sorted columns.
*/
auto matrix_with_sorted_columns(const Eigen::MatrixXd& matrix, const std::vector<Eigen::Index>& indices)
    -> Eigen::MatrixXd;

/*
    Returns a vector of indices such that if the elements in the Eigen::VectorXd `elements` were rearranged
    to these indices, then it would be sorted.
*/
auto indices_to_sort(const Eigen::VectorXd& elements) -> std::vector<Eigen::Index>;

/*
    Calculate the density matrix consistent with the provided Fock matrix.
*/
auto new_density_matrix(
    const Eigen::MatrixXd& fock_mtx,
    const Eigen::MatrixXd& basis_transformation_mtx,
    std::size_t n_electrons
) -> Eigen::MatrixXd;

/*
    Calculate the differences between two density matrices; used in iteration procedures to
    determine how similar the current density matrix is to the density matrix from the previous
    step.
*/
auto density_matrix_difference(const Eigen::MatrixXd& old_density_mtx, const Eigen::MatrixXd& new_density_mtx)
    -> double;

/*
    This function returns the sum of the:
      - electron-electron repulsion energies
      - nuclear-electron attraction energies
*/
auto electron_energy(
    const Eigen::MatrixXd& density_mtx,
    const Eigen::MatrixXd& fock_mtx,
    const Eigen::MatrixXd& core_hamiltonain_mtx
) -> double;

/*
    This function returns the sum of the nuclear-nuclear repulsion energies.
*/
auto nuclear_energy(const std::vector<AtomInfo>& atoms) -> double;

/*
    This function returns the sum of the:
      - nuclear-nuclear repulsion energies
      - electron-electron repulsion energies
      - nuclear-electron attraction energies
*/
auto total_energy(
    const Eigen::MatrixXd& density_mtx,
    const Eigen::MatrixXd& fock_mtx,
    const Eigen::MatrixXd& core_hamiltonain_mtx,
    const std::vector<AtomInfo>& atoms
) -> double;

}  // namespace elec
