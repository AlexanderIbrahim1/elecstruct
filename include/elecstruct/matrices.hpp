#pragma once

#include <vector>

#include <Eigen/Dense>

#include "elecstruct/atoms.hpp"
#include "elecstruct/grids/two_electron_integral_grid.hpp"
#include "elecstruct/basis/basis_sets/sto3g.hpp"


namespace elec
{

auto overlap_matrix(const std::vector<AtomicOrbitalInfoSTO3G>& basis) -> Eigen::MatrixXd;

auto transformation_matrix(const Eigen::MatrixXd& s_overlap) -> Eigen::MatrixXd;

auto kinetic_matrix(const std::vector<AtomicOrbitalInfoSTO3G>& basis) -> Eigen::MatrixXd;

auto nuclear_electron_matrix(const std::vector<AtomicOrbitalInfoSTO3G>& basis, const AtomInfo& atom) -> Eigen::MatrixXd;

/*
    Calculates the core Hamiltonian matrix, which is the sum of the kinetic energy matrix
    and all the nuclear-electron interaction matrices.
*/
auto core_hamiltonian_matrix(
    const std::vector<AtomicOrbitalInfoSTO3G>& basis,
    const std::vector<AtomInfo>& atoms
) -> Eigen::MatrixXd;

auto two_electron_integral_grid(const std::vector<AtomicOrbitalInfoSTO3G>& basis) -> TwoElectronIntegralGrid;

auto density_matrix_restricted_hartree_fock(const Eigen::MatrixXd& coefficient_mtx, std::size_t n_electrons) -> Eigen::MatrixXd;

auto electron_electron_matrix(
    const std::vector<AtomicOrbitalInfoSTO3G>& basis,
    const Eigen::MatrixXd& density_matrix,
    const TwoElectronIntegralGrid& two_electron_integrals
) -> Eigen::MatrixXd;

auto fock_matrix(
    const Eigen::MatrixXd& old_density_mtx,
    const std::vector<AtomicOrbitalInfoSTO3G>& basis,
    const TwoElectronIntegralGrid& two_electron_integrals,
    const Eigen::MatrixXd& core_hamiltonian_mtx
) -> Eigen::MatrixXd;

}  // namespace elec
