#pragma once

#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include "elecstruct/atoms/atoms.hpp"
#include "elecstruct/basis/basis_sets/sto3g.hpp"
#include "elecstruct/grids/grid4d.hpp"

namespace elec
{

/*
    Performs an iteration of the Restricted Hartree-Fock method.

    Returns:
      - the new density matrix
      - the Fock matrix in the original basis (not the orthonormalized basis where S -> I)
      - the eigenvalues of the Fock matrix
*/
inline auto rhf_step(
    const Eigen::MatrixXd& old_density_mtx,
    const std::vector<AtomicOrbitalInfoSTO3G>& basis,
    const std::vector<AtomInfo>& atoms,
    const Eigen::MatrixXd& core_hamiltonian_mtx,
    const Eigen::MatrixXd& basis_transformation_mtx,
    const grid::Grid4D& electron_electron_integrals,
    std::size_t n_electrons
) -> std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd>
{

}

}  // namespace elec
