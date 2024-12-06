#pragma once

#include <Eigen/Dense>

/*
    NOTE: right now we only have the zero matrix, but I may add other ones later.
*/

namespace elec
{

/*
    A poor and non-physical, but very easy guess to make.
*/
auto zero_matrix(std::size_t size) -> Eigen::MatrixXd;

/*
    Approximate the Fock operator using the Extended Hückel Hamiltonian matrix.

    SOURCE:

    title: "Extended Hückel and Slater's rule initial guess for real space grid-based density functional theory"
    authors: M. Lee and K. Leiter and C. Eisner and J. Crone and J. Knap
    journal: Computational and Theoretical Chemistry
    volume: 1062
    year: 2015

    Equation 12, middle of the right column on page 25
*/
auto extended_huckel_guess(const Eigen::MatrixXd& overlap_mtx, const Eigen::MatrixXd& core_hamiltonian_mtx)
    -> Eigen::MatrixXd;

/*
    Approximate the Fock operator using the Core Hamiltonian matrix.
*/
auto core_hamiltonian_guess(const Eigen::MatrixXd& core_hamiltonian_mtx) -> Eigen::MatrixXd;

}  // namespace elec
