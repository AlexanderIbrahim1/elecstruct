#pragma once

#include <iomanip>
#include <iostream>
#include <string_view>
#include <vector>

#include "elecstruct/atoms/atoms.hpp"
#include "elecstruct/basis/basis_sets/sto3g.hpp"
#include "elecstruct/grids/grid4d.hpp"
#include "elecstruct/matrices.hpp"
#include "elecstruct/restricted_hartree_fock/initial_density_matrix.hpp"
#include "elecstruct/restricted_hartree_fock/step.hpp"


namespace elec
{

enum class Verbose
{
    TRUE,
    FALSE
};

namespace impl_elec::restricted_hartree_fock
{

inline void maybe_print(Verbose is_verbose, const Eigen::MatrixXd& matrix, std::string_view name)
{
    if (is_verbose == Verbose::TRUE) {
        std::cout << name << '\n';
        std::cout << matrix << "\n\n";
    }
}

inline void maybe_print(Verbose is_verbose, std::string_view text)
{
    if (is_verbose == Verbose::TRUE) {
        std::cout << text << "\n\n";
    }
}

inline void maybe_print(Verbose is_verbose, const grid::Grid4D& grid, std::string_view name)
{
    if (is_verbose == Verbose::TRUE) {
        const auto [size0, size1, size2, size3] = grid.sizes();

        std::cout << name << '\n';

        for (std::size_t i0 {0}; i0 < size0; ++i0)
        for (std::size_t i1 {0}; i1 < size1; ++i1)
        for (std::size_t i2 {0}; i2 < size2; ++i2)
        for (std::size_t i3 {0}; i3 < size3; ++i3) {
            std::cout << "(" << i0 << ", " << i1 << ", " << i2 << ", " << i3 << ") = ";
            std::cout << grid.get(i0, i1, i2, i3) << '\n';
        }
    }
}

}  // namespace elec::impl_elec::restricted_hartree_fock

/*
    NOTE: what information should I return from the method besides the energy?
*/
inline void perform_restricted_hartree_fock(
    const std::vector<AtomInfo>& atoms,
    const std::vector<AtomicOrbitalInfoSTO3G>& basis,
    std::size_t n_electrons,
    std::size_t n_max_iter,
    double density_mtx_convergence,
    Verbose is_verbose = Verbose::TRUE
)
{
    namespace ierhf = impl_elec::restricted_hartree_fock;

    ierhf::maybe_print(is_verbose, "Calculating 'overlap_mtx'");
    const auto overlap_mtx = overlap_matrix(basis);
    ierhf::maybe_print(is_verbose, overlap_mtx, "overlap_mtx");

    ierhf::maybe_print(is_verbose, "Calculating 'transformation_mtx'");
    const auto transformation_mtx = transformation_matrix(overlap_mtx);
    ierhf::maybe_print(is_verbose, transformation_mtx, "transformation_mtx");

    ierhf::maybe_print(is_verbose, "Calculating 'core_hamiltonian_mtx'");
    const auto core_hamiltonian_mtx = core_hamiltonian_matrix(basis, atoms);
    ierhf::maybe_print(is_verbose, core_hamiltonian_mtx, "core_hamiltonian_mtx");

    std::exit(EXIT_FAILURE);

    ierhf::maybe_print(is_verbose, "Calculating 'two_electron_integrals'");
    const auto two_electron_integrals = two_electron_integral_grid(basis);
    ierhf::maybe_print(is_verbose, two_electron_integrals, "two_electron_integrals");

    auto prev_density_mtx = zero_matrix(basis.size());

    for (std::size_t i_iter {0}; i_iter < n_max_iter; ++i_iter) {
        std::cout << "Performing iteration " << i_iter << '\n';
        const auto fock_mtx = fock_matrix(prev_density_mtx, basis, two_electron_integrals, core_hamiltonian_mtx);
        const auto curr_density_mtx = new_density_matrix(fock_mtx, transformation_mtx, n_electrons);

        const auto tot_energy = total_energy(curr_density_mtx, fock_mtx, core_hamiltonian_mtx, atoms);
        std::cout << "Total energy = " << tot_energy << '\n';

        const auto difference = density_matrix_difference(prev_density_mtx, curr_density_mtx);
        std::cout << "Density matrix difference = " << std::fixed << std::setprecision(12) << difference;

        if (difference < density_mtx_convergence) {
            std::cout << "Converged!\n";
            break;
        }
    }

    std::cout << "Failed to converge!\n";
}

}  // namespace elec
