#pragma once

#include <iomanip>
#include <iostream>
#include <string_view>
#include <vector>

#include "elecstruct/atoms.hpp"
#include "elecstruct/basis/basis.hpp"
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

inline void maybe_print(Verbose is_verbose, const TwoElectronIntegralGrid& grid, std::size_t n_basis_functions, std::string_view name)
{
    if (is_verbose == Verbose::TRUE) {
        std::cout << name << '\n';

        for (std::size_t i0 {0}; i0 < n_basis_functions; ++i0)
        for (std::size_t i1 {0}; i1 < n_basis_functions; ++i1)
        for (std::size_t i2 {0}; i2 < n_basis_functions; ++i2)
        for (std::size_t i3 {0}; i3 < n_basis_functions; ++i3) {
            std::cout << "(" << i0 << ", " << i1 << ", " << i2 << ", " << i3 << ") = ";
            std::cout << grid.get(i0, i1, i2, i3) << '\n';
        }
    }
}

inline void maybe_print_divider(Verbose is_verbose)
{
    if (is_verbose == Verbose::TRUE) {
        std::cout << "----------------------------------------------------------------------------------------------------\n";
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

    // ------------------------------------------------------------------------
    ierhf::maybe_print(is_verbose, "Calculating 'overlap_mtx'");
    const auto overlap_mtx = overlap_matrix(basis);
    ierhf::maybe_print(is_verbose, overlap_mtx, "overlap_mtx");
    ierhf::maybe_print_divider(is_verbose);

    // ------------------------------------------------------------------------
    ierhf::maybe_print(is_verbose, "Calculating 'kinetic_mtx'");
    auto kinetic_mtx = kinetic_matrix(basis);
    ierhf::maybe_print(is_verbose, kinetic_mtx, "kinetic_mtx");
    ierhf::maybe_print_divider(is_verbose);

    // ------------------------------------------------------------------------
    ierhf::maybe_print(is_verbose, "Calculating 'nuclear_electron_mtx'");
    auto nuclear_mtx = Eigen::MatrixXd::Zero(kinetic_mtx.cols(), kinetic_mtx.cols()).eval();

    for (const auto& atom : atoms) {
        if (is_verbose == Verbose::TRUE) {
            const auto atom_name = elec::atom_name_from_label(atom.label);
            std::cout << "Calculating 'nuclear_electron_mtx' term for atom '" << atom_name << "'\n";
        }
        const auto atom_nuclear_electron_mtx = nuclear_electron_matrix(basis, atom);
        ierhf::maybe_print(is_verbose, atom_nuclear_electron_mtx, "atom_nuclear_electron_mtx");

        nuclear_mtx += atom_nuclear_electron_mtx;
    }

    ierhf::maybe_print(is_verbose, nuclear_mtx, "nuclear_mtx");
    ierhf::maybe_print_divider(is_verbose);

    // ------------------------------------------------------------------------
    ierhf::maybe_print(is_verbose, "Calculating 'core_hamiltonian_mtx'");
    const auto core_hamiltonian_mtx = kinetic_mtx + nuclear_mtx;
    ierhf::maybe_print(is_verbose, core_hamiltonian_mtx, "core_hamiltonian_mtx");
    ierhf::maybe_print_divider(is_verbose);

    // ------------------------------------------------------------------------
    ierhf::maybe_print(is_verbose, "Calculating 'transformation_mtx'");
    const auto transformation_mtx = transformation_matrix(overlap_mtx);
    ierhf::maybe_print(is_verbose, transformation_mtx, "transformation_mtx");
    ierhf::maybe_print_divider(is_verbose);

    // ------------------------------------------------------------------------
    ierhf::maybe_print(is_verbose, "Calculating 'two_electron_integrals'");
    const auto two_electron_integrals = two_electron_integral_grid(basis);
    // ierhf::maybe_print(is_verbose, two_electron_integrals, basis.size(), "two_electron_integrals");
    ierhf::maybe_print_divider(is_verbose);

    ierhf::maybe_print_divider(is_verbose);
    ierhf::maybe_print_divider(is_verbose);

    // --- ITERATION 0 ---
    std::cout << "\nPerforming iteration 0\n";

    std::cout << "Calculating the initial Fock matrix\n";
    // auto fock_mtx = core_hamiltonian_mtx.eval();
    const auto huckel_constant = double {1.75};
    auto fock_mtx = huckel_guess(overlap_mtx, core_hamiltonian_mtx, huckel_constant);
    // auto fock_mtx = zero_matrix(basis.size());
    ierhf::maybe_print_divider(is_verbose);

    std::cout << "Calculating initial density matrix\n";
    auto prev_density_mtx = new_density_matrix(fock_mtx, transformation_mtx, n_electrons);

    auto tot_energy = total_energy(prev_density_mtx, fock_mtx, core_hamiltonian_mtx, atoms);
    std::cout << "Total energy = " << tot_energy << '\n';

    // --- REMAINING ITERATIONS ---
    for (std::size_t i_iter {1}; i_iter < n_max_iter; ++i_iter) {
        std::cout << "\nPerforming iteration " << i_iter << '\n';

        std::cout << "Calculating the Fock matrix\n";
        fock_mtx = fock_matrix(prev_density_mtx, basis, two_electron_integrals, core_hamiltonian_mtx);

        std::cout << "Calculating the density matrix\n";
        const auto density_mtx = new_density_matrix(fock_mtx, transformation_mtx, n_electrons);

        tot_energy = total_energy(density_mtx, fock_mtx, core_hamiltonian_mtx, atoms);
        std::cout << "Total energy = " << tot_energy << '\n';

        std::cout << "Calculating the density matrix difference\n";
        const auto difference = density_matrix_difference(prev_density_mtx, density_mtx);
        std::cout << "Density matrix difference = " << std::fixed << std::setprecision(12) << difference << '\n';

        if (difference < density_mtx_convergence) {
            std::cout << "\nConverged!\n";
            std::cout << "The total energy is " << tot_energy << " A.U. \n";
            return;
        }

        prev_density_mtx = density_mtx.eval();

        ierhf::maybe_print_divider(is_verbose);
    }

    std::cout << "Failed to converge!\n";
}

}  // namespace elec
