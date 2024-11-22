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

    ierhf::maybe_print(is_verbose, "Calculating 'kinetic_mtx'");
    auto kinetic_mtx = kinetic_matrix(basis);
    ierhf::maybe_print(is_verbose, kinetic_mtx, "kinetic_mtx");

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

    // auto nuclear_mtx = Eigen::MatrixXd {7, 7};
    // nuclear_mtx << -61.733,  -7.447,   0.000,   0.000,   0.019,  -1.778,  -1.778,
    //                 -7.447, -10.151,   0.000,   0.000,   0.226,  -3.920,  -3.920,
    //                  0.000,   0.000,  -9.926,   0.000,   0.000,   0.000,   0.000,
    //                  0.000,   0.000,   0.000, -10.152,   0.000,  -0.228,   0.228,
    //                  0.019,   0.226,   0.000,   0.000, -10.088,   0.184,   0.184,
    //                 -1.778,  -3.920,   0.000,  -0.228,   0.184,  -5.867,  -1.652,
    //                 -1.778,  -3.920,   0.000,   0.228,   0.184,  -1.652,  -5.867;

    ierhf::maybe_print(is_verbose, nuclear_mtx, "nuclear_mtx");

    const auto core_hamiltonian_mtx = kinetic_mtx + nuclear_mtx;
    ierhf::maybe_print(is_verbose, core_hamiltonian_mtx, "core_hamiltonian_mtx");

    ierhf::maybe_print(is_verbose, "Calculating 'transformation_mtx'");
    const auto transformation_mtx = transformation_matrix(overlap_mtx);
    ierhf::maybe_print(is_verbose, transformation_mtx, "transformation_mtx");

    ierhf::maybe_print(is_verbose, "Calculating 'two_electron_integrals'");
    const auto two_electron_integrals = two_electron_integral_grid(basis);
    ierhf::maybe_print(is_verbose, two_electron_integrals, "two_electron_integrals");

    // --- ITERATION 0 ---
    std::cout << "Performing iteration 0\n";

    // the initial guess is given to the Fock matrix
    std::cout << "Calculating the initial Fock matrix\n";
    // auto fock_mtx = core_hamiltonian_mtx.eval();
    const auto huckel_constant = double {1.75};  // TODO: move somewhere else!
    auto fock_mtx = huckel_guess(overlap_mtx, core_hamiltonian_mtx, huckel_constant);

    // std::cout << "Calculating expected density matrix\n";
    // std::cout << new_density_matrix(fock_mtx, transformation_mtx, n_electrons) << '\n';

    // // specific to this water example, taken from the PDF (since I can't seem to get it right)
    // auto prev_density_mtx = Eigen::MatrixXd {7, 7};
    // prev_density_mtx <<  2.108, -0.456,  0.000,  0.000, -0.104, -0.022, -0.022,
    //                     -0.456,  2.010,  0.000,  0.000,  0.618, -0.059, -0.059,
    //                      0.000,  0.000,  2.000,  0.000,  0.000,  0.000,  0.000,
    //                      0.000,  0.000,  0.000,  0.737,  0.000,  0.539, -0.539,
    //                     -0.104,  0.618,  0.000,  0.000,  1.215, -0.482, -0.482,
    //                     -0.022, -0.059,  0.000,  0.539, -0.482,  0.606, -0.183,
    //                     -0.022, -0.059,  0.000, -0.539, -0.482, -0.183,  0.606;

    // std::cout << "PREV DENSITY MTX\n";
    // std::cout << prev_density_mtx << '\n';

    // std::cout << "Calculating the total energy\n";
    // auto tot_energy = total_energy(prev_density_mtx, fock_mtx, core_hamiltonian_mtx, atoms);
    // std::cout << "Total energy = " << tot_energy << '\n';

    auto tot_energy = double {0.0};
    auto prev_density_mtx = Eigen::MatrixXd::Zero(overlap_mtx.cols(), overlap_mtx.rows()).eval();

    // --- REMAINING ITERATIONS ---
    for (std::size_t i_iter {1}; i_iter < n_max_iter; ++i_iter) {
        std::cout << "Performing iteration " << i_iter << '\n';

        // WORKS ? 
        std::cout << "Calculating the Fock matrix\n";
        fock_mtx = fock_matrix(prev_density_mtx, basis, two_electron_integrals, core_hamiltonian_mtx);
        std::cout << "FOCK MTX\n";
        std::cout << fock_mtx << "\n\n";

        // PROBABLY DOES NOT WORK
        std::cout << "Calculating the new density matrix\n";
        const auto curr_density_mtx = new_density_matrix(fock_mtx, transformation_mtx, n_electrons);
        std::cout << "NEW DENSITY MTX\n";
        std::cout << curr_density_mtx << '\n';

        std::cout << "Calculating the total energy\n";
        tot_energy = total_energy(prev_density_mtx, fock_mtx, core_hamiltonian_mtx, atoms);
        std::cout << "Total energy = " << tot_energy << '\n';

        std::cout << "Calculating the density matrix difference\n";
        const auto difference = density_matrix_difference(prev_density_mtx, curr_density_mtx);
        std::cout << "Density matrix difference = " << std::fixed << std::setprecision(12) << difference;

        if (difference < density_mtx_convergence) {
            std::cout << "Converged!\n";
            break;
        }

        prev_density_mtx = curr_density_mtx.eval();
        std::cout << "DENSITY MTX AFTER CHANGE\n";
        std::cout << prev_density_mtx << "\n\n";
    }

    std::cout << "Failed to converge!\n";
}

}  // namespace elec
