#pragma once

#include <cstddef>
#include <vector>

#include "elecstruct/atoms.hpp"
#include "elecstruct/basis/basis.hpp"
#include "elecstruct/input_file_parser/input_file_options.hpp"


namespace elec
{

void perform_restricted_hartree_fock(
    const std::vector<AtomInfo>& atoms,
    const std::vector<AtomicOrbitalInfoSTO3G>& basis,
    InitialFockGuess initial_fock,
    std::size_t n_electrons,
    std::size_t n_max_iter,
    double tolerance_change_density_matrix,
    Verbose is_verbose
);

}  // namespace elec
