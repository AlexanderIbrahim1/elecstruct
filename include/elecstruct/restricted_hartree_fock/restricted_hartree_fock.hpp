#pragma once

#include <cstddef>
#include <vector>

#include "elecstruct/atoms.hpp"
#include "elecstruct/basis/basis.hpp"


namespace elec
{

enum class Verbose
{
    TRUE,
    FALSE
};

void perform_restricted_hartree_fock(
    const std::vector<AtomInfo>& atoms,
    const std::vector<AtomicOrbitalInfoSTO3G>& basis,
    std::size_t n_electrons,
    std::size_t n_max_iter,
    double density_mtx_convergence,
    Verbose is_verbose = Verbose::TRUE
);

}  // namespace elec
