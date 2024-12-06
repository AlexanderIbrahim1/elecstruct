#pragma once

#include <any>
#include <unordered_map>
#include <vector>

#include "elecstruct/atoms.hpp"
#include "elecstruct/input_file_parser/input_file_options.hpp"

namespace elec
{

enum class InputFileKey
{
    ATOM_INFORMATION,
    INITIAL_FOCK_GUESS,
    MAX_HARTREE_FOCK_ITERATIONS,
    TOL_CHANGE_DENSITY_MATRIX,
    TOL_CHANGE_HARTREE_FOCK_ENERGY,
    N_ELECTRONS,
    VERBOSE
};

class ParsedInformation
{
public:
    auto operator[](InputFileKey key) -> std::any&;
    auto atom_information() const -> std::vector<AtomInfo>;
    auto initial_fock_guess() const -> InitialFockGuess;
    auto max_hartree_fock_iterations() const -> std::size_t;
    auto tol_change_density_matrix() const -> double;
    auto tol_change_hartree_fock_energy() const -> double;
    auto n_electrons() const -> std::size_t;
    auto verbose() const -> Verbose;

private:
    std::unordered_map<InputFileKey, std::any> info_ {};
};

}  // namespace elec
