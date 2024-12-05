#include <any>
#include <unordered_map>
#include <vector>

#include "elecstruct/atoms.hpp"

#include "elecstruct/input_file_parser/parsed_information.hpp"

namespace elec
{

auto ParsedInformation::operator[](InputFileKey key) -> std::any&
{
    return info_[key];
}

auto ParsedInformation::atom_information() const -> std::vector<AtomInfo>
{
    using T = std::vector<AtomInfo>;
    const auto key = InputFileKey::ATOM_INFORMATION;

    if (info_.find(key) == info_.end()) {
        throw std::runtime_error {"ERROR: 'atom_information' has not been parsed."};
    }

    return std::any_cast<T>(info_.at(key));
}

auto ParsedInformation::initial_fock_guess() const -> InitialFockGuess
{
    using T = InitialFockGuess;
    const auto key = InputFileKey::INITIAL_FOCK_GUESS;

    if (info_.find(key) == info_.end()) {
        throw std::runtime_error {"ERROR: 'initial_fock_guess' has not been parsed."};
    }

    return std::any_cast<T>(info_.at(key));
}

auto ParsedInformation::max_hartree_fock_iterations() const -> std::size_t
{
    using T = std::size_t;
    const auto key = InputFileKey::MAX_HARTREE_FOCK_ITERATIONS;

    if (info_.find(key) == info_.end()) {
        throw std::runtime_error {"ERROR: 'max_hartree_fock_iterations' has not been parsed."};
    }

    return std::any_cast<T>(info_.at(key));
}

auto ParsedInformation::tol_change_density_matrix() const -> double
{
    using T = double;
    const auto key = InputFileKey::TOL_CHANGE_DENSITY_MATRIX;

    if (info_.find(key) == info_.end()) {
        throw std::runtime_error {"ERROR: 'tol_change_density_matrix' has not been parsed."};
    }

    return std::any_cast<T>(info_.at(key));
}

auto ParsedInformation::tol_change_hartree_fock_energy() const -> double
{
    using T = double;
    const auto key = InputFileKey::TOL_CHANGE_HARTREE_FOCK_ENERGY;

    if (info_.find(key) == info_.end()) {
        throw std::runtime_error {"ERROR: 'tol_change_hartree_fock_energy' has not been parsed."};
    }

    return std::any_cast<T>(info_.at(key));
}

}  // anonymous namespace
