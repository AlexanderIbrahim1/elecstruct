#pragma once

#include <any>
#include <unordered_map>
#include <vector>

#include "elecstruct/atoms.hpp"


namespace elec
{

enum class InputFileKey
{
    ATOM_INFORMATION,
    INITIAL_FOCK_GUESS,
    N_MAX_HARTREE_FOCK_ITERATIONS,
    TOL_CHANGE_DENSITY_MATRIX,
    TOL_CHANGE_HARTREE_FOCK_ENERGY
};

class ParsedInformation
{
public:
    auto operator[](InputFileKey key) -> std::any&;
    auto atom_information() const -> std::vector<AtomInfo>;

private:
    std::unordered_map<InputFileKey, std::any> info_ {};
};

}  // anonymous namespace
