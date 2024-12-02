#pragma once

#include <string>
#include <vector>

#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/orbitals.hpp"

namespace elec
{

enum class AtomLabel
{
    H,
    He,
    Li,
    Be,
    B,
    C,
    N,
    O,
    F
};

struct AtomInfo
{
    AtomLabel label;
    coord::Cartesian3D position;
    std::vector<AtomicOrbitalLabel> orbitals;
};

auto nuclear_charge(AtomLabel label) -> double;

auto atom_label_from_name(const std::string&) -> AtomLabel;

auto atom_name_from_label(AtomLabel) -> std::string;

}  // namespace elec
