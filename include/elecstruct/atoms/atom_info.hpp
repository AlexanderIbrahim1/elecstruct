#pragma once

#include <vector>

#include "elecstruct/atoms/atom_label.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/orbitals.hpp"

namespace elec
{

struct AtomInfo
{
    AtomLabel label;
    coord::Cartesian3D position;
    std::vector<AtomicOrbitalLabel> orbitals;
};

}  // namespace elec
