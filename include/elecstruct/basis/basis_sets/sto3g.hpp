#pragma once

#include <array>
#include <vector>

#include "elecstruct/atoms.hpp"
#include "elecstruct/basis/gaussian_contraction_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/orbitals.hpp"

namespace elec
{

struct AtomicOrbitalInfoSTO3G
{
    AtomLabel atom_label;
    AtomicOrbitalLabel orbital_label;
    coord::Cartesian3D position;
    AngularMomentumNumbers angular_momentum;
    std::array<GaussianContractionInfo, 3> gaussians;
};

void fill_atomic_orbitals_sto3g(std::vector<AtomInfo>& atom_infos);
auto create_atomic_orbitals_sto3g(const std::vector<AtomInfo>& atom_infos) -> std::vector<AtomicOrbitalInfoSTO3G>;

}  // namespace elec
