#include <vector>

#include "elecstruct/elecstruct.hpp"

/*
    Get positions of an example H2O molecule, matching the example given in Szabo and Ostlund.
*/
auto h2o_positions() -> std::tuple<coord::Cartesian3D, coord::Cartesian3D, coord::Cartesian3D>
{
    const auto bohr_per_angstrom = 0.529177210544;

    const auto y_h = 0.751155 / bohr_per_angstrom;
    const auto z_h = 0.465285 / bohr_per_angstrom;
    const auto z_o = 0.116321 / bohr_per_angstrom;

    const auto pos_h_0 = coord::Cartesian3D {0.0, y_h, -z_h};
    const auto pos_o = coord::Cartesian3D {0.0, 0.0, z_o};
    const auto pos_h_1 = coord::Cartesian3D {0.0, -y_h, -z_h};

    return {pos_h_0, pos_h_1, pos_o};
}

auto main() -> int
{
    const auto [pos_h_0, pos_h_1, pos_o] = h2o_positions();

    using AOL = elec::AtomicOrbitalLabel;

    const auto atoms = std::vector<elec::AtomInfo> {
        elec::AtomInfo {elec::AtomLabel::O, pos_o,   {AOL::S1, AOL::S2, AOL::P2}},
        elec::AtomInfo {elec::AtomLabel::H, pos_h_0, {AOL::S1}                  },
        elec::AtomInfo {elec::AtomLabel::H, pos_h_1, {AOL::S1}                  }
    };

    const auto basis = elec::create_atomic_orbitals_sto3g(atoms);
    const auto initial_fock = elec::InitialFockGuess::ZERO_MATRIX;

    const auto n_electrons = std::size_t {10};
    const auto n_max_iter = std::size_t {100};
    const auto tolerance = double {1.0e-8};
    const auto verbose = elec::Verbose::TRUE;

    elec::perform_restricted_hartree_fock(atoms, basis, initial_fock, n_electrons, n_max_iter, tolerance, verbose);

    return 0;
}
