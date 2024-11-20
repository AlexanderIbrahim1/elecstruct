#pragma once

#include <array>
#include <numeric>
#include <vector>

#include <extern/mapbox/eternal.hpp>

#include "elecstruct/atoms/atoms.hpp"
#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/orbitals.hpp"

namespace elec
{

namespace impl_elec
{

inline auto total_number_of_atomic_orbitals(const std::vector<elec::AtomInfo>& atom_infos) -> std::size_t
{
    auto n_orbitals = std::size_t {0};
    for (const auto& atom_info : atom_infos) {
        n_orbitals += atom_info.orbitals.size();
    }

    return n_orbitals;
}

}  // namespace impl_elec

struct GaussianConstantsSTO3G
{
    double coeff0;
    double coeff1;
    double coeff2;
    double expon0;
    double expon1;
    double expon2;
};

constexpr auto zeta_orbital1_exponents_sto3g = mapbox::eternal::map<AtomLabel, double>({
    {AtomLabel::H,  1.24  },
    {AtomLabel::He, 2.0925},
    {AtomLabel::Li, 2.69  },
    {AtomLabel::Be, 3.68  },
    {AtomLabel::B,  4.68  },
    {AtomLabel::C,  5.67  },
    {AtomLabel::N,  6.67  },
    {AtomLabel::O,  7.66  },
    {AtomLabel::F,  8.65  }
});

constexpr auto zeta_orbital2_exponents_sto3g = mapbox::eternal::map<AtomLabel, double>({
    {AtomLabel::Li, 0.75},
    {AtomLabel::Be, 1.10},
    {AtomLabel::B,  1.45},
    {AtomLabel::C,  1.72},
    {AtomLabel::N,  1.95},
    {AtomLabel::O,  2.25},
    {AtomLabel::F,  2.55}
});

// clang-format off
// coeff0, coeff1, coeff2, expon0, expon1, expon2
constexpr auto gaussian_constants_sto3g = mapbox::eternal::map<AtomicOrbitalLabel, GaussianConstantsSTO3G> ({
    {AtomicOrbitalLabel::S1, { 0.4446345422e+00,  0.5353281423e+00,  0.1543289673e+00,  0.109818e+00,  0.405771e+00,  0.222766e+01}},
    {AtomicOrbitalLabel::S2, { 0.7001154689e+00,  0.3995128261e+00, -0.9996722919e-01,  0.751386e-01,  0.231031e+00,  0.994203e+00}},
    {AtomicOrbitalLabel::P2, { 0.3919573931E+00,  0.6076837186E+00,  0.1559162750E+00,  0.751386e-01,  0.231031e+00,  0.994203e+00}}
});
// clang-format on

struct AtomicOrbitalInfoSTO3G
{
    AtomLabel atom_label;
    AtomicOrbitalLabel orbital_label;
    coord::Cartesian3D position;
    AngularMomentumNumbers angular_momentum;
    std::array<GaussianContractionInfo, 3> gaussians;
};

inline auto create_atomic_orbitals_sto3g(const std::vector<AtomInfo>& atom_infos) -> std::vector<AtomicOrbitalInfoSTO3G>
{
    using AOL = AtomicOrbitalLabel;

    const auto n_orbitals = impl_elec::total_number_of_atomic_orbitals(atom_infos);

    auto output = std::vector<AtomicOrbitalInfoSTO3G> {};
    output.reserve(n_orbitals);

    for (const auto& atom_info : atom_infos) {
        for (auto orbital : atom_info.orbitals) {
            const auto atom_label = atom_info.label;
            const auto angular_momenta = atomic_orbitals_to_angular_momentum_numbers(orbital);

            const auto zeta = [&]()
            {
                if (orbital == AOL::S1) {
                    return zeta_orbital1_exponents_sto3g.at(atom_label);
                }
                else {
                    return zeta_orbital2_exponents_sto3g.at(atom_label);
                }
            }();

            const auto constants = gaussian_constants_sto3g.at(orbital);
            const auto gauss0 = GaussianContractionInfo {constants.coeff0, constants.expon0 * zeta * zeta};
            const auto gauss1 = GaussianContractionInfo {constants.coeff1, constants.expon1 * zeta * zeta};
            const auto gauss2 = GaussianContractionInfo {constants.coeff2, constants.expon2 * zeta * zeta};

            for (const auto& momentum : angular_momenta) {
                output.push_back({
                    atom_label, orbital, atom_info.position, momentum, {gauss0, gauss1, gauss2}
                });
            }
        }
    }

    return output;
}

}  // namespace elec
