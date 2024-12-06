#include <array>
#include <numeric>
#include <vector>

#include <extern/mapbox/eternal.hpp>

#include "elecstruct/atoms.hpp"
#include "elecstruct/basis/gaussian_contraction_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/orbitals.hpp"

#include "elecstruct/basis/basis_sets/sto3g.hpp"

namespace
{

auto total_number_of_atomic_orbitals(const std::vector<elec::AtomInfo>& atom_infos) -> std::size_t
{
    auto n_orbitals = std::size_t {0};
    for (const auto& atom_info : atom_infos) {
        n_orbitals += atom_info.orbitals.size();
    }

    return n_orbitals;
}

struct GaussianConstantsSTO3G
{
    double coeff0;
    double coeff1;
    double coeff2;
    double expon0;
    double expon1;
    double expon2;
};

constexpr auto zeta_orbital1_exponents_sto3g = mapbox::eternal::map<elec::AtomLabel, double>({
    {elec::AtomLabel::H,  1.24  },
    {elec::AtomLabel::He, 2.0925},
    {elec::AtomLabel::Li, 2.69  },
    {elec::AtomLabel::Be, 3.68  },
    {elec::AtomLabel::B,  4.68  },
    {elec::AtomLabel::C,  5.67  },
    {elec::AtomLabel::N,  6.67  },
    {elec::AtomLabel::O,  7.66  },
    {elec::AtomLabel::F,  8.65  }
});

constexpr auto zeta_orbital2_exponents_sto3g = mapbox::eternal::map<elec::AtomLabel, double>({
    {elec::AtomLabel::Li, 0.75},
    {elec::AtomLabel::Be, 1.10},
    {elec::AtomLabel::B,  1.45},
    {elec::AtomLabel::C,  1.72},
    {elec::AtomLabel::N,  1.95},
    {elec::AtomLabel::O,  2.25},
    {elec::AtomLabel::F,  2.55}
});

// clang-format off
// coeff0, coeff1, coeff2, expon0, expon1, expon2
constexpr auto gaussian_constants_sto3g = mapbox::eternal::map<elec::AtomicOrbitalLabel, GaussianConstantsSTO3G> ({
    {elec::AtomicOrbitalLabel::S1, { 0.4446345422e+00,  0.5353281423e+00,  0.1543289673e+00,  0.109818e+00,  0.405771e+00,  0.222766e+01}},
    {elec::AtomicOrbitalLabel::S2, { 0.7001154689e+00,  0.3995128261e+00, -0.9996722919e-01,  0.751386e-01,  0.231031e+00,  0.994203e+00}},
    {elec::AtomicOrbitalLabel::P2, { 0.3919573931E+00,  0.6076837186E+00,  0.1559162750E+00,  0.751386e-01,  0.231031e+00,  0.994203e+00}}
});
// clang-format on

}  // anonymous namespace


namespace elec
{

void fill_atomic_orbitals_sto3g(std::vector<AtomInfo>& atom_infos)
{
    using AL = AtomLabel;
    using AOL = AtomicOrbitalLabel;

    for (auto& atom : atom_infos) {
        switch (atom.label)
        {
            case AL::H:
            case AL::He: {
                atom.orbitals.push_back(AOL::S1);
                break;
            };
            case AL::Li:
            case AL::Be:
            case AL::B:
            case AL::C:
            case AL::N:
            case AL::O:
            case AL::F: {
                atom.orbitals.push_back(AOL::S1);
                atom.orbitals.push_back(AOL::S2);
                atom.orbitals.push_back(AOL::P2);
                break;
            }
            default: {
                throw std::runtime_error {"UNREACHABLE: no atoms beyond 'F' on the periodic table have been implemented yet."};
            }
        }
    }
}

auto create_atomic_orbitals_sto3g(const std::vector<AtomInfo>& atom_infos) -> std::vector<AtomicOrbitalInfoSTO3G>
{
    using AOL = AtomicOrbitalLabel;

    const auto n_orbitals = total_number_of_atomic_orbitals(atom_infos);

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
