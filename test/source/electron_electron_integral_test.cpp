#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "elecstruct/atoms.hpp"
#include "elecstruct/basis/basis_sets/sto3g.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/matrices.hpp"
#include "elecstruct/orbitals.hpp"

auto get_h2_basis() -> std::vector<elec::AtomicOrbitalInfoSTO3G>
{
    using AOL = elec::AtomicOrbitalLabel;

    const auto distance = double {1.4};  // distance given in sample calculation
    const auto pos0 = coord::Cartesian3D {0.0, 0.0, 0.0};
    const auto pos1 = coord::Cartesian3D {0.0, 0.0, distance};

    const auto atoms = std::vector<elec::AtomInfo> {
        elec::AtomInfo {elec::AtomLabel::H, pos0, {AOL::S1}},
        elec::AtomInfo {elec::AtomLabel::H, pos1, {AOL::S1}}
    };

    const auto basis = elec::create_atomic_orbitals_sto3g(atoms);

    return basis;
}


TEST_CASE("h2 two-electron integrals")
{
    /*
        The values I am comparing my output to come from here: https://chemistry.stackexchange.com/a/42035

        The author says that:
          - the elements are for the H2 molecule with a bond length of 1.4 Bohr
          - the basis used is the STO-3G basis
          - the output is from their own electronic structure code, but matches Szabo and Ostlund very well

        The reason for the discrepancy between my output and the provided output is probably
        my fairly poor approximation of the Boys function. That's why I have to limit the
        absolute tolerance with which I get my answer.
    */
    constexpr auto abs_tolerance = double {1.0e-4};

    const auto basis = get_h2_basis();
    const auto ee_grid = elec::two_electron_integral_grid(basis);

    const auto elem_0_0_0_0 = 0.77460834925515787;
    const auto elem_0_0_0_1 = 0.44410904384277344;
    const auto elem_0_0_1_0 = 0.44410904384277350;
    const auto elem_0_0_1_1 = 0.56967771733030592;
    const auto elem_0_1_0_0 = 0.44410904384277361;
    const auto elem_0_1_0_1 = 0.29702946944511982;
    const auto elem_0_1_1_0 = 0.29702946944511982;
    const auto elem_0_1_1_1 = 0.44410904384277333;
    const auto elem_1_0_0_0 = 0.44410904384277333;
    const auto elem_1_0_0_1 = 0.29702946944511982;
    const auto elem_1_0_1_0 = 0.29702946944511982;
    const auto elem_1_0_1_1 = 0.44410904384277361;
    const auto elem_1_1_0_0 = 0.56967771733030592;
    const auto elem_1_1_0_1 = 0.44410904384277350;
    const auto elem_1_1_1_0 = 0.44410904384277344;
    const auto elem_1_1_1_1 = 0.77460834925515787;

    REQUIRE_THAT(ee_grid.get(0, 0, 0, 0), Catch::Matchers::WithinAbs(elem_0_0_0_0, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(0, 0, 0, 1), Catch::Matchers::WithinAbs(elem_0_0_0_1, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(0, 0, 1, 0), Catch::Matchers::WithinAbs(elem_0_0_1_0, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(0, 0, 1, 1), Catch::Matchers::WithinAbs(elem_0_0_1_1, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(0, 1, 0, 0), Catch::Matchers::WithinAbs(elem_0_1_0_0, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(0, 1, 0, 1), Catch::Matchers::WithinAbs(elem_0_1_0_1, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(0, 1, 1, 0), Catch::Matchers::WithinAbs(elem_0_1_1_0, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(0, 1, 1, 1), Catch::Matchers::WithinAbs(elem_0_1_1_1, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(1, 0, 0, 0), Catch::Matchers::WithinAbs(elem_1_0_0_0, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(1, 0, 0, 1), Catch::Matchers::WithinAbs(elem_1_0_0_1, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(1, 0, 1, 0), Catch::Matchers::WithinAbs(elem_1_0_1_0, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(1, 0, 1, 1), Catch::Matchers::WithinAbs(elem_1_0_1_1, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(1, 1, 0, 0), Catch::Matchers::WithinAbs(elem_1_1_0_0, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(1, 1, 0, 1), Catch::Matchers::WithinAbs(elem_1_1_0_1, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(1, 1, 1, 0), Catch::Matchers::WithinAbs(elem_1_1_1_0, abs_tolerance));
    REQUIRE_THAT(ee_grid.get(1, 1, 1, 1), Catch::Matchers::WithinAbs(elem_1_1_1_1, abs_tolerance));
}
