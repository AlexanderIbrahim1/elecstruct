#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/integrals/kinetic_integrals.hpp"
#include "elecstruct/orbitals.hpp"

TEST_CASE("cyclic shifts")
{
    constexpr auto REL_TOLERANCE = double {1.0e-8};

    SECTION("AngularMomentumNumbers")
    {
        const auto angmoms = elec::AngularMomentumNumbers {1, 2, 3};

        const auto directed_orig = elec::DirectedAngularMomentumNumbers {angmoms};
        const auto directed_left = elec::left_cyclic_shift(directed_orig);
        const auto directed_right = elec::right_cyclic_shift(directed_orig);

        REQUIRE(directed_orig.main == angmoms.x);
        REQUIRE(directed_orig.other0 == angmoms.y);
        REQUIRE(directed_orig.other1 == angmoms.z);

        REQUIRE(directed_left.main == angmoms.y);
        REQUIRE(directed_left.other0 == angmoms.z);
        REQUIRE(directed_left.other1 == angmoms.x);

        REQUIRE(directed_right.main == angmoms.z);
        REQUIRE(directed_right.other0 == angmoms.x);
        REQUIRE(directed_right.other1 == angmoms.y);
    }

    SECTION("Cartesian3D")
    {
        const auto positions = coord::Cartesian3D {1.0, 2.0, 3.0};

        const auto directed_orig = elec::DirectedCartesian3D {positions};
        const auto directed_left = elec::left_cyclic_shift(directed_orig);
        const auto directed_right = elec::right_cyclic_shift(directed_orig);

        REQUIRE_THAT(directed_orig.main, Catch::Matchers::WithinRel(positions.x, REL_TOLERANCE));
        REQUIRE_THAT(directed_orig.other0, Catch::Matchers::WithinRel(positions.y, REL_TOLERANCE));
        REQUIRE_THAT(directed_orig.other1, Catch::Matchers::WithinRel(positions.z, REL_TOLERANCE));

        REQUIRE_THAT(directed_left.main, Catch::Matchers::WithinRel(positions.y, REL_TOLERANCE));
        REQUIRE_THAT(directed_left.other0, Catch::Matchers::WithinRel(positions.z, REL_TOLERANCE));
        REQUIRE_THAT(directed_left.other1, Catch::Matchers::WithinRel(positions.x, REL_TOLERANCE));

        REQUIRE_THAT(directed_right.main, Catch::Matchers::WithinRel(positions.z, REL_TOLERANCE));
        REQUIRE_THAT(directed_right.other0, Catch::Matchers::WithinRel(positions.x, REL_TOLERANCE));
        REQUIRE_THAT(directed_right.other1, Catch::Matchers::WithinRel(positions.y, REL_TOLERANCE));
    }
}
