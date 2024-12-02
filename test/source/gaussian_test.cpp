#include <cmath>
#include <cstdint>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/geometry.hpp"
#include "elecstruct/mathtools/gaussian.hpp"
#include "elecstruct/integrals/overlap_integrals.hpp"
#include "elecstruct/orbitals.hpp"

constexpr auto ABS_TOL = double {1.0e-8};

TEST_CASE("gaussian product")
{
    SECTION("example from desmos")
    {
        const auto centre0 = coord::Cartesian3D {1.0, 0.0, 0.0};
        const auto centre1 = coord::Cartesian3D {-1.0, 0.0, 0.0};
        const auto exponent0 = double {1.0 / 2.0};
        const auto exponent1 = double {1.0 / 3.0};

        const auto expected_centre = coord::Cartesian3D {0.2, 0.0, 0.0};
        const auto expected_coefficient = double {0.449328964117}; // from DESMOS

        const auto [centre, coeff] = elec::math::gaussian_product(centre0, centre1, exponent0, exponent1);

        REQUIRE(coord::almost_equals(centre, expected_centre, ABS_TOL));
        REQUIRE_THAT(coeff, Catch::Matchers::WithinAbs(expected_coefficient, ABS_TOL));
    }
}

TEST_CASE("gaussian norm")
{
    struct TestPair
    {
        elec::AngularMomentumNumbers input_ang_mom;
        double input_exponent;
        double expected_output;
    };

    SECTION("s-orbital norms")
    {
        constexpr auto ang_mom_s = elec::AngularMomentumNumbers {0, 0, 0};

        auto pair = GENERATE_COPY(
            TestPair {ang_mom_s, 1.0, std::pow(2.0 * 1.0 / M_PI, 3.0 / 4.0)},
            TestPair {ang_mom_s, 2.0, std::pow(2.0 * 2.0 / M_PI, 3.0 / 4.0)},
            TestPair {ang_mom_s, 3.0, std::pow(2.0 * 3.0 / M_PI, 3.0 / 4.0)},
            TestPair {ang_mom_s, 0.5, std::pow(2.0 * 0.5 / M_PI, 3.0 / 4.0)}
        );

        const auto actual_output = elec::math::gaussian_norm(pair.input_ang_mom, pair.input_exponent);
        REQUIRE_THAT(actual_output, Catch::Matchers::WithinRel(pair.expected_output));
    }
}

TEST_CASE("unnormalized overlap integrals")
{
    // Test case results are compared with the values provided in the body of
    //
    // journal: The Mathematica Journal
    // title:   Evaluation of Gaussian Molecular Integrals I. Overlap Integrals
    // authors: Minhhuy Hô and Julio Manuel Hernández-Pérez

    SECTION("s-function with itself")
    {
        // the point of including these is to show that they should not affect the final result
        // - but maybe I should remove them, since we're running 16 unit tests here?
        const auto exponent0 = GENERATE(1.0, 1.5);
        const auto exponent1 = GENERATE(1.25, 1.75);

        const auto centre0 = coord::Cartesian3D {0.0, 0.0, 0.0};
        const auto centre1 = coord::Cartesian3D {0.0, 0.0, 0.0};
        const auto angmom0 = elec::AngularMomentumNumbers {0, 0, 0};
        const auto angmom1 = elec::AngularMomentumNumbers {0, 0, 0};

        [[maybe_unused]] const auto [new_centre, new_info] = elec::math::gaussian_product(centre0, centre1, exponent0, exponent1);

        // clang-format off
        const auto unorm_overlap_x = elec::unnormalized_overlap_integral_1d(
            {angmom0.x, exponent0, centre0.x},
            {angmom1.x, exponent1, centre1.x},
            new_centre.x
        );

        const auto unorm_overlap_y = elec::unnormalized_overlap_integral_1d(
            {angmom0.y, exponent0, centre0.y},
            {angmom1.y, exponent1, centre1.y},
            new_centre.y
        );

        const auto unorm_overlap_z = elec::unnormalized_overlap_integral_1d(
            {angmom0.z, exponent0, centre0.z},
            {angmom1.z, exponent1, centre1.z},
            new_centre.z
        );
        // clang-format on

        const auto expected_unorm_overlap_s = 1.0;
        REQUIRE_THAT(unorm_overlap_x, Catch::Matchers::WithinRel(expected_unorm_overlap_s));
        REQUIRE_THAT(unorm_overlap_y, Catch::Matchers::WithinRel(expected_unorm_overlap_s));
        REQUIRE_THAT(unorm_overlap_z, Catch::Matchers::WithinRel(expected_unorm_overlap_s));
    }

    SECTION("p-function with itself")
    {
        // the point of including these is to show that they should not affect the final result
        // - but maybe I should remove them, since we're running 16 unit tests here?
        const auto exponent0 = GENERATE(1.0, 1.5);
        const auto exponent1 = GENERATE(1.25, 1.75);

        const auto centre0 = coord::Cartesian3D {0.0, 0.0, 0.0};
        const auto centre1 = coord::Cartesian3D {0.0, 0.0, 0.0};
        const auto angmom0 = elec::AngularMomentumNumbers {1, 0, 0};
        const auto angmom1 = elec::AngularMomentumNumbers {1, 0, 0};

        const auto [new_centre, new_info] = elec::math::gaussian_product(centre0, centre1, exponent0, exponent1);

        // clang-format off
        const auto unorm_overlap_x = elec::unnormalized_overlap_integral_1d(
            {angmom0.x, exponent0, centre0.x},
            {angmom1.x, exponent1, centre1.x},
            new_centre.x
        );

        const auto unorm_overlap_y = elec::unnormalized_overlap_integral_1d(
            {angmom0.y, exponent0, centre0.y},
            {angmom1.y, exponent1, centre1.y},
            new_centre.y
        );

        const auto unorm_overlap_z = elec::unnormalized_overlap_integral_1d(
            {angmom0.z, exponent0, centre0.z},
            {angmom1.z, exponent1, centre1.z},
            new_centre.z
        );
        // clang-format on

        const auto expected_unorm_overlap_x = 1.0 / (2.0 * (exponent0 + exponent1));
        const auto expected_unorm_overlap_y = 1.0;
        const auto expected_unorm_overlap_z = 1.0;

        REQUIRE_THAT(unorm_overlap_x, Catch::Matchers::WithinRel(expected_unorm_overlap_x));
        REQUIRE_THAT(unorm_overlap_y, Catch::Matchers::WithinRel(expected_unorm_overlap_y));
        REQUIRE_THAT(unorm_overlap_z, Catch::Matchers::WithinRel(expected_unorm_overlap_z));
    }
}
