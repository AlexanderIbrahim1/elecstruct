#include <cmath>
#include <cstdint>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/integrals/integrals.hpp"
#include "elecstruct/integrals/overlap_integrals.hpp"
#include "elecstruct/orbitals.hpp"

constexpr auto ABS_TOL = double {1.0e-8};

TEST_CASE("gaussian product")
{
    SECTION("example from desmos")
    {
        const auto centre0 = coord::Cartesian3D {1.0, 0.0, 0.0};
        const auto centre1 = coord::Cartesian3D {-1.0, 0.0, 0.0};
        const auto info0 = elec::GaussianInfo {2.0, 1.0 / 2.0};
        const auto info1 = elec::GaussianInfo {3.0, 1.0 / 3.0};

        const auto expected_centre = coord::Cartesian3D {0.2, 0.0, 0.0};
        const auto expected_coefficient = double {2.6959737847};
        const auto expected_exponent = info0.exponent + info1.exponent;

        const auto [centre, gauss] = elec::gaussian_product(centre0, centre1, info0, info1);

        REQUIRE(coord::almost_equals(centre, expected_centre));
        REQUIRE_THAT(gauss.coefficient, Catch::Matchers::WithinAbs(expected_coefficient, ABS_TOL));
        REQUIRE_THAT(gauss.exponent, Catch::Matchers::WithinAbs(expected_exponent, ABS_TOL));
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

        const auto actual_output = elec::gaussian_norm(pair.input_ang_mom, pair.input_exponent);
        REQUIRE_THAT(actual_output, Catch::Matchers::WithinRel(pair.expected_output));
    }
}

TEST_CASE("overlap integrals")
{
    // Test case results are compared with the values provided in the body of
    //
    // journal: The Mathematica Journal
    // title:   Evaluation of Gaussian Molecular Integrals I. Overlap Integrals
    // authors: Minhhuy Hô and Julio Manuel Hernández-Pérez

    SECTION("s-function with itself")
    {
        const auto exponent = GENERATE(1.0, 1.5, 2.0);
        const auto coefficient = GENERATE(1.0, 2.0, 3.0);

        const auto centre0 = coord::Cartesian3D {0.0, 0.0, 0.0};
        const auto centre1 = coord::Cartesian3D {0.0, 0.0, 0.0};
        const auto gauss0 = elec::GaussianInfo {coefficient, exponent};
        const auto gauss1 = elec::GaussianInfo {coefficient, exponent};
        const auto angmom0 = elec::AngularMomentumNumbers {0, 0, 0};
        const auto angmom1 = elec::AngularMomentumNumbers {0, 0, 0};

        const auto [new_centre, new_info] = elec::gaussian_product(centre0, centre1, gauss0, gauss1);

        // clang-format off
        const auto overlap_x = elec::unnormalized_overlap_integral_1d(
            {angmom0.x, gauss0.exponent, centre0.x},
            {angmom1.x, gauss1.exponent, centre1.x},
            new_centre.x
        );

        const auto overlap_y = elec::unnormalized_overlap_integral_1d(
            {angmom0.y, gauss0.exponent, centre0.y},
            {angmom1.y, gauss1.exponent, centre1.y},
            new_centre.y
        );

        const auto overlap_z = elec::unnormalized_overlap_integral_1d(
            {angmom0.z, gauss0.exponent, centre0.z},
            {angmom1.z, gauss1.exponent, centre1.z},
            new_centre.z
        );
        // clang-format on

        const auto actual_total_overlap = new_info.coefficient * overlap_x * overlap_y * overlap_z;

        const auto expected_overlap_1d = std::sqrt(M_PI / (gauss0.exponent + gauss1.exponent));
        const auto expected_coefficient = gauss0.coefficient * gauss1.coefficient;
        const auto expected_total_overlap = expected_coefficient * std::pow(expected_overlap_1d, 3);

        REQUIRE_THAT(new_info.coefficient, Catch::Matchers::WithinRel(expected_coefficient));
        REQUIRE_THAT(overlap_x, Catch::Matchers::WithinRel(expected_overlap_1d));
        REQUIRE_THAT(overlap_y, Catch::Matchers::WithinRel(expected_overlap_1d));
        REQUIRE_THAT(overlap_z, Catch::Matchers::WithinRel(expected_overlap_1d));
        REQUIRE_THAT(actual_total_overlap, Catch::Matchers::WithinRel(expected_total_overlap));
    }

    SECTION("p-function with itself")
    {
        const auto exponent0 = GENERATE(1.0, 1.5);
        const auto exponent1 = GENERATE(1.25, 1.75);
        const auto coefficient0 = GENERATE(1.0, 2.0);
        const auto coefficient1 = GENERATE(1.5, 2.5);

        const auto centre0 = coord::Cartesian3D {0.0, 0.0, 0.0};
        const auto centre1 = coord::Cartesian3D {0.0, 0.0, 0.0};
        const auto gauss0 = elec::GaussianInfo {coefficient0, exponent0};
        const auto gauss1 = elec::GaussianInfo {coefficient1, exponent1};
        const auto angmom0 = elec::AngularMomentumNumbers {1, 0, 0};
        const auto angmom1 = elec::AngularMomentumNumbers {1, 0, 0};

        const auto [new_centre, new_info] = elec::gaussian_product(centre0, centre1, gauss0, gauss1);

        // clang-format off
        const auto overlap_x = elec::unnormalized_overlap_integral_1d(
            {angmom0.x, gauss0.exponent, centre0.x},
            {angmom1.x, gauss1.exponent, centre1.x},
            new_centre.x
        );

        const auto overlap_y = elec::unnormalized_overlap_integral_1d(
            {angmom0.y, gauss0.exponent, centre0.y},
            {angmom1.y, gauss1.exponent, centre1.y},
            new_centre.y
        );

        const auto overlap_z = elec::unnormalized_overlap_integral_1d(
            {angmom0.z, gauss0.exponent, centre0.z},
            {angmom1.z, gauss1.exponent, centre1.z},
            new_centre.z
        );
        // clang-format on

        const auto actual_total_overlap = new_info.coefficient * overlap_x * overlap_y * overlap_z;

        const auto expected_overlap_x = [&]()
        {
            const auto sum_expon = gauss0.exponent + gauss1.exponent;
            const auto overlap_s = std::sqrt(M_PI / sum_expon);
            const auto overlap_due_to_p = 1.0 / (2.0 * sum_expon);

            return overlap_s * overlap_due_to_p;
        }();

        const auto expected_overlap_y = std::sqrt(M_PI / (gauss0.exponent + gauss1.exponent));
        const auto expected_overlap_z = std::sqrt(M_PI / (gauss0.exponent + gauss1.exponent));

        const auto expected_coefficient = gauss0.coefficient * gauss1.coefficient;
        const auto expected_total_overlap = [&]()
        {
            const auto expected_overlap_xyz = expected_overlap_x * expected_overlap_y * expected_overlap_z;
            return expected_coefficient * expected_overlap_xyz;
        }();

        REQUIRE_THAT(new_info.coefficient, Catch::Matchers::WithinRel(expected_coefficient));
        REQUIRE_THAT(overlap_x, Catch::Matchers::WithinRel(expected_overlap_x));
        REQUIRE_THAT(overlap_y, Catch::Matchers::WithinRel(expected_overlap_y));
        REQUIRE_THAT(overlap_z, Catch::Matchers::WithinRel(expected_overlap_z));
        REQUIRE_THAT(actual_total_overlap, Catch::Matchers::WithinRel(expected_total_overlap));
    }
}
