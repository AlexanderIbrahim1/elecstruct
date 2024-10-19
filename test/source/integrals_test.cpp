#include <cstdint>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/basis/gaussian_info.hpp"
#include "elecstruct/integrals/integrals.hpp"

constexpr auto ABS_TOL = double {1.0e-8};

TEST_CASE("gaussian product")
{
    SECTION("example from desmos")
    {
        const auto centre0 = coord::Cartesian3D {1.0, 0.0, 0.0};
        const auto centre1 = coord::Cartesian3D {-1.0, 0.0, 0.0};
        const auto info0 = elec::GaussianInfo {2.0, 1.0/2.0};
        const auto info1 = elec::GaussianInfo {3.0, 1.0/3.0};

        const auto expected_centre = coord::Cartesian3D {0.2, 0.0, 0.0};
        const auto expected_coefficient = double {2.6959737847};
        const auto expected_exponent = info0.exponent + info1.exponent;

        const auto [centre, gauss] = elec::gaussian_product(centre0, centre1, info0, info1);

        REQUIRE(coord::almost_equals(centre, expected_centre));
        REQUIRE_THAT(gauss.coefficient, Catch::Matchers::WithinAbs(expected_coefficient, ABS_TOL));
        REQUIRE_THAT(gauss.exponent, Catch::Matchers::WithinAbs(expected_exponent, ABS_TOL));
    }
}

TEST_CASE("double factorial")
{
    struct TestPair
    {
        std::uint64_t input;
        std::uint64_t expected_output;
    };

    auto pair = GENERATE(
        TestPair {0, 1},
        TestPair {1, 1},
        TestPair {2, 2},
        TestPair {3, 3 * 1},
        TestPair {4, 4 * 2},
        TestPair {5, 5 * 3},
        TestPair {6, 6 * 4 * 2},
        TestPair {7, 7 * 5 * 3 * 1},
        TestPair {8, 8 * 6 * 4 * 2},
        TestPair {9, 9 * 7 * 5 * 3 * 1},
        TestPair {10, 10 * 8 * 6 * 4 * 2}
    );

    const auto actual_output = impl_elec::double_factorial(pair.input);
    REQUIRE(actual_output == pair.expected_output);
}
