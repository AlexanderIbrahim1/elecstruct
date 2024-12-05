#include <cmath>
#include <cstdint>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "elecstruct/mathtools/misc.hpp"

constexpr auto ABS_TOL = double {1.0e-8};

TEST_CASE("double factorial")
{
    struct TestPair
    {
        std::int64_t input;
        std::int64_t expected_output;
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

    const auto actual_output = elec::math::double_factorial(pair.input);
    REQUIRE(actual_output == pair.expected_output);
}
