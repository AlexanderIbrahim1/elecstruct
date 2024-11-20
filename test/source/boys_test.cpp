#include <cstdint>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "elecstruct/integrals/boys.hpp"


TEST_CASE("Boys function comparison with Python code")
{
    /*
        NOTE: all the expected outputs were created by a Python function
        where the comparison against the true Boys function was easier
        to verify
    */

    constexpr auto REL_TOLERANCE = double {1.0e-8};

    const auto n_terms_small = std::int64_t {11};

    struct TestPair
    {
        double input;
        double expected;
    };

    SECTION("order = 0")
    {
        const auto order = std::int64_t {0};

        const auto pair = GENERATE(
            TestPair {0.01, 0.9966766429033636},
            TestPair {1.00, 0.7468241338237177},
            TestPair {2.00, 0.5981459383100741},
            TestPair {5.00, 0.3963327297606011},
            TestPair {10.00, 0.2802495608198964}
        );

        const auto actual = elec::boys_mix_small_large(pair.input, order, n_terms_small);
        REQUIRE_THAT(actual, Catch::Matchers::WithinRel(pair.expected, REL_TOLERANCE));
    }

    SECTION("order = 1")
    {
        const auto order = std::int64_t {1};

        const auto pair = GENERATE(
            TestPair {0.01, 0.3313404577097725},
            TestPair {1.00, 0.1894723467504448},
            TestPair {2.00, 0.1157039563971642},
            TestPair {5.00, 0.0396332729760601},
            TestPair {10.00, 0.0140124780409948}
        );

        const auto actual = elec::boys_mix_small_large(pair.input, order, n_terms_small);
        REQUIRE_THAT(actual, Catch::Matchers::WithinRel(pair.expected, REL_TOLERANCE));
    }
}
