#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

#include "elecstruct/matrices/overlap.hpp"


TEST_CASE("transformation matrix")
{
    // this test performs a comparison with a result from the Python code
    const auto tolerance = 1.0e-5;

    const auto input = [&]() {
        auto input = Eigen::MatrixXd {2, 2};
        input << 5, 1, 1, 4;

        return input;
    }();

    const auto expected_output = [&]() {
        // values taken directly from Python code
        auto output = Eigen::MatrixXd {2, 2};
        // output << 0.35888817, -0.22180508, 0.2858769, 0.46255854;
        output << -0.2858769, -0.46255854, -0.35888817, 0.22180508;

        return output;
    }();

    const auto actual_output = elec::transformation_matrix(input);

    REQUIRE_THAT(actual_output(0, 0), Catch::Matchers::WithinAbs(expected_output(0, 0), tolerance));
    REQUIRE_THAT(actual_output(0, 1), Catch::Matchers::WithinAbs(expected_output(0, 1), tolerance));
    REQUIRE_THAT(actual_output(1, 0), Catch::Matchers::WithinAbs(expected_output(1, 0), tolerance));
    REQUIRE_THAT(actual_output(1, 1), Catch::Matchers::WithinAbs(expected_output(1, 1), tolerance));
}
