#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <Eigen/Dense>

#include "elecstruct/restricted_hartree_fock/step.hpp"

constexpr auto ABS_TOLERANCE = double {1.0e-6};

TEST_CASE("sorted_indices")
{
    auto input = Eigen::VectorXd(5);
    input << 3.0, 1.0, 0.0, 2.0, 4.0;
    const auto sorted = elec::indices_to_sort(input);

    const auto expected = std::vector<Eigen::Index> {2, 1, 3, 0, 4};

    REQUIRE_THAT(sorted, Catch::Matchers::Equals(expected));
}

TEST_CASE("matrix_with_sorted_columns")
{
    const auto indices = std::vector<Eigen::Index> {2, 1, 0};

    auto input = Eigen::MatrixXd {3, 3};
    input << 0.0, 0.1, 0.2, 1.0, 1.1, 1.2, 2.0, 2.1, 2.2;

    const auto output = elec::matrix_with_sorted_columns(input, indices);

    REQUIRE_THAT(output(0, 0), Catch::Matchers::WithinAbs(0.2, ABS_TOLERANCE));
    REQUIRE_THAT(output(1, 0), Catch::Matchers::WithinAbs(1.2, ABS_TOLERANCE));
    REQUIRE_THAT(output(2, 0), Catch::Matchers::WithinAbs(2.2, ABS_TOLERANCE));
    REQUIRE_THAT(output(0, 1), Catch::Matchers::WithinAbs(0.1, ABS_TOLERANCE));
    REQUIRE_THAT(output(1, 1), Catch::Matchers::WithinAbs(1.1, ABS_TOLERANCE));
    REQUIRE_THAT(output(2, 1), Catch::Matchers::WithinAbs(2.1, ABS_TOLERANCE));
    REQUIRE_THAT(output(0, 2), Catch::Matchers::WithinAbs(0.0, ABS_TOLERANCE));
    REQUIRE_THAT(output(1, 2), Catch::Matchers::WithinAbs(1.0, ABS_TOLERANCE));
    REQUIRE_THAT(output(2, 2), Catch::Matchers::WithinAbs(2.0, ABS_TOLERANCE));
}
