#include <cmath>
#include <cstdint>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Eigen/Dense>

#include "elecstruct/matrices.hpp"


auto approx_equal(const Eigen::Vector2d& eig0, const Eigen::Vector2d& eig1, double tolerance) -> bool
{
    const auto diff0 = eig0(0) - eig1(0);
    const auto diff1 = eig0(1) - eig1(1);

    return std::sqrt(diff0 * diff0 + diff1 * diff1) < tolerance;
}


auto are_columns_equal(
    const std::vector<Eigen::Vector2d>& vec0,
    const std::vector<Eigen::Vector2d>& vec1,
    double tolerance
) -> bool
{
    // putting elements into a set, or performing a sort before comparing, gets a bit finicky
    // with the floating-point numbers; so I perform an O(n^2) check instead
    if (vec0.size() != vec1.size()) {
        return false;
    }

    const auto size = vec0.size();

    auto matches = std::vector<std::uint8_t> (size, 0);

    for (std::size_t i0 {0}; i0 < size; ++i0) {
        bool flag_found = false;
        for (std::size_t i1 {0}; i1 < size; ++i1) {
            if (approx_equal(vec0[i0], vec1[i1], tolerance) && matches[i1] != 1) {
                matches[i1] = 1;
                flag_found = true;
                break;
            }
        }

        if (!flag_found) {
            return false;
        }
    }

    return true;
}


TEST_CASE("transformation matrix")
{
    // this test performs a comparison with a result from the Python code
    const auto tolerance = 1.0e-5;

    const auto input = [&]() {
        auto input_ = Eigen::MatrixXd {2, 2};
        input_(0, 0) = 5.0;
        input_(1, 0) = 1.0;
        input_(0, 1) = 1.0;
        input_(1, 1) = 4.0;

        return input_;
    }();

    const auto expected_columns = std::vector<Eigen::Vector2d> {
        {-0.2858769, 0.46255854},
        {0.35888817, 0.22180508}
    };

    const auto actual_output = elec::transformation_matrix(input);
    const auto actual_columns = std::vector<Eigen::Vector2d> {
        {actual_output(0, 0), actual_output(1, 0)},
        {actual_output(0, 1), actual_output(1, 1)}
    };

    REQUIRE(are_columns_equal(expected_columns, actual_columns, tolerance));
}
