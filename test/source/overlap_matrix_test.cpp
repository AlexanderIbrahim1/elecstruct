// TODO: remove later
#include <iostream>

#include <cmath>
#include <cstdint>
#include <tuple>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Eigen/Dense>

#include "elecstruct/atoms.hpp"
#include "elecstruct/basis/basis.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/matrices.hpp"

constexpr auto ABS_TOLERANCE = double {1.0e-6};
constexpr auto REL_TOLERANCE = double {1.0e-6};

auto is_column_equal(const Eigen::VectorXd& vec0, const Eigen::VectorXd& vec1, double tolerance) -> bool
{
    if (vec0.size() != vec1.size()) {
        return false;
    }

    const auto size = vec0.size();

    for (Eigen::Index i {0}; i < size; ++i) {
        if (std::fabs(vec0(i) - vec1(i)) > tolerance) {
            return false;
        }
    }

    return true;
}

auto is_column_equal_within_sign(const Eigen::VectorXd& vec0, const Eigen::VectorXd& vec1, double tolerance) -> bool
{
    return is_column_equal(vec0, vec1, tolerance) || is_column_equal(vec0, -vec1, tolerance);
}

auto are_columns_equal_within_sign(
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

    for (std::size_t i {0}; i < vec0.size(); ++i) {
        if (!is_column_equal_within_sign(vec0[i], vec1[i], tolerance)) {
            return false;
        }
    }

    return true;
}

TEST_CASE("overlap matrix : 1s orbital")
{
    using AOL = elec::AtomicOrbitalLabel;

    const auto atoms = std::vector<elec::AtomInfo> {
        elec::AtomInfo {elec::AtomLabel::H, coord::Cartesian3D {0.0, 0.0, 0.0}, {AOL::S1}}
    };

    const auto basis = elec::create_atomic_orbitals_sto3g(atoms);

    const auto overlap_mtx = elec::overlap_matrix(basis);

    REQUIRE(overlap_mtx.cols() == 1);
    REQUIRE(overlap_mtx.rows() == 1);
    REQUIRE_THAT(overlap_mtx(0, 0), Catch::Matchers::WithinRel(1.0, REL_TOLERANCE));
}

TEST_CASE("transformation matrix")
{
    // this test performs a comparison with a result from the Python code
    const auto tolerance = 1.0e-5;

    const auto input = [&]()
    {
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

    // std::cout << expected_columns[0] << '\n';
    // std::cout << expected_columns[1] << '\n';
    // std::cout << actual_columns[0] << '\n';
    // std::cout << actual_columns[1] << '\n';

    REQUIRE(are_columns_equal_within_sign(expected_columns, actual_columns, tolerance));
}

TEST_CASE("transformation matrix applied to overlap matrix gives identity")
{
    auto overlap_mtx = Eigen::MatrixXd {3, 3};
    overlap_mtx << 1.0, 0.2, 0.3, 0.2, 1.0, 0.1, 0.3, 0.1, 1.0;

    const auto transform_mtx = elec::transformation_matrix(overlap_mtx);

    const auto result = transform_mtx.transpose() * overlap_mtx * transform_mtx;

    // check the diagonals
    REQUIRE_THAT(result(0, 0), Catch::Matchers::WithinAbs(1.0, ABS_TOLERANCE));
    REQUIRE_THAT(result(1, 1), Catch::Matchers::WithinAbs(1.0, ABS_TOLERANCE));
    REQUIRE_THAT(result(2, 2), Catch::Matchers::WithinAbs(1.0, ABS_TOLERANCE));

    // check the off-diagonals
    REQUIRE_THAT(result(1, 0), Catch::Matchers::WithinAbs(0.0, ABS_TOLERANCE));
    REQUIRE_THAT(result(2, 0), Catch::Matchers::WithinAbs(0.0, ABS_TOLERANCE));
    REQUIRE_THAT(result(0, 1), Catch::Matchers::WithinAbs(0.0, ABS_TOLERANCE));
    REQUIRE_THAT(result(2, 1), Catch::Matchers::WithinAbs(0.0, ABS_TOLERANCE));
    REQUIRE_THAT(result(0, 2), Catch::Matchers::WithinAbs(0.0, ABS_TOLERANCE));
    REQUIRE_THAT(result(1, 2), Catch::Matchers::WithinAbs(0.0, ABS_TOLERANCE));
}

TEST_CASE("example overlap matrix to transformation matrix")
{
    const auto tolerance = 1.0e-3;

    /*
        example taken for water from online
        SOURCE: https://www.han-sur-lesse-winterschool.nl/downloads/2021/slides_filot.pdf
    */

    auto overlap_mtx = Eigen::MatrixXd {7, 7};

    overlap_mtx << 1.0000, 0.2367, 0.0000, 0.0000, 0.0000, 0.1584, 0.1584, 0.2367, 1.0000, 0.0000, 0.0000, 0.0000,
        0.8098, 0.8098, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 0.3714, -0.3714, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000,
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, -0.2322, -0.2322, 0.1584, 0.8098, 0.3714, 0.0000,
        -0.2322, 1.0000, 0.6158, 0.1584, 0.8098, -0.3714, 0.0000, -0.2322, 0.6158, 1.0000;

    auto expected_transform_mtx = Eigen::MatrixXd {7, 7};
    expected_transform_mtx << -0.1769, -0.0000, 0.8737, -0.0000, -0.5128, -0.0000, 0.1205, 2.7839, 0.0000, -0.2202,
        -0.0000, -0.1818, -0.0000, 0.3609, 0.0000, -1.7232, 0.0000, 0.0000, -0.0000, 0.7607, 0.0000, -0.0000, 0.0000,
        -0.0000, -1.0000, 0.0000, 0.0000, 0.0000, -0.7851, -0.0000, -0.5239, -0.0000, -0.8095, -0.0000, -0.0983,
        -1.5636, 2.1266, -0.1140, 0.0000, 0.0704, 0.3082, 0.3391, -1.5636, -2.1266, -0.1140, -0.0000, 0.0704, -0.3082,
        0.3391;

    const auto actual_transform_mtx = elec::transformation_matrix(overlap_mtx);

    for (Eigen::Index i {0}; i < 7; ++i) {
        REQUIRE(is_column_equal_within_sign(expected_transform_mtx.col(i), actual_transform_mtx.col(i), tolerance));
    }
}
