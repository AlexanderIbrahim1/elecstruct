#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/geometry.hpp"

constexpr auto ABS_TOL = double {1.0e-8};

TEST_CASE("Dot Product")
{
    SECTION("Dot product of parallel vectors")
    {
        const auto point0 = coord::Cartesian3D {1.0, 2.0, 3.0};
        const auto point1 = coord::Cartesian3D {2.0, 4.0, 6.0};
        const auto result = coord::dot_product(point0, point1);
        REQUIRE_THAT(result, Catch::Matchers::WithinRel(28.0));  // 1*2 + 2*4 + 3*6 = 28
    }

    SECTION("Dot product of orthogonal vectors")
    {
        const auto point0 = coord::Cartesian3D {1.0, 0.0, 0.0};
        const auto point1 = coord::Cartesian3D {0.0, 1.0, 0.0};
        const auto result = coord::dot_product(point0, point1);
        REQUIRE_THAT(result, Catch::Matchers::WithinRel(0.0));  // 1*0 + 0*1 + 0*0 = 0
    }

    SECTION("Dot product of zero vector")
    {
        const auto point0 = coord::Cartesian3D {0.0, 0.0, 0.0};
        const auto point1 = coord::Cartesian3D {1.0, 2.0, 3.0};
        const auto result = coord::dot_product(point0, point1);
        REQUIRE_THAT(result, Catch::Matchers::WithinRel(0.0));  // 0*1 + 0*2 + 0*3 = 0
    }
}

TEST_CASE("Cross Product")
{
    SECTION("Cross product of parallel vectors")
    {
        const auto point0 = coord::Cartesian3D {1.0, 2.0, 3.0};
        const auto point1 = coord::Cartesian3D {2.0, 4.0, 6.0};
        const auto result = coord::cross_product(point0, point1);
        REQUIRE(coord::almost_equals(result, coord::Cartesian3D {0.0, 0.0, 0.0}, ABS_TOL));
    }

    SECTION("Cross product of orthogonal vectors")
    {
        const auto point0 = coord::Cartesian3D {1.0, 0.0, 0.0};
        const auto point1 = coord::Cartesian3D {0.0, 1.0, 0.0};
        const auto result = coord::cross_product(point0, point1);
        REQUIRE(coord::almost_equals(result, coord::Cartesian3D {0.0, 0.0, 1.0}, ABS_TOL));
    }

    SECTION("Cross product of arbitrary vectors")
    {
        const auto point0 = coord::Cartesian3D {1.0, 2.0, 3.0};
        const auto point1 = coord::Cartesian3D {4.0, 5.0, 6.0};
        const auto result = coord::cross_product(point0, point1);
        REQUIRE(coord::almost_equals(result, coord::Cartesian3D {-3.0, 6.0, -3.0}, ABS_TOL));
    }
}
