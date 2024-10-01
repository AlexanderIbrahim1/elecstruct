#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "elecstruct/cartesian3d.hpp"

TEST_CASE("Cartesian3D operations", "[Cartesian3D]")
{
    auto point1 = coord::Cartesian3D {1.0, 2.0, 3.0};
    auto point2 = coord::Cartesian3D {4.0, 5.0, 6.0};

    SECTION("Addition")
    {
        const auto result = point1 + point2;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(5.0));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(7.0));
        REQUIRE_THAT(result.z, Catch::Matchers::WithinRel(9.0));

        point1 += point2;
        REQUIRE_THAT(point1.x, Catch::Matchers::WithinRel(5.0));
        REQUIRE_THAT(point1.y, Catch::Matchers::WithinRel(7.0));
        REQUIRE_THAT(point1.z, Catch::Matchers::WithinRel(9.0));
    }

    SECTION("Subtraction")
    {
        const auto result = point1 - point2;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(-3.0));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(-3.0));
        REQUIRE_THAT(result.z, Catch::Matchers::WithinRel(-3.0));

        point1 -= point2;
        REQUIRE_THAT(point1.x, Catch::Matchers::WithinRel(-3.0));
        REQUIRE_THAT(point1.y, Catch::Matchers::WithinRel(-3.0));
        REQUIRE_THAT(point1.z, Catch::Matchers::WithinRel(-3.0));
    }

    SECTION("Multiplication with integers")
    {
        const auto result = 2 * point1;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(2.0));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(4.0));
        REQUIRE_THAT(result.z, Catch::Matchers::WithinRel(6.0));

        point1 *= 3;
        REQUIRE_THAT(point1.x, Catch::Matchers::WithinRel(3.0));
        REQUIRE_THAT(point1.y, Catch::Matchers::WithinRel(6.0));
        REQUIRE_THAT(point1.z, Catch::Matchers::WithinRel(9.0));
    }

    SECTION("Multiplication with floats")
    {
        const auto result = 1.5f * point1;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(1.5));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(3.0));
        REQUIRE_THAT(result.z, Catch::Matchers::WithinRel(4.5));

        point1 *= 0.5f;
        REQUIRE_THAT(point1.x, Catch::Matchers::WithinRel(0.5));
        REQUIRE_THAT(point1.y, Catch::Matchers::WithinRel(1.0));
        REQUIRE_THAT(point1.z, Catch::Matchers::WithinRel(1.5));
    }

    SECTION("Multiplication with doubles")
    {
        const auto result = 2.5 * point1;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(2.5));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(5.0));
        REQUIRE_THAT(result.z, Catch::Matchers::WithinRel(7.5));

        point1 *= 4.0;
        REQUIRE_THAT(point1.x, Catch::Matchers::WithinRel(4.0));
        REQUIRE_THAT(point1.y, Catch::Matchers::WithinRel(8.0));
        REQUIRE_THAT(point1.z, Catch::Matchers::WithinRel(12.0));
    }

    SECTION("Division with non-zero integers")
    {
        const auto result = point1 / 4;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(0.25));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(0.50));
        REQUIRE_THAT(result.z, Catch::Matchers::WithinRel(0.75));

        point1 /= 2;
        REQUIRE_THAT(point1.x, Catch::Matchers::WithinRel(0.5));
        REQUIRE_THAT(point1.y, Catch::Matchers::WithinRel(1.0));
        REQUIRE_THAT(point1.z, Catch::Matchers::WithinRel(1.5));
    }

    SECTION("Division with non-zero floats")
    {
        const auto result = point1 / 2.0f;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(0.5));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(1.0));
        REQUIRE_THAT(result.z, Catch::Matchers::WithinRel(1.5));

        point1 /= 2.0f;
        REQUIRE_THAT(point1.x, Catch::Matchers::WithinRel(0.5));
        REQUIRE_THAT(point1.y, Catch::Matchers::WithinRel(1.0));
        REQUIRE_THAT(point1.z, Catch::Matchers::WithinRel(1.5));
    }

    SECTION("Division with non-zero doubles")
    {
        const auto result = point1 / 0.5;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(2.0));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(4.0));
        REQUIRE_THAT(result.z, Catch::Matchers::WithinRel(6.0));

        point1 /= 4.0;
        REQUIRE_THAT(point1.x, Catch::Matchers::WithinRel(0.25));
        REQUIRE_THAT(point1.y, Catch::Matchers::WithinRel(0.50));
        REQUIRE_THAT(point1.z, Catch::Matchers::WithinRel(0.75));
    }
}
