#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "elecstruct/grids/grid4d.hpp"

TEST_CASE("Grid4D construction with sizes")
{
    SECTION("Create a 2x2x2x2 Grid4D and verify sizes")
    {
        const auto grid = elec::grid::Grid4D {2, 2, 2, 2};
        const auto [size0, size1, size2, size3] = grid.sizes();
        REQUIRE(size0 == 2);
        REQUIRE(size1 == 2);
        REQUIRE(size2 == 2);
        REQUIRE(size3 == 2);
    }

    SECTION("Create a 3x4x5x6 Grid4D and verify sizes")
    {
        const auto grid = elec::grid::Grid4D {3, 4, 5, 6};
        const auto [size0, size1, size2, size3] = grid.sizes();
        REQUIRE(size0 == 3);
        REQUIRE(size1 == 4);
        REQUIRE(size2 == 5);
        REQUIRE(size3 == 6);
    }
}

TEST_CASE("Grid4D set and get in bounds")
{
    auto grid = elec::grid::Grid4D {2, 2, 2, 2};

    SECTION("Set and get values in bounds using set() and get()")
    {
        grid.set(0, 0, 0, 0, 1.0);
        grid.set(1, 1, 1, 1, 4.0);
        REQUIRE_THAT(grid.get(0, 0, 0, 0), Catch::Matchers::WithinRel(1.0));
        REQUIRE_THAT(grid.get(1, 1, 1, 1), Catch::Matchers::WithinRel(4.0));
    }

    SECTION("Set and get values in bounds using set_at() and get_at()")
    {
        grid.set_at(0, 1, 1, 0, 2.0);
        grid.set_at(1, 0, 1, 1, 3.0);
        REQUIRE_THAT(grid.get_at(0, 1, 1, 0), Catch::Matchers::WithinRel(2.0));
        REQUIRE_THAT(grid.get_at(1, 0, 1, 1), Catch::Matchers::WithinRel(3.0));
    }
}

TEST_CASE("Grid4D constructed from std::vector<double>")
{
    // clang-format off
    const auto data = std::vector<double> {
        0.0,  1.0,  2.0,  3.0,  4.0,
        5.0,  6.0,  7.0,  8.0,  9.0,
        10.0, 11.0, 12.0, 13.0, 14.0,
        15.0, 16.0, 17.0, 18.0, 19.0,

        20.0, 21.0, 22.0, 23.0, 24.0,
        25.0, 26.0, 27.0, 28.0, 29.0,
        30.0, 31.0, 32.0, 33.0, 34.0,
        35.0, 36.0, 37.0, 38.0, 39.0,

        40.0, 41.0, 42.0, 43.0, 44.0,
        45.0, 46.0, 47.0, 48.0, 49.0,
        50.0, 51.0, 52.0, 53.0, 54.0,
        55.0, 56.0, 57.0, 58.0, 59.0,



        60.0, 61.0, 62.0, 63.0, 64.0,
        65.0, 66.0, 67.0, 68.0, 69.0,
        70.0, 71.0, 72.0, 73.0, 74.0,
        75.0, 76.0, 77.0, 78.0, 79.0,

        80.0, 81.0, 82.0, 83.0, 84.0,
        85.0, 86.0, 87.0, 88.0, 89.0,
        90.0, 91.0, 92.0, 93.0, 94.0,
        95.0, 96.0, 97.0, 98.0, 99.0,

        100.0, 101.0, 102.0, 103.0, 104.0,
        105.0, 106.0, 107.0, 108.0, 109.0,
        110.0, 111.0, 112.0, 113.0, 114.0,
        115.0, 116.0, 117.0, 118.0, 119.0
    };
    // clang-format on

    const auto grid = elec::grid::Grid4D {data, 2, 3, 4, 5};

    std::size_t index = 0;
    for (std::size_t i0 = 0; i0 < 2; ++i0) {
        for (std::size_t i1 = 0; i1 < 3; ++i1) {
            for (std::size_t i2 = 0; i2 < 4; ++i2) {
                for (std::size_t i3 = 0; i3 < 5; ++i3) {
                    REQUIRE_THAT(grid.get(i0, i1, i2, i3), Catch::Matchers::WithinRel(data[index]));
                    ++index;
                }
            }
        }
    }
}

TEST_CASE("Grid4D out-of-bounds access")
{
    auto grid = elec::grid::Grid4D {2, 2, 2, 2};

    SECTION("Calling get_at out-of-bounds should throw")
    {
        REQUIRE_THROWS_AS(grid.get_at(2, 2, 2, 2), std::runtime_error);
        REQUIRE_THROWS_AS(grid.get_at(3, 3, 3, 3), std::runtime_error);
    }

    SECTION("Calling set_at out-of-bounds should throw")
    {
        REQUIRE_THROWS_AS(grid.set_at(2, 2, 2, 2, 6.0), std::runtime_error);
        REQUIRE_THROWS_AS(grid.set_at(3, 3, 3, 3, 6.0), std::runtime_error);
    }
}
