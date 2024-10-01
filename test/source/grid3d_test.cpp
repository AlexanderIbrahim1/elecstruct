#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "elecstruct/grids/grid3d.hpp"


TEST_CASE("Grid3D construction with sizes") {
    SECTION("Create a 3x3x3 Grid3D and verify sizes") {
        const auto grid = grid::Grid3D{3, 3, 3};
        const auto [size0, size1, size2] = grid.sizes();
        REQUIRE(size0 == 3);
        REQUIRE(size1 == 3);
        REQUIRE(size2 == 3);
    }

    SECTION("Create a 3x4x5 Grid3D and verify sizes") {
        const auto grid = grid::Grid3D{3, 4, 5};
        const auto [size0, size1, size2] = grid.sizes();
        REQUIRE(size0 == 3);
        REQUIRE(size1 == 4);
        REQUIRE(size2 == 5);
    }
}

TEST_CASE("Grid3D set and get in bounds") {
    auto grid = grid::Grid3D{3, 3, 3};

    SECTION("Set and get values in bounds using set() and get()") {
        grid.set(0, 0, 0, 1.0);
        grid.set(2, 2, 2, 3.0);
        REQUIRE_THAT(grid.get(0, 0, 0), Catch::Matchers::WithinRel(1.0));
        REQUIRE_THAT(grid.get(2, 2, 2), Catch::Matchers::WithinRel(3.0));
    }

    SECTION("Set and get values in bounds using set_at() and get_at()") {
        grid.set_at(1, 1, 1, 2.0);
        grid.set_at(2, 1, 0, 4.0);
        REQUIRE_THAT(grid.get_at(1, 1, 1), Catch::Matchers::WithinRel(2.0));
        REQUIRE_THAT(grid.get_at(2, 1, 0), Catch::Matchers::WithinRel(4.0));
    }
}

TEST_CASE("Grid3D constructed from std::vector<double>") {
    // clang-format off
    const auto data = std::vector<double> {
        0.0,  1.0,  2.0,  3.0,
        4.0,  5.0,  6.0,  7.0,
        8.0,  9.0, 10.0, 11.0,

        12.0, 13.0, 14.0, 15.0,
        16.0, 17.0, 18.0, 19.0,
        20.0, 21.0, 22.0, 23.0
    };
    // clang-format on

    const auto grid = grid::Grid3D{data, 2, 3, 4};

    std::size_t index = 0;
    for (std::size_t i0 = 0; i0 < 2; ++i0) {
        for (std::size_t i1 = 0; i1 < 3; ++i1) {
            for (std::size_t i2 = 0; i2 < 4; ++i2) {
                REQUIRE_THAT(grid.get(i0, i1, i2), Catch::Matchers::WithinRel(data[index]));
                ++index;
            }
        }
    }
}


TEST_CASE("Grid3D out-of-bounds access") {
    auto grid = grid::Grid3D{3, 3, 3};

    SECTION("Calling get_at out-of-bounds should throw") {
        REQUIRE_THROWS_AS(grid.get_at(3, 3, 3), std::runtime_error);
        REQUIRE_THROWS_AS(grid.get_at(3, 4, 2), std::runtime_error);
    }

    SECTION("Calling set_at out-of-bounds should throw") {
        REQUIRE_THROWS_AS(grid.set_at(3, 3, 3, 6.0), std::runtime_error);
        REQUIRE_THROWS_AS(grid.set_at(2, 3, 4, 6.0), std::runtime_error);
    }
}
