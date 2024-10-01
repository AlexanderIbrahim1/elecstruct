#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "elecstruct/grids/grid2d.hpp"

TEST_CASE("Grid2D construction with sizes") {
    SECTION("Create a 3x3 Grid2D and verify sizes") {
        const auto grid = grid::Grid2D{3, 3};
        const auto [size0, size1] = grid.sizes();
        REQUIRE(size0 == 3);
        REQUIRE(size1 == 3);
    }

    SECTION("Create a 3x5 Grid2D and verify sizes") {
        const auto grid = grid::Grid2D{3, 5};
        const auto [size0, size1] = grid.sizes();
        REQUIRE(size0 == 3);
        REQUIRE(size1 == 5);
    }
}

TEST_CASE("Grid2D set and get in bounds") {
    auto grid = grid::Grid2D{3, 3};

    SECTION("Set and get values in bounds using set() and get()") {
        grid.set(0, 0, 1.0);
        grid.set(2, 2, 3.0);
        REQUIRE_THAT(grid.get(0, 0), Catch::Matchers::WithinRel(1.0));
        REQUIRE_THAT(grid.get(2, 2), Catch::Matchers::WithinRel(3.0));
    }

    SECTION("Set and get values in bounds using set_at() and get_at()") {
        grid.set_at(1, 1, 2.0);
        grid.set_at(2, 1, 4.0);
        REQUIRE_THAT(grid.get(1, 1), Catch::Matchers::WithinRel(2.0));
        REQUIRE_THAT(grid.get(2, 1), Catch::Matchers::WithinRel(4.0));
    }
}


TEST_CASE("Grid2D constructed from std::vector<double>") {
    // clang-format off
    const auto data = std::vector<double> {
        0.0,  1.0,  2.0,  3.0,  4.0,
        5.0,  6.0,  7.0,  8.0,  9.0,
        10.0, 11.0, 12.0, 13.0, 14.0
    };
    // clang-format on
    
    const auto grid = grid::Grid2D{data, 3, 5};
    
    std::size_t index = 0;
    for (std::size_t i0 = 0; i0 < 3; ++i0) {
        for (std::size_t i1 = 0; i1 < 5; ++i1) {
            REQUIRE_THAT(grid.get(i0, i1), Catch::Matchers::WithinRel(data[index]));
            ++index;
        }
    }
}


TEST_CASE("Grid2D out-of-bounds access") {
    auto grid = grid::Grid2D{3, 3};

    SECTION("Calling get_at out-of-bounds should throw") {
        REQUIRE_THROWS_AS(grid.get_at(3, 3), std::runtime_error);
        REQUIRE_THROWS_AS(grid.get_at(3, 5), std::runtime_error);
    }

    SECTION("Calling set_at out-of-bounds should throw") {
        REQUIRE_THROWS_AS(grid.set_at(3, 3, 6.0), std::runtime_error);
        REQUIRE_THROWS_AS(grid.set_at(3, 5, 6.0), std::runtime_error);
    }
}
