#include <tuple>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "elecstruct/integrals/two_electron_integral_grid.hpp"

using Indices = std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>;

constexpr auto swap_0_and_1(const Indices& indices) noexcept -> Indices
{
    return {std::get<1>(indices), std::get<0>(indices), std::get<2>(indices), std::get<3>(indices)};
}

constexpr auto swap_2_and_3(const Indices& indices) noexcept -> Indices
{
    return {std::get<0>(indices), std::get<1>(indices), std::get<3>(indices), std::get<2>(indices)};
}

constexpr auto swap_01_and_23(const Indices& indices) noexcept -> Indices
{
    return {std::get<2>(indices), std::get<3>(indices), std::get<0>(indices), std::get<1>(indices)};
}

enum class SwapInstruction
{
    INDICES_01,
    INDICES_23,
    INDICES_01_23
};

constexpr auto swap_indices(const std::vector<SwapInstruction>& swaps, const Indices& original) noexcept -> Indices
{
    using SI = SwapInstruction;

    auto indices = original;

    for (auto instr : swaps) {
        if (instr == SI::INDICES_01) {
            indices = swap_0_and_1(indices);
        }
        else if (instr == SI::INDICES_23) {
            indices = swap_2_and_3(indices);
        }
        else {
            indices = swap_01_and_23(indices);
        }
    }

    return indices;
}

auto yoshimine_helper(const Indices& indices) noexcept -> std::size_t
{
    const auto yoshi = elec::yoshimine_sort;

    return yoshi(std::get<0>(indices), std::get<1>(indices), std::get<2>(indices), std::get<3>(indices));
}

TEST_CASE("yoshimine indices")
{
    /* values are taken from handdone calculations */
    struct TestPair
    {
        Indices input;
        std::size_t expected;
    };

    const auto pair = GENERATE(
        TestPair {
            {1, 1, 1, 1},
            5
    },
        TestPair {{2, 1, 1, 1}, 12},
        TestPair {{2, 2, 1, 1}, 17},
        TestPair {{2, 1, 2, 1}, 14},
        TestPair {{2, 2, 2, 1}, 19},
        TestPair {{2, 2, 2, 2}, 20}
    );

    const auto actual = yoshimine_helper(pair.input);

    REQUIRE(pair.expected == actual);
}

TEST_CASE("yoshimine sort on {0, 1, 2, 3}")
{
    using SI = SwapInstruction;

    const auto original_indices = Indices {0, 1, 2, 3};
    const auto yoshimine_original = yoshimine_helper(original_indices);

    SECTION("swap 01")
    {
        const auto instructions = std::vector<SwapInstruction> {SI::INDICES_01};
        const auto indices = swap_indices(instructions, original_indices);

        REQUIRE(std::get<0>(indices) == 1);
        REQUIRE(std::get<1>(indices) == 0);
        REQUIRE(std::get<2>(indices) == 2);
        REQUIRE(std::get<3>(indices) == 3);

        const auto yoshimine_new = yoshimine_helper(indices);
        REQUIRE(yoshimine_original == yoshimine_new);
    }

    SECTION("swap 23")
    {
        const auto instructions = std::vector<SwapInstruction> {SI::INDICES_23};
        const auto indices = swap_indices(instructions, original_indices);

        REQUIRE(std::get<0>(indices) == 0);
        REQUIRE(std::get<1>(indices) == 1);
        REQUIRE(std::get<2>(indices) == 3);
        REQUIRE(std::get<3>(indices) == 2);

        const auto yoshimine_new = yoshimine_helper(indices);
        REQUIRE(yoshimine_original == yoshimine_new);
    }

    SECTION("swap 01 and 23")
    {
        const auto instructions = std::vector<SwapInstruction> {SI::INDICES_01_23};
        const auto indices = swap_indices(instructions, original_indices);

        REQUIRE(std::get<0>(indices) == 2);
        REQUIRE(std::get<1>(indices) == 3);
        REQUIRE(std::get<2>(indices) == 0);
        REQUIRE(std::get<3>(indices) == 1);

        const auto yoshimine_new = yoshimine_helper(indices);
        REQUIRE(yoshimine_original == yoshimine_new);
    }
}

TEST_CASE("yoshimine sort multiple swaps")
{
    using SI = SwapInstruction;

    const auto original_indices = Indices {0, 1, 2, 3};
    const auto yoshimine_original = yoshimine_helper(original_indices);

    const auto instructions = GENERATE(
        std::vector<SI> {SI::INDICES_01},
        std::vector<SI> {SI::INDICES_01_23, SI::INDICES_23},
        std::vector<SI> {SI::INDICES_01, SI::INDICES_23},
        std::vector<SI> {SI::INDICES_23, SI::INDICES_01},
        std::vector<SI> {SI::INDICES_23, SI::INDICES_01, SI::INDICES_01_23},
        std::vector<SI> {SI::INDICES_01, SI::INDICES_01_23, SI::INDICES_23}
    );

    const auto indices = swap_indices(instructions, original_indices);
    const auto yoshimine_new = yoshimine_helper(indices);
    REQUIRE(yoshimine_original == yoshimine_new);
}

TEST_CASE("two_electron_integral_grid")
{
    SECTION("before and after setting")
    {
        const auto i0 = std::size_t {0};
        const auto i1 = std::size_t {1};
        const auto i2 = std::size_t {2};
        const auto i3 = std::size_t {3};
        const auto value = double {123.456};

        auto grid = elec::TwoElectronIntegralGrid {};

        REQUIRE(!grid.exists(i0, i1, i2, i3));

        grid.set(i0, i1, i2, i3, value);

        REQUIRE(grid.exists(i0, i1, i2, i3));
        REQUIRE_THAT(grid.get(i0, i1, i2, i3), Catch::Matchers::WithinRel(value));
    }

    SECTION("swapping does not change existence")
    {
        const auto original = Indices {4, 3, 5, 1};
        const auto swapped_01 = Indices {3, 4, 5, 1};
        const auto swapped_23 = Indices {4, 3, 1, 5};
        const auto swapped_01_23 = Indices {5, 1, 4, 3};

        auto grid = elec::TwoElectronIntegralGrid {};

        REQUIRE(!grid.exists(std::get<0>(original), std::get<1>(original), std::get<2>(original), std::get<3>(original))
        );
        REQUIRE(!grid.exists(
            std::get<0>(swapped_01), std::get<1>(swapped_01), std::get<2>(swapped_01), std::get<3>(swapped_01)
        ));
        REQUIRE(!grid.exists(
            std::get<0>(swapped_23), std::get<1>(swapped_23), std::get<2>(swapped_23), std::get<3>(swapped_23)
        ));
        REQUIRE(!grid.exists(
            std::get<0>(swapped_01_23),
            std::get<1>(swapped_01_23),
            std::get<2>(swapped_01_23),
            std::get<3>(swapped_01_23)
        ));

        grid.set(std::get<0>(original), std::get<1>(original), std::get<2>(original), std::get<3>(original), 1.0);

        REQUIRE(grid.exists(std::get<0>(original), std::get<1>(original), std::get<2>(original), std::get<3>(original))
        );
        REQUIRE(grid.exists(
            std::get<0>(swapped_01), std::get<1>(swapped_01), std::get<2>(swapped_01), std::get<3>(swapped_01)
        ));
        REQUIRE(grid.exists(
            std::get<0>(swapped_23), std::get<1>(swapped_23), std::get<2>(swapped_23), std::get<3>(swapped_23)
        ));
        REQUIRE(grid.exists(
            std::get<0>(swapped_01_23),
            std::get<1>(swapped_01_23),
            std::get<2>(swapped_01_23),
            std::get<3>(swapped_01_23)
        ));
    }
}
