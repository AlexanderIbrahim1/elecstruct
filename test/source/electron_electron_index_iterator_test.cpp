#include <cstdint>
#include <tuple>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include "elecstruct/integrals/electron_electron_index_iterator.hpp"


using Indices = std::tuple<std::int64_t, std::int64_t, std::int64_t, std::int64_t, std::int64_t>;


auto indices_via_standard_nested_loops(
    std::int64_t angmom_0,
    std::int64_t angmom_1,
    std::int64_t angmom_2,
    std::int64_t angmom_3
) -> std::vector<Indices>
{
    auto output = std::vector<Indices> {};

    for (std::int64_t idx_l_01 {0}; idx_l_01 < angmom_0 + angmom_1 + 1; ++idx_l_01)
    for (std::int64_t idx_r_01 {0}; idx_r_01 < static_cast<std::int64_t>(idx_l_01 / 2) + 1; ++idx_r_01)
    for (std::int64_t idx_l_23 {0}; idx_l_23 < angmom_2 + angmom_3 + 1; ++idx_l_23)
    for (std::int64_t idx_r_23 {0}; idx_r_23 < static_cast<std::int64_t>(idx_l_23 / 2) + 1; ++idx_r_23)
    for (std::int64_t idx_i {0}; idx_i < static_cast<std::int64_t>((idx_l_01 + idx_l_23 - 2 * (idx_r_01 + idx_r_23)) / 2) + 1; ++idx_i)
    {
        output.emplace_back(idx_l_01, idx_r_01, idx_l_23, idx_r_23, idx_i);
    }

    return output;
}


auto indices_via_custom_iterator(
    std::int64_t angmom_0,
    std::int64_t angmom_1,
    std::int64_t angmom_2,
    std::int64_t angmom_3
) -> std::vector<Indices>
{
    auto output = std::vector<Indices> {};

    const auto index_generator = elec::ElectronElectronIndexGenerator(angmom_0, angmom_1, angmom_2, angmom_3);

    for (const auto [idx_l_01, idx_r_01, idx_l_23, idx_r_23, idx_i] : index_generator) {
        output.emplace_back(idx_l_01, idx_r_01, idx_l_23, idx_r_23, idx_i);
    }

    return output;
}


TEST_CASE("electron-electron integral index iterator")
{
    struct AngularMomentumNumberPairs
    {
        std::int64_t angmom_0;
        std::int64_t angmom_1;
        std::int64_t angmom_2;
        std::int64_t angmom_3;
    };

    auto pair = GENERATE(
        AngularMomentumNumberPairs {0, 0, 0, 0},
        AngularMomentumNumberPairs {1, 0, 0, 0},
        AngularMomentumNumberPairs {0, 1, 0, 0},
        AngularMomentumNumberPairs {0, 0, 1, 0},
        AngularMomentumNumberPairs {0, 0, 0, 1},
        AngularMomentumNumberPairs {1, 1, 0, 0},
        AngularMomentumNumberPairs {1, 0, 1, 0},
        AngularMomentumNumberPairs {1, 0, 0, 1},
        AngularMomentumNumberPairs {0, 1, 1, 0},
        AngularMomentumNumberPairs {0, 1, 0, 1},
        AngularMomentumNumberPairs {0, 0, 1, 1},
        AngularMomentumNumberPairs {1, 2, 0, 0},
        AngularMomentumNumberPairs {2, 3, 0, 0}
    );

    const auto via_loop = indices_via_standard_nested_loops(pair.angmom_0, pair.angmom_1, pair.angmom_2, pair.angmom_3);
    const auto via_iter = indices_via_custom_iterator(pair.angmom_0, pair.angmom_1, pair.angmom_2, pair.angmom_3);

    REQUIRE_THAT(via_loop, Catch::Matchers::Equals(via_iter));
}
