#include <cstdint>
#include <tuple>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include "elecstruct/integrals/nuclear_electron_index_iterator.hpp"


using Indices = std::tuple<std::int64_t, std::int64_t, std::int64_t>;


auto indices_via_standard_nested_loops(std::int64_t angmom_0, std::int64_t angmom_1) -> std::vector<Indices>
{
    auto output = std::vector<Indices> {};

    for (std::int64_t idx_l {0}; idx_l < angmom_0 + angmom_1 + 1; ++idx_l)
    {
        for (std::int64_t idx_r {0}; idx_r < static_cast<std::int64_t>(idx_l / 2) + 1; ++idx_r)
        {
            for (std::int64_t idx_i {0}; idx_i < static_cast<std::int64_t>((idx_l - 2 * idx_r) / 2) + 1; ++idx_i)
            {
                output.push_back({idx_l, idx_r, idx_i});
            }
        }
    }

    return output;
}


auto indices_via_custom_iterator(std::int64_t angmom_0, std::int64_t angmom_1) -> std::vector<Indices>
{
    auto output = std::vector<Indices> {};

    for (const auto [idx_l, idx_r, idx_i] : elec::NuclearElectronIndexGenerator(angmom_0, angmom_1)) {
        output.push_back({idx_l, idx_r, idx_i});
    }

    return output;
}


TEST_CASE("nuclear-electron integral index iterator")
{
    struct AngularMomentumNumberPairs
    {
        std::int64_t angmom_0;
        std::int64_t angmom_1;
    };

    auto pair = GENERATE(
        AngularMomentumNumberPairs {0, 0},
        AngularMomentumNumberPairs {0, 1},
        AngularMomentumNumberPairs {1, 0},
        AngularMomentumNumberPairs {1, 1},
        AngularMomentumNumberPairs {1, 2},
        AngularMomentumNumberPairs {1, 3},
        AngularMomentumNumberPairs {2, 3},
        AngularMomentumNumberPairs {3, 3}
    );

    const auto via_loop = indices_via_standard_nested_loops(pair.angmom_0, pair.angmom_1);
    const auto via_iter = indices_via_custom_iterator(pair.angmom_0, pair.angmom_1);

    REQUIRE_THAT(via_loop, Catch::Matchers::Equals(via_iter));
}
