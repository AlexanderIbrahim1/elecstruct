#include <algorithm>
#include <cmath>
#include <cstdint>
#include <tuple>

#include "elecstruct/mathtools/n_choose_k.hpp"

#include "elecstruct/integrals/f_coefficient.hpp"

namespace
{

auto range_limits(std::int64_t angmom_j, std::int64_t angmom_l, std::int64_t angmom_m)
    -> std::tuple<std::int64_t, std::int64_t>
{
    const auto minimum = [&]()
    {
        if (angmom_j > angmom_m) {
            return angmom_j - angmom_m;
        }
        else {
            return std::int64_t {0};
        }
    }();

    const auto maximum = std::min(angmom_j, angmom_l) + 1;

    return {minimum, maximum};
}

}  // anonymous namespace

namespace elec
{

auto f_coefficient(
    std::int64_t angmom_j,
    std::int64_t angmom_l,
    std::int64_t angmom_m,
    double separation0,
    double separation1
) -> double
{
    auto result = double {0.0};

    const auto [minimum, maximum] = range_limits(angmom_j, angmom_l, angmom_m);
    for (std::int64_t angmom_k {minimum}; angmom_k < maximum; ++angmom_k) {
        const auto binom0 = elec::math::N_CHOOSE_K_GRID.at(angmom_l, angmom_k);
        const auto coeff0 = std::pow(separation0, angmom_l - angmom_k);

        const auto binom1 = elec::math::N_CHOOSE_K_GRID.at(angmom_m, angmom_j - angmom_k);
        const auto coeff1 = std::pow(separation1, angmom_m - angmom_j + angmom_k);

        const auto contribution = static_cast<double>(binom0 * binom1) * coeff0 * coeff1;

        result += contribution;
    }

    return result;
}

}  // namespace elec
