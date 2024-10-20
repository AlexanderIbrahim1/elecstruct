#pragma once

#include <array>
#include <cstdint>

#include "elecstruct/integrals/integrals.hpp"

namespace impl_elec
{

template <typename T, std::size_t SIZE0, std::size_t SIZE1>
class CompileTimeSquareGrid2D
{
public:
    using Data = std::array<T, SIZE0 * SIZE1>;

    constexpr CompileTimeSquareGrid2D(const Data& data)
        : data_ {data}
    {}

    constexpr auto at(std::size_t i0, std::size_t i1) const noexcept -> const T&
    {
        return data_[index(i0, i1)];
    }

    constexpr auto at(std::size_t i0, std::size_t i1) noexcept -> T&
    {
        return data_[index(i0, i1)];
    }

    constexpr auto size0() const noexcept -> std::size_t
    {
        return SIZE0;
    }

    constexpr auto size1() const noexcept -> std::size_t
    {
        return SIZE1;
    }

private:
    constexpr auto index(std::size_t i0, std::size_t i1) const noexcept -> std::size_t
    {
        return i1 + SIZE0 * i0;
    }

    Data data_;
};

/*
    Angular momenta don't get very large in electronic structure theory, and having a grid
    to calculate N-choose-K instead of an equation is both faster and probably simpler to
    understand than having a recursive equation.
*/

// clang-format off
constexpr auto N_CHOOSE_K_GRID = CompileTimeSquareGrid2D<std::uint64_t, 11, 11> {
    {
    1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    1,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,
    1,   3,   3,   1,   0,   0,   0,   0,   0,   0,   0,
    1,   4,   6,   4,   1,   0,   0,   0,   0,   0,   0,
    1,   5,  10,  10,   5,   1,   0,   0,   0,   0,   0,
    1,   6,  15,  20,  15,   6,   1,   0,   0,   0,   0,
    1,   7,  21,  35,  35,  21,   7,   1,   0,   0,   0,
    1,   8,  28,  56,  70,  56,  28,   8,   1,   0,   0,
    1,   9,  36,  84, 126, 126,  84,  36,   9,   1,   0,
    1,  10,  45, 120, 210, 252, 210, 120,  45,  10,   1
    },
};
// clang-format on

}  // namespace impl_elec

namespace elec
{

struct OverlapIntegralInfo1D
{
    std::uint64_t angular_momentum;
    double gaussian_exponent;
    double centre_coordinate;
};

inline auto overlap_integral_1d(
    const OverlapIntegralInfo1D& gaussian0,
    const OverlapIntegralInfo1D& gaussian1,
    double total_centre
) -> double
{
    if (gaussian0.angular_momentum >= impl_elec::N_CHOOSE_K_GRID.size0()) {
        throw std::runtime_error {"Encountered an angular momentum beyond the N-choose-K grid size bounds."};
    }

    if (gaussian1.angular_momentum >= impl_elec::N_CHOOSE_K_GRID.size1()) {
        throw std::runtime_error {"Encountered an angular momentum beyond the N-choose-K grid size bounds."};
    }

    auto overlap = double {0.0};

    for (std::uint64_t sub_ang_mom0 {0}; sub_ang_mom0 < gaussian0.angular_momentum + 1; ++sub_ang_mom0) {
        for (std::uint64_t sub_ang_mom1 {0}; sub_ang_mom1 < gaussian1.angular_momentum + 1; ++sub_ang_mom1) {
            if ((sub_ang_mom0 + sub_ang_mom1) % 2 != 0) {
                continue;
            }

            
        }
    }
}

}  // namespace elec
