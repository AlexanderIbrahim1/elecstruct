#pragma once

#include <cstdint>

namespace elec::math
{

inline auto double_factorial(std::int64_t n) -> std::int64_t
{
    // right now, -1 is the only negative value we account for separately
    if (n == -1) {
        return 1;
    }

    auto result = std::int64_t {1};

    for (std::int64_t value {n}; value >= 2; value -= 2) {
        result *= value;
    }

    return result;
}

inline auto double_factorial_minus_1(std::uint64_t n) -> std::uint64_t
{
    if (n == 0) {
        return 1;
    }

    auto result = std::uint64_t {1};

    // guaranteed that n - 1 >= 0 here
    for (std::uint64_t value {n - 1}; value >= 2; value -= 2) {
        result *= value;
    }

    return result;
}

}  // namespace elec::math
