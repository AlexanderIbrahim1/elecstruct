#pragma once

#include <cstdint>

namespace elec::math
{

// TODO: maybe create a table that this function can wrap around?
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

// TODO: maybe create a table that this function can wrap around?
inline auto factorial(std::int64_t n) -> std::int64_t
{
    if (n == 0) {
        return 1;
    }

    auto result = std::int64_t {1};

    for (std::int64_t value {n}; value >= 2; value -= 1) {
        result *= value;
    }

    return result;
}

}  // namespace elec::math
