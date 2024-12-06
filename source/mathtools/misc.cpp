#include <cstdint>

#include "elecstruct/mathtools/misc.hpp"

/*
    This file contains code for various math functions that are used throughout
    the codebase but don't fall into any specific category.
*/

namespace elec::math
{

auto neg_1_power(std::int64_t arg) -> std::int64_t
{
    if (arg % 2 == 0) {
        return 1;
    }
    else {
        return -1;
    }
}

auto factorial(std::int64_t n) -> std::int64_t
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

auto double_factorial(std::int64_t n) -> std::int64_t
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

}  // namespace elec::math
