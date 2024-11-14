#pragma once

#include <cstdint>

/*
    This file contains code for various math functions that are used throughout
    the codebase but don't fall into any specific category.
*/

namespace elec::math
{

inline auto neg_1_power(std::int64_t arg) -> std::int64_t
{
    if (arg % 2 == 0) {
        return 1;
    } else {
        return -1;
    }
}

}  // namespace elec::math
