#pragma once

#include <cstdint>

/*
    This file contains code for various math functions that are used throughout
    the codebase but don't fall into any specific category.
*/

namespace elec::math
{

/*
    Calculate (-1)^(arg)
*/
auto neg_1_power(std::int64_t arg) -> std::int64_t;

/*
    Calculate n!
*/
auto factorial(std::int64_t n) -> std::int64_t;

/*
    Calculate n!!
    NOTE: this function uses the convention that (-1)!! == 1
*/
auto double_factorial(std::int64_t n) -> std::int64_t;

}  // namespace elec::math
