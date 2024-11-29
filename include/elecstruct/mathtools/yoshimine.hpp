#pragma once

#include <cmath>
#include <cstdint>

namespace elec
{

/*
    Simplify the four indices in the two-electron integral into a single index.
*/
inline auto yoshimine_sort(std::int64_t a, std::int64_t b, std::int64_t c, std::int64_t d) -> std::int64_t
{
    const auto ab = [&]() {
        if (a > b) {
            return a * (a + 1) / 2 + b;
        } else {
            return b * (b + 1) / 2 + a;
        }
    }();

    const auto cd = [&]() {
        if (c > d) {
            return c * (c + 1) / 2 + d;
        } else {
            return d * (d + 1) / 2 + c;
        }
    }();

    const auto abcd = [&]() {
        if (ab > cd) {
            return ab * (ab + 1) / 2 + cd;
        } else {
            return cd * (cd + 1) / 2 + ab;
        }
    }();

    return abcd;
}

}  // namespace elec
