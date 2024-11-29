#pragma once

#include <cmath>
#include <cstdint>
#include <unordered_map>

namespace elec
{

namespace impl_two_electron
{

/*
    Simplify the four indices in the two-electron integral into a single index.
*/
inline auto yoshimine_sort(std::size_t a, std::size_t b, std::size_t c, std::size_t d) -> std::size_t
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

}  // namespace elec::impl_two_electron


class TwoElectronIntegralGrid
{
public:

private:
    std::unordered_map<std::size_t, double> grid_;
};

}  // namespace elec

