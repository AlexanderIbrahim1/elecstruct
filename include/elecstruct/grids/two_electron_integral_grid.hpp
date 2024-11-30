#pragma once

#include <algorithm>
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
inline auto yoshimine_sort(std::size_t a, std::size_t b, std::size_t c, std::size_t d) noexcept -> std::size_t
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


/*
    The TwoElectronIntegralGrid is a wrapper around a hashmap that converts the four indices
    into a Yoshimine index before interfacing with the hashmap.
*/
class TwoElectronIntegralGrid
{
public:
    auto exists(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3) noexcept -> bool
    {
        const auto i_yoshimine = impl_two_electron::yoshimine_sort(i0, i1, i2, i3);

        return integrals_.find(i_yoshimine) != integrals_.end();
    } 

    auto set(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3, double value) noexcept -> void
    {
        const auto i_yoshimine = impl_two_electron::yoshimine_sort(i0, i1, i2, i3);
        integrals_[i_yoshimine] = value;
    }

    auto get(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3) noexcept -> double
    {
        const auto i_yoshimine = impl_two_electron::yoshimine_sort(i0, i1, i2, i3);

        return integrals_.at(i_yoshimine);
    }

private:
    std::unordered_map<std::size_t, double> integrals_;
};

}  // namespace elec

