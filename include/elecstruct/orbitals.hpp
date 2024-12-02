#pragma once

#include <cstdint>
#include <vector>

namespace elec
{

/*
    The labels for the hydrogen atom orbitals.

    Due to the naming rules for the language, as well as naming conventions for enums, we
    label the orbital using a capital letter, and have the letter come first.

    For example, orbital '1s' is expressed as 'S1', and so on.
*/
enum class AtomicOrbitalLabel
{
    S1,
    S2,
    P2,
};

/*
    This object holds the number of angular momentum quanta for a given orbital.

    NOTE: there are many cases in the code where we perform an operation on an angular quantum
    number, with a negative result, which is then used in a function. To avoid constantly calling
    `static_cast()` before using these results, we can make it a signed integer to begin with.
*/
struct AngularMomentumNumbers
{
    std::int64_t x;
    std::int64_t y;
    std::int64_t z;
};

auto atomic_orbitals_to_angular_momentum_numbers(AtomicOrbitalLabel label) -> std::vector<AngularMomentumNumbers>;

/*
    The total angular momentum is required frequently enough that we may as well create
    a convenience function to calculate it.
*/
auto total_angular_momentum(const AngularMomentumNumbers& ang_mom_nums) -> std::int64_t;

}  // namespace elec
