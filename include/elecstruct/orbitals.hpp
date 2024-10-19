#pragma once

#include <algorithm>
#include <cstddef>
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
*/
struct AngularMomentumNumbers
{
    std::size_t x;
    std::size_t y;
    std::size_t z;
};

inline auto atomic_orbitals_to_angular_momentum_numbers(AtomicOrbitalLabel label) -> std::vector<AngularMomentumNumbers>
{
    using AOL = AtomicOrbitalLabel;

    if (label == AOL::S1 || label == AOL::S2) {
        return {
            {0, 0, 0}
        };
    }
    else if (label == AOL::P2) {
        return {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        };
    }
    else {
        throw std::runtime_error {"UNREACHABLE ERROR: invalid AtomicOrbitalLabel found"};
    }
}

}  // namespace elec
