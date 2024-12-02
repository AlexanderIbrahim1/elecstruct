#include <cstdint>
#include <stdexcept>
#include <vector>

#include "elecstruct/orbitals.hpp"

namespace elec
{

auto atomic_orbitals_to_angular_momentum_numbers(AtomicOrbitalLabel label) -> std::vector<AngularMomentumNumbers>
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

/*
    The total angular momentum is required frequently enough that we may as well create
    a convenience function to calculate it.
*/
auto total_angular_momentum(const AngularMomentumNumbers& ang_mom_nums) -> std::int64_t
{
    return ang_mom_nums.x + ang_mom_nums.y + ang_mom_nums.z;
}

}  // namespace elec
