#pragma once

#include <cstdint>
#include "elecstruct/mathtools/compile_time_grid2d.hpp"

namespace elec::math
{

/*
    Angular momenta don't get very large in electronic structure theory, and having a grid
    to calculate N-choose-K instead of an equation is both faster and probably simpler to
    understand than having a recursive equation.
*/

// clang-format off
constexpr auto N_CHOOSE_K_GRID = elec::grid::CompileTimeGrid2D<std::int64_t, 11, 11> {
    {
    1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    1,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,
    1,   3,   3,   1,   0,   0,   0,   0,   0,   0,   0,
    1,   4,   6,   4,   1,   0,   0,   0,   0,   0,   0,
    1,   5,  10,  10,   5,   1,   0,   0,   0,   0,   0,
    1,   6,  15,  20,  15,   6,   1,   0,   0,   0,   0,
    1,   7,  21,  35,  35,  21,   7,   1,   0,   0,   0,
    1,   8,  28,  56,  70,  56,  28,   8,   1,   0,   0,
    1,   9,  36,  84, 126, 126,  84,  36,   9,   1,   0,
    1,  10,  45, 120, 210, 252, 210, 120,  45,  10,   1
    },
};
// clang-format on

}  // namespace elec::math
