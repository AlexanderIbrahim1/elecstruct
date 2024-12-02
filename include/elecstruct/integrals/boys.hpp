#pragma once

/*
    A reproduction of the iemplementation of the fast Boys algorithm provided in

    title: "A fast algorithm for computing the Boys function"
    authors: G. Beylkin and S. Sharma
    journal: J. Chem. Phys. 155, 174117 (2021)
    link: https://doi.org/10.1063/5.0062444
*/

namespace elec
{

auto boys_beylkin_sharma(double x, std::size_t n) -> double;

}  // namespace elec
