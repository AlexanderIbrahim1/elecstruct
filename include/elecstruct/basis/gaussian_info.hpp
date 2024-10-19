#pragma once

namespace elec
{

/*
    In the STO-NG expansions, the Slater type orbital is expressed as a linear combination of
    Gaussians, and each Gaussian is expressed using two values: the coefficient and exponent.

    The GaussianExpansionPairInfo struct holds this information.
*/
struct GaussianInfo
{
    double coefficient;
    double exponent;
};

}  // namespace elec
