#pragma once

#include "extern/mapbox/eternal.hpp"


#include "elecstruct/atoms.hpp"


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

constexpr auto zeta1_exponents_sto3g = mapbox::eternal::map<AtomLabel, double> ({
    {AtomLabel::H, 1.24},
    {AtomLabel::He, 2.0925},
    {AtomLabel::Li, 2.69},
    {AtomLabel::Be, 3.68},
    {AtomLabel::B, 4.68},
    {AtomLabel::C, 5.67},
    {AtomLabel::N, 6.67},
    {AtomLabel::O, 7.66},
    {AtomLabel::F, 8.65}
});

constexpr auto zeta2_exponents_sto3g = mapbox::eternal::map<AtomLabel, double> ({
    {AtomLabel::Li, 0.75},
    {AtomLabel::Be, 1.10},
    {AtomLabel::B,  1.45},
    {AtomLabel::C,  1.72},
    {AtomLabel::N,  1.95},
    {AtomLabel::O,  2.25},
    {AtomLabel::F,  2.55}
});

/*
    In the STO-NG expansions, the Slater type orbital is expressed as a linear combination of
    Gaussians, and each Gaussian is expressed using two values: the coefficient and exponent.

    The GaussianExpansionPairInfo struct holds this information.
*/
struct GaussianExpansionPairInfo
{
    double coefficient;
    double exponent;
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

struct GaussianConstantsSTO3G
{
    double coeff0;
    double coeff1;
    double coeff2;
    double expon0;
    double expon1;
    double expon2;
};

// clang-format off
constexpr auto gaussian_constants_sto3g = mapbox::eternal::map<AtomicOrbitalLabel, GaussianConstantsSTO3G> ({
    {AtomicOrbitalLabel::S1, { 0.4446345422e+00,  0.5353281423e+00,  0.1543289673e+00,  0.109818e+00,  0.405771e+00,  0.222766e+01}},
    {AtomicOrbitalLabel::S2, { 0.7001154689e+00,  0.3995128261e+00, -0.9996722919e-01,  0.751386e-01,  0.231031e+00,  0.994203e+00}},
    {AtomicOrbitalLabel::P2, { 0.3919573931E+00,  0.6076837186E+00,  0.1559162750E+00,  0.751386e-01,  0.231031e+00,  0.994203e+00}}
});
// clang-format on

}  // namespace elec
