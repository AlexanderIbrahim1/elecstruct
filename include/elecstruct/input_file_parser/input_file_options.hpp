#pragma once

/*
    This header file contains the enum classes for various options that the InputFileParser
    object can accommodate.
*/

namespace elec
{

enum class InitialFockGuess
{
    ZERO_MATRIX,
    EXTENDED_HUCKEL_MATRIX,
    CORE_HAMILTONIAN_MATRIX
};

enum class Verbose
{
    TRUE,
    FALSE
};

}  // namespace elec
