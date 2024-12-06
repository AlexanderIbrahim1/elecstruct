#include "extern/tomlplusplus/toml.hpp"

namespace
{

[[maybe_unused]]
auto parse_tol_change_hartree_fock_energy(const toml::table& table) -> double
{
    const auto tol_change_toml = table["tol_change_hartree_fock_energy"].as_floating_point();
    if (!tol_change_toml) {
        throw std::runtime_error {"Failed to parse 'tol_change_hartree_fock_energy'\n"};
    }

    const auto tol_change = *tol_change_toml->value_exact<double>();

    if (tol_change <= 0.0) {
        throw std::runtime_error {"'tol_change_hartree_fock_energy' must be positive\n"};
    }

    return tol_change;
}

}  // anonymous namespace
