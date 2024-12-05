#include "extern/tomlplusplus/toml.hpp"

namespace
{

auto parse_tol_change_density_matrix(const toml::table& table) -> double
{
    const auto tol_change_toml = table["tol_change_density_matrix"].as_floating_point();
    if (!tol_change_toml) {
        throw std::runtime_error {"Failed to parse 'tol_change_density_matrix'\n"};
    }

    const auto tol_change = *tol_change_toml->value_exact<double>();

    if (tol_change <= 0.0) {
        throw std::runtime_error {"'tol_change_density_matrix' must be positive\n"};
    }

    return tol_change;
}

}  // anonymous namespace
