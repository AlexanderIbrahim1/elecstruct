#include <cstdint>

#include "extern/tomlplusplus/toml.hpp"

namespace
{

auto parse_max_hartree_fock_iterations(const toml::table& table) -> std::size_t
{
    const auto n_max_toml = table["max_hartree_fock_iterations"].as_integer();
    if (!n_max_toml) {
        throw std::runtime_error {"Failed to parse 'max_hartree_fock_iterations'\n"};
    }

    const auto n_max = *n_max_toml->value_exact<std::int64_t>();

    if (n_max <= 0) {
        throw std::runtime_error {"'max_hartree_fock_iterations' must be positive\n"};
    }

    return static_cast<std::size_t>(n_max);
}

}  // anonymous namespace
