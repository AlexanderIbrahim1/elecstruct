#include <cstdint>

#include "extern/tomlplusplus/toml.hpp"

namespace
{

auto parse_n_electrons(const toml::table& table) -> std::size_t
{
    const auto n_electrons_toml = table["n_electrons"].as_integer();
    if (!n_electrons_toml) {
        throw std::runtime_error {"Failed to parse 'n_electrons'\n"};
    }

    const auto n_electrons = *n_electrons_toml->value_exact<std::int64_t>();

    if (n_electrons < 0) {
        throw std::runtime_error {"'n_electrons' must be non-negative\n"};
    }

    return static_cast<std::size_t>(n_electrons);
}

}  // anonymous namespace
