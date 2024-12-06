#include <sstream>
#include <stdexcept>
#include <string_view>

#include "extern/mapbox/eternal.hpp"
#include "extern/tomlplusplus/toml.hpp"

#include "elecstruct/input_file_parser/input_file_options.hpp"

namespace
{

using estr = mapbox::eternal::string;
using IFG = elec::InitialFockGuess;

constexpr auto map_string_to_initial_fock_guess = mapbox::eternal::map<estr, IFG>({
    {"zero",             IFG::ZERO_MATRIX            },
    {"extended_huckel",  IFG::EXTENDED_HUCKEL_MATRIX },
    {"core_hamiltonian", IFG::CORE_HAMILTONIAN_MATRIX}
});

auto parse_initial_fock_guess(const toml::table& table) -> IFG
{
    const auto ifg_string = table["initial_fock_guess"].as_string();
    if (!ifg_string) {
        throw std::runtime_error {"Failed to parse 'initial_fock_guess'\n"};
    }

    const auto str = *ifg_string->value_exact<std::string>();

    if (map_string_to_initial_fock_guess.find(str.c_str()) == map_string_to_initial_fock_guess.end()) {
        auto err_msg = std::stringstream {};
        err_msg << "ERROR: Invalid choice of initial fock guess.\n";
        err_msg << "Allowed options: \n";
        for (const auto& item : map_string_to_initial_fock_guess) {
            err_msg << "  - " << item.first.c_str() << '\n';
        }
        err_msg << '\n';

        throw std::runtime_error {err_msg.str()};
    }

    return map_string_to_initial_fock_guess.at(str.c_str());
}

}  // anonymous namespace
