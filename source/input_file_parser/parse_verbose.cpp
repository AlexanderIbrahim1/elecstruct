#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string_view>

#include "extern/mapbox/eternal.hpp"
#include "extern/tomlplusplus/toml.hpp"

#include "elecstruct/input_file_parser/input_file_options.hpp"

namespace
{

constexpr auto map_bool_to_verbose = mapbox::eternal::map<bool, elec::Verbose>({
    {true,  elec::Verbose::TRUE },
    {false, elec::Verbose::FALSE}
});

auto parse_verbose(const toml::table& table) -> elec::Verbose
{
    const auto verbose_toml = table["verbose"].as_boolean();
    if (!verbose_toml) {
        throw std::runtime_error {"Failed to parse 'verbose'\n"};
    }

    const auto verbose = *verbose_toml->value_exact<bool>();

    if (map_bool_to_verbose.find(verbose) == map_bool_to_verbose.end()) {
        auto err_msg = std::stringstream {};
        err_msg << "ERROR: Invalid choice of 'verbose'.\n";
        err_msg << "Allowed options: \n";
        for (const auto& item : map_bool_to_verbose) {
            err_msg << "  - " << std::boolalpha << item.first << '\n';
        }
        err_msg << '\n';

        throw std::runtime_error {err_msg.str()};
    }

    return map_bool_to_verbose.at(verbose);
}

}  // anonymous namespace
