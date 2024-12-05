#include <sstream>
#include <stdexcept>
#include <string_view>

#include "extern/tomlplusplus/toml.hpp"
#include "extern/mapbox/eternal.hpp"

#include "elecstruct/input_file_parser/input_file_options.hpp"

namespace
{

using estr = mapbox::eternal::string;

constexpr auto map_string_to_verbose = mapbox::eternal::map<estr, elec::Verbose>({
    {"true", elec::Verbose::TRUE},
    {"false", elec::Verbose::FALSE}
});

auto parse_verbose(const toml::table& table) -> elec::Verbose
{
    const auto verbose_string = table["verbose"].as_string();
    if (!verbose_string) {
        throw std::runtime_error {"Failed to parse 'verbose'\n"};
    }

    const auto str = *verbose_string->value_exact<std::string>();

    if (map_string_to_verbose.find(str.c_str()) == map_string_to_verbose.end()) {
        auto err_msg = std::stringstream {};
        err_msg << "ERROR: Invalid choice of 'verbose'.\n";
        err_msg << "Allowed options: \n";
        for (const auto& item : map_string_to_verbose) {
            err_msg << "  - " << item.first.c_str() << '\n';
        }
        err_msg << '\n';

        throw std::runtime_error {err_msg.str()};
    }

    return map_string_to_verbose.at(str.c_str());
}

}  // anonymous namespace
