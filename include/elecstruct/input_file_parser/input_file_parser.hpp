#pragma once

#include <any>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>

#include "extern/tomlplusplus/toml.hpp"

#include "elecstruct/atoms.hpp"
#include "elecstruct/input_file_parser/parsed_information.hpp"

/*
    This source file contains the InputFileParser class, which parses a toml file to get
    the input for the electronic structure calculation.
*/

/*
TODO:
    - add initial Fock matrix guess as option
    - add maximum number of iterations as an option
    - add density matrix difference tolerance as an option
    - add energy difference tolerance as an option
*/

namespace elec
{

class InputFileParser
{
public:
    explicit InputFileParser(std::istream& toml_stream);

    auto error_message() const noexcept -> std::string;
    auto is_valid() const noexcept -> bool;

    auto parse(InputFileKey key) -> void;
    auto parse_all() -> void;

    auto parsed_information() const -> const ParsedInformation&;

private:
    bool parse_success_flag_ {};
    std::string error_message_ {};
    ParsedInformation parsed_information_ {};
    std::optional<toml::table> table_;

    auto parse_stream_to_table_(std::istream& toml_stream) noexcept -> std::optional<toml::table>;
};

}  // namespace elec
