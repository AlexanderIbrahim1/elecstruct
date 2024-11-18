#pragma once

#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "extern/tomlplusplus/toml.hpp"

#include "elecstruct/atoms/atoms.hpp"

/*
This source file contains the InputFileParser class, which parses a toml file to get
the input for the electronic structure calculation.
*/

namespace elec
{

class InputFileParser
{
public:
    explicit InputFileParser(std::istream& toml_stream);

    auto error_message() const noexcept -> std::string;
    auto is_valid() const noexcept -> bool;

    std::vector<AtomInfo> atom_information;

private:
    bool parse_success_flag_ {};
    std::string error_message_ {};

    void parse_helper_(std::istream& toml_stream);
    void parse_atoms_(const toml::table& table);
};

}  // namespace elec
