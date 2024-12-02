/*
This source file contains the InputFileParser class, which parses a toml file to get
the input for the electronic structure calculation.
*/

#include <iostream>

#include <filesystem>
#include <iosfwd>
#include <stdexcept>
#include <string>
#include <vector>

#include "extern/tomlplusplus/toml.hpp"

#include "elecstruct/atoms.hpp"
#include "elecstruct/input_file_parser.hpp"

namespace elec
{

InputFileParser::InputFileParser(std::istream& toml_stream)
{
    parse_helper_(toml_stream);
}

auto InputFileParser::is_valid() const noexcept -> bool
{
    return parse_success_flag_;
}

auto InputFileParser::error_message() const noexcept -> std::string
{
    return error_message_;
}

void InputFileParser::parse_helper_(std::istream& toml_stream)
{
    try {
        const auto table = toml::parse(toml_stream);
        parse_atoms_(table);
        parse_success_flag_ = true;
    }
    catch (const toml::parse_error& err) {
        parse_success_flag_ = false;
        error_message_ = err.what();
    }
    catch (const std::runtime_error& err) {
        parse_success_flag_ = false;
        error_message_ = err.what();
    }
}

void InputFileParser::parse_atoms_(const toml::table& table)
{
    using EmptyOrbitalVector = std::vector<AtomicOrbitalLabel>;

    const auto atom_info = table["positions"].as_array();

    if (!atom_info) {
        throw std::runtime_error("Failed to parse array of atom names and positions.\n");
    }

    if (atom_info->size() == 0) {
        throw std::runtime_error("Found no atom names and positions to parse.\n");
    }

    atom_information.reserve(atom_info->size());

    for (auto&& atom_info_node : *atom_info) {
        const auto atom_and_positions = atom_info_node.as_array();

        if (!atom_and_positions) {
            throw std::runtime_error("Failed to parse row in array of atom names and positions.\n");
        }

        if (atom_and_positions->size() != 4) {
            throw std::runtime_error("Each row needs 4 elements; the atom ID, and the three Cartesian coordinates.");
        }

        const auto atom_name = atom_and_positions->at(0).value<std::string>();
        if (!atom_name) {
            throw std::runtime_error("Failed to parse the atom name.");
        }

        const auto atom_label = atom_label_from_name(atom_name->c_str());

        const auto atom_x_position = atom_and_positions->at(1).value<double>();
        if (!atom_x_position) {
            throw std::runtime_error("Failed to parse the x-position of the atom.");
        }

        const auto atom_y_position = atom_and_positions->at(2).value<double>();
        if (!atom_y_position) {
            throw std::runtime_error("Failed to parse the y-position of the atom.");
        }

        const auto atom_z_position = atom_and_positions->at(3).value<double>();
        if (!atom_z_position) {
            throw std::runtime_error("Failed to parse the z-position of the atom.");
        }

        const auto position = coord::Cartesian3D {*atom_x_position, *atom_y_position, *atom_z_position};

        atom_information.emplace_back(atom_label, position, EmptyOrbitalVector {});
    }
}

}  // namespace elec
