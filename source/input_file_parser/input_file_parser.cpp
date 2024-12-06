/*
This source file contains the InputFileParser class, which parses a toml file to get
the input for the electronic structure calculation.
*/

#include <filesystem>
#include <iosfwd>
#include <stdexcept>
#include <string>
#include <vector>

#include "extern/tomlplusplus/toml.hpp"

#include "elecstruct/atoms.hpp"
#include "elecstruct/input_file_parser/input_file_parser.hpp"

#include "parse_atom_information.cpp"
#include "parse_initial_fock_guess.cpp"
#include "parse_max_hartree_fock_iterations.cpp"
#include "parse_tol_change_density_matrix.cpp"
#include "parse_tol_change_hartree_fock_energy.cpp"
#include "parse_n_electrons.cpp"
#include "parse_verbose.cpp"

namespace elec
{

InputFileParser::InputFileParser(std::istream& toml_stream)
    : table_ {parse_stream_to_table_(toml_stream)}
{
    parse_stream_to_table_(toml_stream);
}

auto InputFileParser::is_valid() const noexcept -> bool
{
    return parse_success_flag_;
}

auto InputFileParser::error_message() const noexcept -> std::string
{
    return error_message_;
}

auto InputFileParser::parse_stream_to_table_(std::istream& toml_stream) noexcept -> std::optional<toml::table>
{
    try {
        const auto table = toml::parse(toml_stream);
        parse_success_flag_ = true;

        return std::optional<toml::table> {table};
    }
    catch (const toml::parse_error& err) {
        parse_success_flag_ = false;
        error_message_ = err.what();

        return std::nullopt;
    }
    catch (const std::runtime_error& err) {
        parse_success_flag_ = false;
        error_message_ = err.what();

        return std::nullopt;
    }
}

auto InputFileParser::parse(InputFileKey key) -> void
{
    using IFK = InputFileKey;

    if (!table_) {
        throw std::runtime_error {"ERROR: cannot parse data from table, if table was not parsed from toml file."};
    }

    const auto& table = table_.value();

    switch (key)
    {
        case IFK::ATOM_INFORMATION: {
            parsed_information_[key] = parse_atoms(table);
            break;
        }
        case IFK::INITIAL_FOCK_GUESS: {
            parsed_information_[key] = parse_initial_fock_guess(table);
            break;
        }
        case IFK::MAX_HARTREE_FOCK_ITERATIONS: {
            parsed_information_[key] = parse_max_hartree_fock_iterations(table);
            break;
        }
        case IFK::TOL_CHANGE_DENSITY_MATRIX: {
            parsed_information_[key] = parse_tol_change_density_matrix(table);
            break;
        }
        case IFK::TOL_CHANGE_HARTREE_FOCK_ENERGY: {
            parsed_information_[key] = parse_tol_change_hartree_fock_energy(table);
            break;
        }
        case IFK::N_ELECTRONS: {
            parsed_information_[key] = parse_n_electrons(table);
            break;
        }
        case IFK::VERBOSE: {
            parsed_information_[key] = parse_verbose(table);
            break;
        }
        default: {
            // unreachable; maybe upgrade to C++23 to get the unreachable attribute
        }
    }
}

auto InputFileParser::parse_all() -> void
{
    using IFK = InputFileKey;

    parse(IFK::ATOM_INFORMATION);
    parse(IFK::INITIAL_FOCK_GUESS);
    parse(IFK::MAX_HARTREE_FOCK_ITERATIONS);
    parse(IFK::TOL_CHANGE_DENSITY_MATRIX);
    parse(IFK::TOL_CHANGE_HARTREE_FOCK_ENERGY);
    parse(IFK::N_ELECTRONS);
    parse(IFK::VERBOSE);
}

auto InputFileParser::parsed_information() const -> const ParsedInformation&
{
    return parsed_information_;
}

}  // namespace elec
