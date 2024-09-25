/*
This source file contains the InputFileParser class, which parses a toml file to get
the input for the electronic structure calculation.
*/

#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <tomlplusplus/toml.hpp>

namespace elec
{

enum class AtomLabel
{
    H,
    C,
    N,
    O
};

struct AtomInfo
{
    std::string label;
    double x;
    double y;
    double z;
};

class InputFileParser
{
public:
    explicit InputFileParser(std::istream& toml_stream);

    constexpr auto error_message() const noexcept -> std::string;
    constexpr auto is_valid() const noexcept -> bool;

    std::vector<AtomInfo> atom_information;

private:
    bool parse_success_flag_ {};
    std::string error_message_ {};

    void parse_helper_(std::istream& toml_stream);
    void parse_atoms_(const toml::table& table);
};

}  // namespace elec
