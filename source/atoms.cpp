#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>

#include "elecstruct/atoms.hpp"

namespace
{

using AL = elec::AtomLabel;

/*
    I wanted a map between the enum and string versions of the atom names.

    Unfortunately, neither `std::string` nor `std::unordered_map` are constexpr, and so I'm
    using a constexpr global C array that holds entries of pairs of `const char*` and the enum
    class that I want.

    The linear search is not a performance blocker whatsoever for this program, but the
    code I need to deal with this issue is awkward.
*/
struct AtomMapEntry
{
    const char* name;
    elec::AtomLabel label;
};

constexpr AtomMapEntry atom_names_and_labels[] = {
    {"H",  AL::H },
    {"He", AL::He},
    {"Li", AL::Li},
    {"Be", AL::Be},
    {"B",  AL::B },
    {"C",  AL::C },
    {"N",  AL::N },
    {"O",  AL::O },
    {"F",  AL::F }
};

constexpr auto atom_charge_map = mapbox::eternal::map<AL, double>({
    {AL::H,  1.0},
    {AL::He, 2.0},
    {AL::Li, 3.0},
    {AL::Be, 4.0},
    {AL::B,  5.0},
    {AL::C,  6.0},
    {AL::N,  7.0},
    {AL::O,  8.0},
    {AL::F,  9.0}
});

}  // anonymous namespace

namespace elec
{

auto nuclear_charge(AtomLabel label) -> double
{
    return atom_charge_map.at(label);
}

auto atom_label_from_name(const std::string& name) -> AtomLabel
{
    for (const auto& entry : atom_names_and_labels) {
        if (std::strcmp(name.c_str(), entry.name)) {
            return entry.label;
        }
    }

    auto err_msg = std::stringstream {};
    err_msg << "Invalid atom name provided: '" << name << "'";
    throw std::runtime_error {err_msg.str()};
}

auto atom_name_from_label(AtomLabel label) -> std::string
{
    for (const auto& entry : atom_names_and_labels) {
        if (label == entry.label) {
            return std::string {entry.name};
        }
    }

    throw std::runtime_error {"Invalid atom label provided."};
}

}  // namespace elec
