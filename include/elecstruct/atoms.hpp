#pragma once

#include <string>

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
    AtomLabel label;
    double x;
    double y;
    double z;
};

auto atom_label_from_name(const std::string&) -> AtomLabel;

auto atom_name_from_label(AtomLabel) -> std::string;

}  // namespace elec
