#pragma once

#include <string>

namespace elec
{

enum class AtomLabel
{
    H,
    He,
    Li,
    Be,
    B,
    C,
    N,
    O,
    F
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
