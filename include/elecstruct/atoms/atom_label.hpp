#pragma once

#include <string>

#include <extern/mapbox/eternal.hpp>

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

constexpr auto atom_charge_map = mapbox::eternal::map<AtomLabel, double>({
    {AtomLabel::H,  1.0},
    {AtomLabel::He, 2.0},
    {AtomLabel::Li, 3.0},
    {AtomLabel::Be, 4.0},
    {AtomLabel::B,  5.0},
    {AtomLabel::C,  6.0},
    {AtomLabel::N,  7.0},
    {AtomLabel::O,  8.0},
    {AtomLabel::F,  9.0}
});

auto atom_label_from_name(const std::string&) -> AtomLabel;

auto atom_name_from_label(AtomLabel) -> std::string;

}  // namespace elec
