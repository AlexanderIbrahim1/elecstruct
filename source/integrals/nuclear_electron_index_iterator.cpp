#include <cstdint>
#include <tuple>

#include "elecstruct/integrals/nuclear_electron_index_iterator.hpp"

/*
    NOTE: Using a custom iterator object takes up a lot of additional code, adds some
    complications, might even be a little bit slower, and these loops don't appear in
    enough places in the code for this encapulation to be worth it.

    However, part of the purpose of writing this code is to improve my general C++
    skills, and so I decided to implement the nested loops in the body of the integral
    for the nuclear-electron attraction using custom iterators. This reduces the amount
    of nesting in the body of the function from 9 to 3.
*/

namespace
{

auto calculate_idx_i_max_(std::int64_t idx_l, std::int64_t idx_r) noexcept -> std::int64_t
{
    return static_cast<std::int64_t>((idx_l - 2 * idx_r) / 2) + 1;
}

auto calculate_idx_r_max_(std::int64_t idx_l) noexcept -> std::int64_t
{
    return static_cast<std::int64_t>(idx_l / 2) + 1;
}

}  // anonymous namespace


namespace elec
{

// --- NuclearElectronIndexIterator

NuclearElectronIndexIterator::NuclearElectronIndexIterator(std::int64_t idx_l_max)
    : idx_l_max_ {idx_l_max}
{}

auto NuclearElectronIndexIterator::operator*() const noexcept -> value_type
{
    return {idx_l_, idx_r_, idx_i_};
}

auto NuclearElectronIndexIterator::operator->() const noexcept -> value_type
{
    return **this;
}

auto NuclearElectronIndexIterator::operator++() noexcept -> NuclearElectronIndexIterator&
{
    ++idx_i_;
    if (idx_i_ == idx_i_max_) {
        idx_i_ = 0;
        ++idx_r_;

        idx_i_max_ = calculate_idx_i_max_(idx_l_, idx_r_);
        if (idx_r_ == idx_r_max_) {
            idx_r_ = 0;
            ++idx_l_;

            idx_r_max_ = calculate_idx_r_max_(idx_l_);
            idx_i_max_ = calculate_idx_i_max_(idx_l_, idx_r_);
        }
    }

    return *this;
}

auto NuclearElectronIndexIterator::operator++(int) noexcept -> NuclearElectronIndexIterator
{
    auto temp = *this;
    ++(*this);
    return temp;
}

void NuclearElectronIndexIterator::set_to_end() noexcept
{
    idx_l_ = idx_l_max_;
    idx_r_ = 0;
    idx_i_ = 0;

    idx_r_max_ = calculate_idx_r_max_(idx_l_);
    idx_i_max_ = calculate_idx_i_max_(idx_l_, idx_r_);
}

bool NuclearElectronIndexIterator::operator==(const NuclearElectronIndexIterator&) const noexcept = default;

auto operator!=(const NuclearElectronIndexIterator& left, const NuclearElectronIndexIterator& right) noexcept -> bool
{
    return !(left == right);
}


// --- NuclearElectronIndexGenerator

NuclearElectronIndexGenerator::NuclearElectronIndexGenerator(std::int64_t angmom_0, std::int64_t angmom_1)
    : idx_l_max_ {angmom_0 + angmom_1 + 1}
{}

auto NuclearElectronIndexGenerator::begin() noexcept -> NuclearElectronIndexIterator
{
    return NuclearElectronIndexIterator {idx_l_max_};
}

auto NuclearElectronIndexGenerator::end() noexcept -> NuclearElectronIndexIterator
{
    auto iterator = NuclearElectronIndexIterator {idx_l_max_};
    iterator.set_to_end();

    return iterator;
}

auto NuclearElectronIndexGenerator::begin() const noexcept -> NuclearElectronIndexIterator
{
    return NuclearElectronIndexIterator {idx_l_max_};
}

auto NuclearElectronIndexGenerator::end() const noexcept -> NuclearElectronIndexIterator
{
    auto iterator = NuclearElectronIndexIterator {idx_l_max_};
    iterator.set_to_end();

    return iterator;
}

}  // namespace elec
