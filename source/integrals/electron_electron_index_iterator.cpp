#include <cstdint>
#include <tuple>

#include "elecstruct/integrals/electron_electron_index_iterator.hpp"

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

auto calculate_idx_i_max_(
    std::int64_t idx_l_01,
    std::int64_t idx_l_23,
    std::int64_t idx_r_01,
    std::int64_t idx_r_23
) noexcept -> std::int64_t
{
    return static_cast<std::int64_t>((idx_l_01 + idx_l_23 - 2 * (idx_r_01 + idx_r_23)) / 2) + 1;
}

auto calculate_idx_r_max_(std::int64_t idx_l) noexcept -> std::int64_t
{
    return static_cast<std::int64_t>(idx_l / 2) + 1;
}

}  // anonymous namespace

namespace elec
{

// --- ElectronElectronIndexIterator

ElectronElectronIndexIterator::ElectronElectronIndexIterator(std::int64_t idx_l_01_max, std::int64_t idx_l_23_max)
    : idx_l_01_max_ {idx_l_01_max}
    , idx_l_23_max_ {idx_l_23_max}
{}

auto ElectronElectronIndexIterator::operator*() const noexcept -> value_type
{
    return ElectronElectronIntegralIndices {idx_l_01_, idx_r_01_, idx_l_23_, idx_r_23_, idx_i_};
}

auto ElectronElectronIndexIterator::operator->() const noexcept -> value_type
{
    return **this;
}

auto ElectronElectronIndexIterator::operator++() noexcept -> ElectronElectronIndexIterator&
{
    ++idx_i_;
    if (idx_i_ == idx_i_max_) {
        idx_i_ = 0;
        ++idx_r_23_;

        idx_i_max_ = calculate_idx_i_max_(idx_l_01_, idx_l_23_, idx_r_01_, idx_r_23_);

        if (idx_r_23_ == idx_r_23_max_) {
            idx_r_23_ = 0;
            ++idx_l_23_;

            idx_r_23_max_ = calculate_idx_r_max_(idx_l_23_);
            idx_i_max_ = calculate_idx_i_max_(idx_l_01_, idx_l_23_, idx_r_01_, idx_r_23_);

            if (idx_l_23_ == idx_l_23_max_) {
                idx_l_23_ = 0;
                ++idx_r_01_;

                idx_r_23_max_ = calculate_idx_r_max_(idx_l_23_);
                idx_i_max_ = calculate_idx_i_max_(idx_l_01_, idx_l_23_, idx_r_01_, idx_r_23_);

                if (idx_r_01_ == idx_r_01_max_) {
                    idx_r_01_ = 0;
                    ++idx_l_01_;

                    idx_r_01_max_ = calculate_idx_r_max_(idx_l_01_);
                    idx_r_23_max_ = calculate_idx_r_max_(idx_l_23_);
                    idx_i_max_ = calculate_idx_i_max_(idx_l_01_, idx_l_23_, idx_r_01_, idx_r_23_);
                }
            }
        }
    }

    return *this;
}

auto ElectronElectronIndexIterator::operator++(int) noexcept -> ElectronElectronIndexIterator
{
    auto temp = *this;
    ++(*this);
    return temp;
}

void ElectronElectronIndexIterator::set_to_end() noexcept
{
    idx_l_01_ = idx_l_01_max_;
    idx_l_23_ = 0;
    idx_r_01_ = 0;
    idx_r_23_ = 0;
    idx_i_ = 0;

    idx_r_01_max_ = calculate_idx_r_max_(idx_l_01_);
    idx_r_23_max_ = calculate_idx_r_max_(idx_l_23_);
    idx_i_max_ = calculate_idx_i_max_(idx_l_01_, idx_l_23_, idx_r_01_, idx_r_23_);
}

bool ElectronElectronIndexIterator::operator==(const ElectronElectronIndexIterator&) const noexcept = default;

auto operator!=(const ElectronElectronIndexIterator& left, const ElectronElectronIndexIterator& right) noexcept -> bool
{
    return !(left == right);
}

// --- ElectronElectronIndexGenerator

ElectronElectronIndexGenerator::ElectronElectronIndexGenerator(
    std::int64_t angmom_gauss0,
    std::int64_t angmom_gauss1,
    std::int64_t angmom_gauss2,
    std::int64_t angmom_gauss3
)
    : idx_l_01_max_ {angmom_gauss0 + angmom_gauss1 + 1}
    , idx_l_23_max_ {angmom_gauss2 + angmom_gauss3 + 1}
{}

ElectronElectronIndexGenerator::ElectronElectronIndexGenerator(const AngularMomenta1D& angmoms)
    : idx_l_01_max_ {angmoms.angmom_0 + angmoms.angmom_1 + 1}
    , idx_l_23_max_ {angmoms.angmom_2 + angmoms.angmom_3 + 1}
{}

auto ElectronElectronIndexGenerator::begin() noexcept -> ElectronElectronIndexIterator
{
    return ElectronElectronIndexIterator {idx_l_01_max_, idx_l_23_max_};
}

auto ElectronElectronIndexGenerator::end() noexcept -> ElectronElectronIndexIterator
{
    auto iterator = ElectronElectronIndexIterator {idx_l_01_max_, idx_l_23_max_};
    iterator.set_to_end();

    return iterator;
}

auto ElectronElectronIndexGenerator::begin() const noexcept -> ElectronElectronIndexIterator
{
    return ElectronElectronIndexIterator {idx_l_01_max_, idx_l_23_max_};
}

auto ElectronElectronIndexGenerator::end() const noexcept -> ElectronElectronIndexIterator
{
    auto iterator = ElectronElectronIndexIterator {idx_l_01_max_, idx_l_23_max_};
    iterator.set_to_end();

    return iterator;
}

}  // namespace elec
