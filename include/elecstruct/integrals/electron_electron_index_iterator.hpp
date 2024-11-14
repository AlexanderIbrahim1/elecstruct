#pragma once

#include <iterator>
#include <cstdint>
#include <tuple>

/*
    NOTE: Using a custom iterator object takes up a lot of additional code, adds some
    complications, might even be a little bit slower, and these loops don't appear in
    enough places in the code for this encapulation to be worth it.

    However, I part of the purpose of writing this code is to improve my general C++
    skills, and so I decided to implement the nested loops in the body of the integral
    for the nuclear-electron attraction using custom iterators. This reduces the amount
    of nesting in the body of the function from 9 to 3.
*/

namespace elec
{

struct ElectronElectronIntegralIndices
{
    std::int64_t idx_l_01;
    std::int64_t idx_r_01;
    std::int64_t idx_l_23;
    std::int64_t idx_r_23;
    std::int64_t idx_i;

    constexpr bool operator==(const ElectronElectronIntegralIndices&) const noexcept = default;
};

struct AngularMomenta1D
{
    std::int64_t angmom_0;
    std::int64_t angmom_1;
    std::int64_t angmom_2;
    std::int64_t angmom_3;
};

class ElectronElectronIndexIterator
{
public:
    // clang-format off
    using iterator_category = std::forward_iterator_tag;
    using value_type        = ElectronElectronIntegralIndices;
    using size_type         = std::size_t;
    using difference_type   = std::ptrdiff_t;
    using pointer           = value_type*;
    using const_pointer     = const value_type*;
    using iterator          = pointer;
    using const_iterator    = const_pointer;
    // clang-format on

    constexpr explicit ElectronElectronIndexIterator(
        std::int64_t idx_l_01_max,
        std::int64_t idx_l_23_max
    )
        : idx_l_01_max_ {idx_l_01_max}
        , idx_l_23_max_ {idx_l_23_max}
    {}

    constexpr auto operator*() const noexcept -> value_type
    {
        return ElectronElectronIntegralIndices {idx_l_01_, idx_r_01_, idx_l_23_, idx_r_23_, idx_i_};
    }

    constexpr auto operator->() const noexcept -> value_type
    {
        return **this;
    }

    constexpr auto operator++() noexcept -> ElectronElectronIndexIterator&
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

    constexpr auto operator++(int) noexcept -> ElectronElectronIndexIterator
    {
        auto temp = *this;
        ++(*this);
        return temp;
    }

    constexpr void set_to_end() noexcept
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

    constexpr bool operator==(const ElectronElectronIndexIterator&) const noexcept = default;

    friend constexpr auto operator!=(const ElectronElectronIndexIterator& left, const ElectronElectronIndexIterator& right) noexcept -> bool
    {
        return !(left == right);
    }

private:
    std::int64_t idx_l_01_max_;
    std::int64_t idx_l_23_max_;
    std::int64_t idx_r_01_max_ {1};
    std::int64_t idx_r_23_max_ {1};
    std::int64_t idx_i_max_ {1};

    std::int64_t idx_l_01_ {0};
    std::int64_t idx_l_23_ {0};
    std::int64_t idx_r_01_ {0};
    std::int64_t idx_r_23_ {0};
    std::int64_t idx_i_ {0};

    constexpr auto calculate_idx_i_max_(
        std::int64_t idx_l_01,
        std::int64_t idx_l_23,
        std::int64_t idx_r_01,
        std::int64_t idx_r_23
    ) const noexcept -> std::int64_t
    {
        return static_cast<std::int64_t>((idx_l_01 + idx_l_23 - 2 * (idx_r_01 + idx_r_23)) / 2) + 1;
    }

    constexpr auto calculate_idx_r_max_(std::int64_t idx_l) const noexcept -> std::int64_t
    {
        return static_cast<std::int64_t>(idx_l / 2) + 1;
    }
};


class ElectronElectronIndexGenerator
{
public:
    constexpr explicit ElectronElectronIndexGenerator(
        std::int64_t angmom_gauss0,
        std::int64_t angmom_gauss1,
        std::int64_t angmom_gauss2,
        std::int64_t angmom_gauss3
    )
        : idx_l_01_max_ {angmom_gauss0 + angmom_gauss1 + 1}
        , idx_l_23_max_ {angmom_gauss2 + angmom_gauss3 + 1}
    {}

    constexpr explicit ElectronElectronIndexGenerator(const AngularMomenta1D& angmoms)
        : idx_l_01_max_ {angmoms.angmom_0 + angmoms.angmom_1 + 1}
        , idx_l_23_max_ {angmoms.angmom_2 + angmoms.angmom_3 + 1}
    {}

    constexpr auto begin() noexcept -> ElectronElectronIndexIterator
    {
        return ElectronElectronIndexIterator {idx_l_01_max_, idx_l_23_max_};
    }

    constexpr auto end() noexcept -> ElectronElectronIndexIterator
    {
        auto iterator = ElectronElectronIndexIterator {idx_l_01_max_, idx_l_23_max_};
        iterator.set_to_end();

        return iterator;
    }

    constexpr auto begin() const noexcept -> ElectronElectronIndexIterator
    {
        return ElectronElectronIndexIterator {idx_l_01_max_, idx_l_23_max_};
    }

    constexpr auto end() const noexcept -> ElectronElectronIndexIterator
    {
        auto iterator = ElectronElectronIndexIterator {idx_l_01_max_, idx_l_23_max_};
        iterator.set_to_end();

        return iterator;
    }

private:
    std::int64_t idx_l_01_max_;
    std::int64_t idx_l_23_max_;
};

}  // namespace elec
