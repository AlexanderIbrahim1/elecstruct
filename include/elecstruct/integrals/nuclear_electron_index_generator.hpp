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

class NuclearElectronIndexIterator
{
public:
    // clang-format off
    using iterator_category      = std::forward_iterator_tag;
    using value_type             = std::tuple<std::int64_t, std::int64_t, std::int64_t>;
    using size_type              = std::size_t;
    using difference_type        = std::ptrdiff_t;
    using pointer                = value_type*;
    using const_pointer          = const value_type*;
    using iterator               = pointer;
    using const_iterator         = const_pointer;
    // clang-format on

    constexpr explicit NuclearElectronIndexIterator(std::int64_t angmom_0, std::int64_t angmom_1)
        : idx_l_max_ {angmom_0 + angmom_1 + 1}
    {}

    constexpr auto operator*() const noexcept -> value_type
    {
        return {idx_l_, idx_r_, idx_i_};
    }

    constexpr auto operator->() const noexcept -> value_type
    {
        return **this;
    }

    constexpr auto operator++() noexcept -> NuclearElectronIndexIterator&
    {
        ++idx_i_;
        if (idx_i_ == idx_i_max_) {
            idx_i_ = 0;
            ++idx_r_;

            idx_i_max_ = update_idx_i_max_(idx_l_, idx_r_);
            if (idx_r_ == idx_r_max_) {
                idx_r_ = 0;
                ++idx_l_;

                idx_r_max_ = update_idx_r_max_(idx_l_);
                idx_i_max_ = update_idx_i_max_(idx_l_, idx_r_);
            }
        }

        return *this;
    }

    constexpr auto operator++(int) noexcept -> NuclearElectronIndexIterator&
    {
        auto temp = *this;
        ++(*this);
        return temp;
    }

    constexpr void set_to_end() noexcept
    {
        idx_l_ = idx_l_max_;
        idx_r_ = 0;
        idx_i_ = 0;

        idx_r_max_ = update_idx_r_max_(idx_l_);
        idx_i_max_ = update_idx_i_max_(idx_l_, idx_r_);
    }

    friend constexpr auto operator==(const NuclearElectronIndexIterator& left, const NuclearElectronIndexIterator& right) noexcept -> bool
    {
        // clang-format off
        return \
            left.idx_l_max_ == right.idx_l_max_ && \
            left.idx_r_max_ == right.idx_r_max_ && \
            left.idx_i_max_ == right.idx_i_max_ && \
            left.idx_l_ == right.idx_l_ && \
            left.idx_r_ == right.idx_r_ && \
            left.idx_i_ == right.idx_i_;
        // clang-format on
    }

    friend constexpr auto operator!=(const NuclearElectronIndexIterator& left, const NuclearElectronIndexIterator& right) noexcept -> bool
    {
        return !(left == right);
    }

private:
    std::int64_t idx_l_max_;
    std::int64_t idx_r_max_ {1};
    std::int64_t idx_i_max_ {1};

    std::int64_t idx_l_ {0};
    std::int64_t idx_r_ {0};
    std::int64_t idx_i_ {0};

    constexpr auto update_idx_i_max_(std::int64_t idx_l, std::int64_t idx_r) const noexcept -> std::int64_t
    {
        return static_cast<std::int64_t>((idx_l - 2 * idx_r) / 2) + 1;
    }

    constexpr auto update_idx_r_max_(std::int64_t idx_l) const noexcept -> std::int64_t
    {
        return static_cast<std::int64_t>(idx_l / 2) + 1;
    }
};

}  // namespace elec
