#pragma once

#include <iterator>
#include <cstdint>
#include <tuple>

namespace elec
{

class NuclearElectronIndexIterator
{
public:
    // clang-format off
    using iterator_category = std::forward_iterator_tag;
    using value_type        = std::tuple<std::int64_t, std::int64_t, std::int64_t>;
    using size_type         = std::size_t;
    using difference_type   = std::ptrdiff_t;
    using pointer           = value_type*;
    using const_pointer     = const value_type*;
    using iterator          = pointer;
    using const_iterator    = const_pointer;
    // clang-format on

    explicit NuclearElectronIndexIterator(std::int64_t idx_l_max);

    auto operator*() const noexcept -> value_type;
    auto operator->() const noexcept -> value_type;
    auto operator++() noexcept -> NuclearElectronIndexIterator&;
    auto operator++(int) noexcept -> NuclearElectronIndexIterator;
    bool operator==(const NuclearElectronIndexIterator&) const noexcept;
    friend auto operator!=(const NuclearElectronIndexIterator& left, const NuclearElectronIndexIterator& right) noexcept -> bool;

    void set_to_end() noexcept;

private:
    std::int64_t idx_l_max_;
    std::int64_t idx_r_max_ {1};
    std::int64_t idx_i_max_ {1};

    std::int64_t idx_l_ {0};
    std::int64_t idx_r_ {0};
    std::int64_t idx_i_ {0};
};


class NuclearElectronIndexGenerator
{
public:
    explicit NuclearElectronIndexGenerator(std::int64_t angmom_0, std::int64_t angmom_1);

    auto begin() noexcept -> NuclearElectronIndexIterator;
    auto end() noexcept -> NuclearElectronIndexIterator;
    auto begin() const noexcept -> NuclearElectronIndexIterator;
    auto end() const noexcept -> NuclearElectronIndexIterator;

private:
    std::int64_t idx_l_max_;
};

}  // namespace elec
