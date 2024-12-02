#pragma once

#include <iterator>
#include <cstdint>

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

    explicit ElectronElectronIndexIterator(
        std::int64_t idx_l_01_max,
        std::int64_t idx_l_23_max
    );

    auto operator*() const noexcept -> value_type;
    auto operator->() const noexcept -> value_type;
    auto operator++() noexcept -> ElectronElectronIndexIterator&;
    auto operator++(int) noexcept -> ElectronElectronIndexIterator;
    bool operator==(const ElectronElectronIndexIterator&) const noexcept;
    friend auto operator!=(const ElectronElectronIndexIterator& left, const ElectronElectronIndexIterator& right) noexcept -> bool;
    void set_to_end() noexcept;

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
};


class ElectronElectronIndexGenerator
{
public:
    explicit ElectronElectronIndexGenerator(
        std::int64_t angmom_gauss0,
        std::int64_t angmom_gauss1,
        std::int64_t angmom_gauss2,
        std::int64_t angmom_gauss3
    );

    explicit ElectronElectronIndexGenerator(const AngularMomenta1D& angmoms);

    auto begin() noexcept -> ElectronElectronIndexIterator;
    auto end() noexcept -> ElectronElectronIndexIterator;
    auto begin() const noexcept -> ElectronElectronIndexIterator;
    auto end() const noexcept -> ElectronElectronIndexIterator;

private:
    std::int64_t idx_l_01_max_;
    std::int64_t idx_l_23_max_;
};

}  // namespace elec
