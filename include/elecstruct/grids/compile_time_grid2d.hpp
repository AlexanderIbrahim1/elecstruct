#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

namespace elec::grid
{

template <typename T, std::size_t SIZE0, std::size_t SIZE1>
class CompileTimeGrid2D
{
public:
    using Data = std::array<T, SIZE0 * SIZE1>;

    constexpr CompileTimeGrid2D(const Data& data)
        : data_ {data}
    {}

    constexpr auto at_unsigned(std::size_t i0, std::size_t i1) const noexcept -> const T&
    {
        return data_[index(i0, i1)];
    }

    constexpr auto at_unsigned(std::size_t i0, std::size_t i1) noexcept -> T&
    {
        return data_[index(i0, i1)];
    }

    constexpr auto at(std::int64_t i0, std::int64_t i1) const noexcept -> const T&
    {
        const auto i0_ = static_cast<std::size_t>(i0);
        const auto i1_ = static_cast<std::size_t>(i1);
        return data_[index(i0_, i1_)];
    }

    constexpr auto at(std::int64_t i0, std::int64_t i1) noexcept -> T&
    {
        const auto i0_ = static_cast<std::size_t>(i0);
        const auto i1_ = static_cast<std::size_t>(i1);
        return data_[index(i0_, i1_)];
    }

    constexpr auto size0() const noexcept -> std::size_t
    {
        return SIZE0;
    }

    constexpr auto size1() const noexcept -> std::size_t
    {
        return SIZE1;
    }

private:
    constexpr auto index(std::size_t i0, std::size_t i1) const noexcept -> std::size_t
    {
        return i1 + SIZE0 * i0;
    }

    Data data_;
};

}  // namespace elec::grid
