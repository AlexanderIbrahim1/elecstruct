#pragma once

#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

/*
    The Grid2D class stores elements of type double in a 2D grid, with a contiguous 1D
    vector as the underlying data type.
*/

namespace elec::grid
{

class Grid2D
{
public:
    constexpr Grid2D(std::vector<double> data, std::size_t size0, std::size_t size1)
        : data_ {std::move(data)}
        , size0_ {size0}
        , size1_ {size1}
        , n_elements_ {size0 * size1}
    {
        check_grid_sizes_();
    }

    constexpr Grid2D(std::size_t size0, std::size_t size1)
        : data_ {}
        , size0_ {size0}
        , size1_ {size1}
        , n_elements_ {size0 * size1}
    {
        data_.resize(n_elements_);
    }

    constexpr auto get(std::size_t i0, std::size_t i1) const noexcept -> double
    {
        const auto idx = index_(i0, i1);
        return data_[idx];
    }

    constexpr void set(std::size_t i0, std::size_t i1, double value) noexcept
    {
        const auto idx = index_(i0, i1);
        data_[idx] = value;
    }

    auto get_at(std::size_t i0, std::size_t i1) const -> double
    {
        const auto idx = index_(i0, i1);

        if (idx >= n_elements_) {
            throw_index_bounds_error_(i0, i1);
        }

        return data_[idx];
    }

    void set_at(std::size_t i0, std::size_t i1, double value)
    {
        const auto idx = index_(i0, i1);

        if (idx >= n_elements_) {
            throw_index_bounds_error_(i0, i1);
        }

        data_[idx] = value;
    }

    constexpr auto sizes() const noexcept -> std::tuple<std::size_t, std::size_t>
    {
        return {size0_, size1_};
    }

private:
    std::vector<double> data_;
    std::size_t size0_;
    std::size_t size1_;
    std::size_t n_elements_;

    constexpr auto index_(std::size_t i0, std::size_t i1) const noexcept -> std::size_t
    {
        return i1 + i0 * size1_;
    }

    constexpr void check_grid_sizes_() const
    {
        const auto actual_n_elements = data_.size();
        const auto expected_n_elements = size0_ * size1_;

        if (actual_n_elements != expected_n_elements) {
            throw std::runtime_error {"Attempt to create Grid2D instance with inconsistent sizes.\n"};
        }
    }

    void throw_index_bounds_error_(std::size_t i0, std::size_t i1) const
    {
        auto err_msg = std::stringstream {};
        err_msg << "(" << i0 << ", " << i1 << ") is outside of the bounds of the grid.\n";
        err_msg << "size0 = " << size0_ << '\n';
        err_msg << "size1 = " << size1_ << '\n';

        throw std::runtime_error {err_msg.str()};
    }
};

}  // namespace elec::grid
