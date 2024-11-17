#pragma once

#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

/*
    The Grid4D class stores elements of type double in a 4D grid, with a contiguous 1D
    vector as the underlying data type.
*/

namespace elec::grid
{

class Grid4D
{
public:
    Grid4D(std::vector<double> data, std::size_t size0, std::size_t size1, std::size_t size2, std::size_t size3)
        : data_ {std::move(data)}
        , size0_ {size0}
        , size1_ {size1}
        , size2_ {size2}
        , size3_ {size3}
        , n_elements_ {size0 * size1 * size2 * size3}
    {
        check_grid_sizes_();
    }

    Grid4D(std::size_t size0, std::size_t size1, std::size_t size2, std::size_t size3)
        : data_ {}
        , size0_ {size0}
        , size1_ {size1}
        , size2_ {size2}
        , size3_ {size3}
        , n_elements_ {size0 * size1 * size2 * size3}
    {
        data_.resize(n_elements_);
    }

    auto get(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3) const noexcept -> double
    {
        const auto idx = index_(i0, i1, i2, i3);
        return data_[idx];
    }

    void set(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3, double value) noexcept
    {
        const auto idx = index_(i0, i1, i2, i3);
        data_[idx] = value;
    }

    auto get_at(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3) const -> double
    {
        const auto idx = index_(i0, i1, i2, i3);

        if (idx >= n_elements_) {
            throw_index_bounds_error_(i0, i1, i2, i3);
        }

        return data_[idx];
    }

    void set_at(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3, double value)
    {
        const auto idx = index_(i0, i1, i2, i3);

        if (idx >= n_elements_) {
            throw_index_bounds_error_(i0, i1, i2, i3);
        }

        data_[idx] = value;
    }

    auto sizes() const noexcept -> std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>
    {
        return {size0_, size1_, size2_, size3_};
    }

private:
    std::vector<double> data_;
    std::size_t size0_;
    std::size_t size1_;
    std::size_t size2_;
    std::size_t size3_;
    std::size_t n_elements_;

    auto index_(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3) const noexcept -> std::size_t
    {
        return i3 + i2 * size3_ + i1 * size2_ * size3_ + i0 * size1_ * size2_ * size3_;
    }

    void check_grid_sizes_() const
    {
        const auto actual_n_elements = data_.size();
        const auto expected_n_elements = size0_ * size1_ * size2_ * size3_;

        if (actual_n_elements != expected_n_elements) {
            auto err_msg = std::stringstream {};
            err_msg << "Attempt to create Grid4D instance with inconsistent sizes.\n";
            err_msg << "actual number of elements in data: " << actual_n_elements << '\n';
            err_msg << "expected number of elements based on sizes passed: " << expected_n_elements << '\n';

            throw std::runtime_error {err_msg.str()};
        }
    }

    void throw_index_bounds_error_(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3) const
    {
        auto err_msg = std::stringstream {};
        err_msg << "(" << i0 << ", " << i1 << ", " << i2 << ", " << i3 << ") is outside of the bounds of the grid.\n";
        err_msg << "size0 = " << size0_ << '\n';
        err_msg << "size1 = " << size1_ << '\n';
        err_msg << "size2 = " << size2_ << '\n';
        err_msg << "size3 = " << size3_ << '\n';

        throw std::runtime_error {err_msg.str()};
    }
};

}  // namespace elec::grid
