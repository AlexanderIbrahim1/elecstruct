#pragma once

#include <cstddef>
#include <unordered_map>

namespace elec
{

auto yoshimine_sort(std::size_t a, std::size_t b, std::size_t c, std::size_t d) noexcept -> std::size_t;

/*
    The TwoElectronIntegralGrid is a wrapper around a hashmap that converts the four indices
    into a Yoshimine index before interfacing with the hashmap.
*/
class TwoElectronIntegralGrid
{
public:
    /*
        Check if these four indices corresponds to a Yoshimine index that has already been set in
        the internal mapping.
    */
    auto exists(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3) const noexcept -> bool;

    /*
        Map the Yoshimine index corresponding to the four basis function indices to the value.
    */
    auto set(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3, double value) noexcept -> void;

    /*
        Get the value mapped to from the  Yoshimine index corresponding to the four basis function indices.
    */
    auto get(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3) const noexcept -> double;

private:
    std::unordered_map<std::size_t, double> integrals_;
};

}  // namespace elec

