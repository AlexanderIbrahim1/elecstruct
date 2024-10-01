#pragma once

#include <cmath>

#include "cartesian3d.hpp"

namespace coord
{

auto distance_squared(const Cartesian3D& point0, const Cartesian3D point1) noexcept -> double
{
    const auto dx = point0.x - point1.x;
    const auto dy = point0.y - point1.y;
    const auto dz = point0.z - point1.z;

    return dx * dx + dy * dy + dz * dz;
}

auto distance(const Cartesian3D point0, const Cartesian3D point1) -> double
{
    return std::sqrt(distance_squared(point0, point1));
}

auto norm_squared(const Cartesian3D point) noexcept -> double
{
    const auto x_sq = point.x * point.x;
    const auto y_sq = point.y * point.y;
    const auto z_sq = point.z * point.z;

    return x_sq + y_sq + z_sq;
}

auto norm(const Cartesian3D point) -> double
{
    return std::sqrt(norm_squared(point));
}

auto dot_product(const Cartesian3D point0, const Cartesian3D point1) -> double
{
    const auto dot_x = point0.x * point1.x;
    const auto dot_y = point0.y * point1.y;
    const auto dot_z = point0.z * point1.z;

    return dot_x + dot_y + dot_z;
}

auto unit_vector(const Cartesian3D point) -> Cartesian3D
{
    return point / norm(point);
}

}  // namespace coord
