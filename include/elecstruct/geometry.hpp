#pragma once

#include <cmath>

#include "elecstruct/cartesian3d.hpp"

namespace coord
{

auto distance_squared(const Cartesian3D& point0, const Cartesian3D point1) noexcept -> double;

auto distance(const Cartesian3D& point0, const Cartesian3D& point1) -> double;

auto norm_squared(const Cartesian3D& point) noexcept -> double;

auto norm(const Cartesian3D& point) -> double;

auto dot_product(const Cartesian3D& point0, const Cartesian3D& point1) -> double;

auto cross_product(const Cartesian3D& point0, const Cartesian3D& point1) -> Cartesian3D;

auto unit_vector(const Cartesian3D& point) -> Cartesian3D;

auto bond_angle(const Cartesian3D& point0, const Cartesian3D& point1, const Cartesian3D& point2) -> double;

auto almost_equals(const Cartesian3D& point0, const Cartesian3D& point1, double abs_tol) -> bool;

}  // namespace coord
