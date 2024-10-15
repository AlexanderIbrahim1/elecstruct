#pragma once

#include <concepts>
#include <iomanip>
#include <ostream>

/*
    This header file contains the Cartesian3D type, which represents a point in 3D
    Cartesian space.
*/

namespace coord
{

template <typename T>
concept Numeric = std::integral<T> || std::floating_point<T>;

constexpr static auto CARTESIAN3D_OSTREAM_PRECISION = int {14};
constexpr static auto ALMOST_EQUALS_DISTANCE_SQ_TOLERANCE = double {1.0e-6};

struct Cartesian3D
{
    double x;
    double y;
    double z;

    auto operator+=(const Cartesian3D& other) noexcept -> Cartesian3D&
    {
        x += other.x;
        y += other.y;
        z += other.z;

        return *this;
    }

    auto operator+() const noexcept -> Cartesian3D
    {
        return Cartesian3D {x, y, z};
    }

    auto operator-=(const Cartesian3D& other) noexcept -> Cartesian3D&
    {
        x -= other.x;
        y -= other.y;
        z -= other.z;

        return *this;
    }

    auto operator-() const noexcept -> Cartesian3D
    {
        return Cartesian3D {-x, -y, -z};
    }

    template <Numeric Number>
    auto operator*=(Number other) noexcept -> Cartesian3D&
    {
        const auto other_double = static_cast<double>(other);

        x *= other_double;
        y *= other_double;
        z *= other_double;

        return *this;
    }

    template <Numeric Number>
    auto operator/=(Number other) -> Cartesian3D&
    {
        const auto other_double = static_cast<double>(other);

        x /= other_double;
        y /= other_double;
        z /= other_double;

        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const Cartesian3D& point)
    {
        os << "(";

        if (point.x >= 0.0) {
            os << ' ';
        }
        os << std::fixed << std::setprecision(CARTESIAN3D_OSTREAM_PRECISION) << point.x;
        os << ", ";

        if (point.y >= 0.0) {
            os << ' ';
        }
        os << std::fixed << std::setprecision(CARTESIAN3D_OSTREAM_PRECISION) << point.y;
        os << ", ";

        if (point.z >= 0.0) {
            os << ' ';
        }
        os << std::fixed << std::setprecision(CARTESIAN3D_OSTREAM_PRECISION) << point.z;

        os << ")";

        return os;
    }
};

inline auto operator+(Cartesian3D left, const Cartesian3D& right) noexcept -> Cartesian3D
{
    left += right;
    return left;
}

inline auto operator-(Cartesian3D left, const Cartesian3D& right) noexcept -> Cartesian3D
{
    left -= right;
    return left;
}

template <Numeric Number>
auto operator*(Cartesian3D point, Number other) noexcept -> Cartesian3D
{
    point *= other;
    return point;
}

template <Numeric Number>
auto operator*(Number other, Cartesian3D point) noexcept -> Cartesian3D
{
    point *= other;
    return point;
}

template <Numeric Number>
auto operator/(Cartesian3D point, Number other) -> Cartesian3D
{
    point /= other;
    return point;
}

inline auto almost_equals(
    const Cartesian3D& point0,
    const Cartesian3D& point1,
    double distance_sq_tolerance = ALMOST_EQUALS_DISTANCE_SQ_TOLERANCE
) -> bool
{
    const auto dx = point1.x - point0.x;
    const auto dy = point1.y - point0.y;
    const auto dz = point1.z - point0.z;
    const auto distance_sq = dx * dx + dy * dy + dz * dz;

    return distance_sq < distance_sq_tolerance;
}

}  // namespace coord
