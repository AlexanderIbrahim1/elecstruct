#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>

/*
    A reproduction of the iemplementation of the fast Boys algorithm given by https://doi.org/10.1063/5.0062444
*/


namespace elec
{

namespace impl_constants
{

constexpr auto BOYS_ORDER_UPPER_BOUND_ = std::size_t {13};

constexpr auto LARGE_CUTOFF_ = double {0.45425955121971775e01};
constexpr auto LARGE_SQRT_PI_O_2_ = double {0.886226925452758014};

constexpr auto SMALL_TOL_ = double {1.0e-3};
constexpr auto SMALL_RZZ_ = double {-0.96321934290343840e01};
constexpr auto SMALL_RFACT_ = double {0.15247844519077540e05};
constexpr auto SMALL_RWW_ = double {0.18995875677635889e-04};

constexpr auto T_DATA_ = std::array<double, 12> {
    0.20000000000000000e01,
    0.66666666666666663e00,
    0.40000000000000002e00,
    0.28571428571428570e00,
    0.22222222222222221e00,
    0.18181818181818182e00,
    0.15384615384615385e00,
    0.13333333333333333e00,
    0.11764705882352941e00,
    0.10526315789473684e00,
    0.95238095238095233e-01,
    0.86956521739130432e-01
};

constexpr auto ZZ_DATA_ = std::array<std::complex<double>, 10> {
    std::complex {0.64304020652330500e01, 0.18243694739308491e02},
    std::complex {0.64304020652330500e01, -0.18243694739308491e02},
    std::complex {-0.12572081889410178e01, 0.14121366415342502e02},
    std::complex {-0.12572081889410178e01, -0.14121366415342502e02},
    std::complex {-0.54103079551670268e01, 0.10457909575828442e02},
    std::complex {-0.54103079551670268e01, -0.10457909575828442e02},
    std::complex {-0.78720025594983341e01, 0.69309284623985663e01},
    std::complex {-0.78720025594983341e01, -0.69309284623985663e01},
    std::complex {-0.92069621609035313e01, 0.34559308619699376e01},
    std::complex {-0.92069621609035313e01, -0.34559308619699376e01}
};

constexpr auto FACT_DATA_ = std::array<std::complex<double>, 10> {
    std::complex {0.13249210991966042e-02, 0.91787356295447745e-03},
    std::complex {0.13249210991966042e-02, -0.91787356295447745e-03},
    std::complex {0.55545905103006735e-01, -0.35151540664451613e01},
    std::complex {0.55545905103006735e-01, 0.35151540664451613e01},
    std::complex {-0.11456407675096416e03, 0.19213789620924834e03},
    std::complex {-0.11456407675096416e03, -0.19213789620924834e03},
    std::complex {0.20915556220686653e04, -0.15825742912360638e04},
    std::complex {0.20915556220686653e04, 0.15825742912360638e04},
    std::complex {-0.94779394228935325e04, 0.30814443710192086e04},
    std::complex {-0.94779394228935325e04, -0.30814443710192086e04}
};

// WW_DATA: list[complex] = [
constexpr auto WW_DATA_ = std::array<std::complex<double>, 10> {
    std::complex {-0.83418049867878959e-08, -0.70958810331788253e-08},
    std::complex {-0.83418050437598581e-08, 0.70958810084577824e-08},
    std::complex {0.82436739552884774e-07, -0.27704117936134414e-06},
    std::complex {0.82436739547688584e-07, 0.27704117938414886e-06},
    std::complex {0.19838416382728666e-05, 0.78321058613942770e-06},
    std::complex {0.19838416382681279e-05, -0.78321058613180811e-06},
    std::complex {-0.47372729839268780e-05, 0.58076919074212929e-05},
    std::complex {-0.47372729839287016e-05, -0.58076919074154416e-05},
    std::complex {-0.68186014282131608e-05, -0.13515261354290787e-04},
    std::complex {-0.68186014282138385e-05, 0.13515261354295612e-04}
};

}  // namespace elec::impl_constants

namespace impl_boys
{

/*
    NOTE: the check to make sure that n <= 12 should be done outside of this function.
*/
inline auto boys_fast_large_(double x, std::size_t n) -> double
{
    namespace constants = elec::impl_constants;

    const auto halfy = std::exp(-x) / 2.0;
    const auto yy = std::sqrt(x);

    std::array<double, constants::BOYS_ORDER_UPPER_BOUND_> values;
    values[0] = constants::LARGE_SQRT_PI_O_2_ * std::erf(yy) / yy;

    if (n == 0) {
        return values[0];
    }

    for (std::size_t i {1}; i < constants::BOYS_ORDER_UPPER_BOUND_; ++i) {
        const auto i_shifted = static_cast<double>(i) - 0.5;
        values[i] = (i_shifted * values[i - 1] - halfy) / x;
        if (i == n) {
            return values[i];
        }
    }

    // technically this function only accepts 0 <= n <= constants::BOYS_ORDER_UPPER_BOUND_, and
    // so should return from the loop; but the compiler raises a warning if we don't write this
    return values[constants::BOYS_ORDER_UPPER_BOUND_ - 1];
}

inline auto boys_fast_small_rtmp_(double x, double y) -> double
{
    namespace cts = elec::impl_constants;

    auto rtmp_ = double {0.0};
    for (std::size_t i {1}; i < 11; i += 2) {
        const auto component = cts::WW_DATA_[i] * (1.0 - cts::FACT_DATA_[i] * y) / (x + cts::ZZ_DATA_[i]);
        rtmp_ += component.real();
    }

    return rtmp_;
}

inline auto boys_fast_small_tmp_(double x, double y) -> double
{
    namespace cts = elec::impl_constants;

    const auto q = x + cts::SMALL_RZZ_;

    if (std::fabs(x + cts::SMALL_RZZ_) >= cts::SMALL_TOL_) {
        return cts::SMALL_RWW_ * (1.0 - cts::SMALL_RFACT_ * y) / q;
    }
    else {
        const auto p0 = double {1.0};
        const auto p1 = (1.0 / 2.0) * q;
        const auto p2 = p1 * q / 3.0;
        const auto p3 = p2 * q / 4.0;
        const auto p4 = p3 * q / 5.0;

        return cts::SMALL_RWW_ * (p0 + p1 + p2 + p3 + p4);
    }
}

inline auto boys_fast_small_(double x, std::size_t n) -> double
{
    namespace cts = elec::impl_constants;

    const auto y = std::exp(-x);
    const auto rtmp = boys_fast_small_rtmp_(x, y);
    const auto tmp = boys_fast_small_tmp_(x, y);

    std::array<double, cts::BOYS_ORDER_UPPER_BOUND_> values;
    values[12] = 2.0 * rtmp + tmp;
    
    if (n == 12) {
        return values[12];
    }
    else {
        const auto halfy = y / 2.0;

        // need to do an index shift because the terminating condition ends up with i == -1,
        // but `i` is an instance of `std::size_t`
        //
        // for (std::size_t i {11}; i >= 0; --i) {
        for (std::size_t j {12}; j >= 1; --j) {
            std::size_t i = j - 1;
            values[i] = (x * values[i + 1] + halfy) * cts::T_DATA_[i];
            if (i == n) {
                return values[i];
            }
        }

        // technically this function only accepts 0 <= n <= constants::BOYS_ORDER_UPPER_BOUND_, and
        // so should return from the loop; but the compiler raises a warning if we don't write this
        return values[0];
    }
}

}  // namespace elec::impl_boys

inline auto boys_fast(double x, std::size_t n) -> double
{
    namespace cts = elec::impl_constants;
    namespace fun = elec::impl_boys;

    if (n >= cts::BOYS_ORDER_UPPER_BOUND_) {
        throw std::runtime_error {
            "This implementation only works for the Boys function up to and including order 12."
        };
    }

    if (std::fabs(x) >= cts::LARGE_CUTOFF_) {
        return fun::boys_fast_large_(x, n);
    } else {
        return fun::boys_fast_small_(x, n);
    }
}

}  // namespace elec
