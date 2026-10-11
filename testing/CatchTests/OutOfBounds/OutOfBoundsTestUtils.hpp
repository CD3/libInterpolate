#pragma once
// Shared helpers for the out-of-bounds tests.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <limits>

namespace oob {

// The value a default-configured interpolator returns outside its domain.
// This must stay 0 after extrapolation and sentinel options are added. New
// modes get their own SECTIONs that pass their own expected value (e.g. NaN)
// to CHECK_OOB.
template <typename T>
constexpr T default_value() {
    return T(0);
}

template <typename T>
T next_below(T v) {
    return std::nextafter(v, -std::numeric_limits<T>::infinity());
}
template <typename T>
T next_above(T v) {
    return std::nextafter(v, std::numeric_limits<T>::infinity());
}

// Compares exactly, and treats NaN == NaN as a match, so the same helper
// works for the default (0) and for a future quiet_NaN sentinel.
template <typename T>
bool same_value(T actual, T expected) {
    if (std::isnan(expected)) return std::isnan(actual);
    return actual == expected;
}

}  // namespace oob

// Macro rather than function so Catch reports the caller's line and the query.
#define CHECK_OOB(call, expected)                        \
    do {                                                 \
        const auto oob_actual_ = (call);                 \
        INFO(#call << " = " << oob_actual_);             \
        CHECK(oob::same_value(oob_actual_, (expected))); \
    } while (0)
