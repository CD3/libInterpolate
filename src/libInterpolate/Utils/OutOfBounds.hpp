#pragma once

/** @file OutOfBounds.hpp
 * @brief Out-of-bounds settings shared by the 1D and 2D interpolators.
 */

#include <limits>
#include <type_traits>

namespace libInterpolate {
namespace Utils {

/**
 * Detects whether an interpolator can extrapolate. Interpolators opt in with
 * `static constexpr bool supports_extrapolation = true;`.
 */
template <typename T, typename = void>
struct supports_extrapolation : std::false_type {};
template <typename T>
struct supports_extrapolation<T,
                              std::void_t<decltype(T::supports_extrapolation)>>
    : std::integral_constant<bool, T::supports_extrapolation> {};

/**
 * Controls what an interpolator returns for queries outside of its data. By
 * default it returns 0 and does not extrapolate.
 *
 * The 1D and 2D interpolator base classes inherit from this class, so these
 * methods are available on every interpolator.
 */
template <class Derived, typename Real>
class OutOfBoundsSettings {
   public:
    /**
     * Return v for queries outside of the data (the default is 0). This turns
     * extrapolation off.
     */
    void setOutOfBoundsValue(Real v) {
        m_value = v;
        m_extrapolate = false;
    }

    /**
     * Return quiet NaN for queries outside of the data. This turns
     * extrapolation off.
     */
    void setOutOfBoundsNaN() {
        static_assert(std::numeric_limits<Real>::has_quiet_NaN,
                      "setOutOfBoundsNaN() requires a type with a quiet NaN. "
                      "Use setOutOfBoundsValue() instead.");
        setOutOfBoundsValue(std::numeric_limits<Real>::quiet_NaN());
    }

    /**
     * Extrapolate for queries outside of the data. Call setOutOfBoundsValue()
     * to turn it off again.
     *
     * Not every interpolator can extrapolate. Calling this on one that cannot
     * is a compile error.
     */
    void enableExtrapolation() {
        static_assert(supports_extrapolation<Derived>::value,
                      "This interpolator does not support extrapolation.");
        m_extrapolate = true;
    }

    Real getOutOfBoundsValue() const { return m_value; }
    bool extrapolationEnabled() const { return m_extrapolate; }

   private:
    Real m_value = Real(0);
    bool m_extrapolate = false;
};

}  // namespace Utils
}  // namespace libInterpolate
