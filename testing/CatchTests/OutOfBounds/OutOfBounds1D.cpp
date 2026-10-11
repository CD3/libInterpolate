#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <libInterpolate/AnyInterpolator.hpp>
#include <libInterpolate/Interpolators/_1D/CubicSplineInterpolator.hpp>
#include <libInterpolate/Interpolators/_1D/LinearInterpolator.hpp>
#include <libInterpolate/Interpolators/_1D/MonotonicInterpolator.hpp>
#include <vector>

#include "OutOfBoundsTestUtils.hpp"
using namespace Catch;

// Data is y = 2x + 1 on [1, 4]. It is never 0 inside the domain, so a 0
// result is always the out-of-bounds value and never a real data value. All
// three 1D methods reproduce a straight line exactly.
TEMPLATE_PRODUCT_TEST_CASE("1D out-of-bounds queries return 0 by default",
                           "[out-of-bounds][1D]",
                           (_1D::LinearInterpolator, _1D::CubicSplineInterpolator,
                            _1D::MonotonicInterpolator),
                           (float, double)) {
    using Real = typename TestType::VectorType::Scalar;
    using L = std::numeric_limits<Real>;
    auto f = [](Real x) { return 2 * x + 1; };
    std::vector<Real> x{1, 2, 3, 4}, y{3, 5, 7, 9};
    const Real lo = x.front(), hi = x.back();
    const Real oob = oob::default_value<Real>();

    TestType interp;
    interp.setData(x, y);

    SECTION("below the lower bound") {
        CHECK_OOB(interp(lo - 1), oob);
        CHECK_OOB(interp(oob::next_below(lo)), oob);
        CHECK_OOB(interp(Real(-1e30)), oob);
        CHECK_OOB(interp(L::lowest()), oob);
        CHECK_OOB(interp(-L::infinity()), oob);
    }
    SECTION("above the upper bound") {
        CHECK_OOB(interp(hi + 1), oob);
        CHECK_OOB(interp(oob::next_above(hi)), oob);
        CHECK_OOB(interp(Real(1e30)), oob);
        CHECK_OOB(interp(L::max()), oob);
        CHECK_OOB(interp(L::infinity()), oob);
    }
    SECTION("boundary knots are inside the domain") {
        CHECK(interp(lo) == Approx(f(lo)));
        CHECK(interp(hi) == Approx(f(hi)));
        CHECK(interp(oob::next_above(lo)) == Approx(f(lo)));
        CHECK(interp(oob::next_below(hi)) == Approx(f(hi)));
    }
    SECTION("batch() matches operator() out of bounds") {
        std::vector<Real> xs{lo - 1, oob::next_below(lo), lo, Real(2.5), hi,
                             oob::next_above(hi), hi + 1};
        auto ys = interp.batch(xs);
        REQUIRE(ys.size() == xs.size());
        CHECK_OOB(ys[0], oob);
        CHECK_OOB(ys[1], oob);
        CHECK(ys[2] == Approx(f(lo)));
        CHECK(ys[3] == Approx(f(Real(2.5))));
        CHECK(ys[4] == Approx(f(hi)));
        CHECK_OOB(ys[5], oob);
        CHECK_OOB(ys[6], oob);
    }
    SECTION("through _1D::AnyInterpolator") {
        _1D::AnyInterpolator<Real> any = TestType();
        any.setData(static_cast<int>(x.size()), x.data(), y.data());
        CHECK_OOB(any(lo - 1), oob);
        CHECK_OOB(any(hi + 1), oob);
        CHECK(any(lo) == Approx(f(lo)));
        CHECK(any(hi) == Approx(f(hi)));
    }
    SECTION("through setUnsafeDataReference") {
        TestType ref;
        ref.setUnsafeDataReference(x.size(), x.data(), y.data());
        CHECK_OOB(ref(lo - 1), oob);
        CHECK_OOB(ref(hi + 1), oob);
        CHECK(ref(lo) == Approx(f(lo)));
    }
    SECTION("copies keep the out-of-bounds behavior") {
        TestType copy(interp);
        TestType assigned;
        assigned = interp;
        CHECK_OOB(copy(lo - 1), oob);
        CHECK_OOB(assigned(hi + 1), oob);
    }
    SECTION("NaN query propagates NaN (current behavior, not a bounds check)") {
        // NaN fails both `x < lo` and `x > hi`, so it reaches the interpolation
        // formula. If a sentinel mode treats NaN as out of bounds, change this
        // on purpose.
        CHECK_THAT(interp(L::quiet_NaN()), Matchers::IsNaN());
    }
}

TEMPLATE_PRODUCT_TEST_CASE(
    "1D descending x data: every query is treated as out of bounds",
    "[out-of-bounds][1D][unsupported-input]",
    (_1D::LinearInterpolator, _1D::CubicSplineInterpolator,
     _1D::MonotonicInterpolator),
    (float, double)) {
    // Pins current behavior. Descending data is not validated. The bounds
    // check uses x.front() and x.back(), so no x satisfies
    // front <= x <= back.
    using Real = typename TestType::VectorType::Scalar;
    std::vector<Real> x{4, 3, 2, 1}, y{9, 7, 5, 3};
    TestType interp;
    interp.setData(x, y);
    for (Real q : {Real(0), Real(1), Real(2.5), Real(4), Real(5)})
        CHECK_OOB(interp(q), oob::default_value<Real>());
}

TEMPLATE_TEST_CASE("CubicSpline derivative/integral outside the domain",
                   "[out-of-bounds][1D][spline]", float, double) {
    using Real = TestType;
    std::vector<Real> x{1, 2, 3, 4}, y{1, 4, 9, 16};  // y = x^2
    _1D::CubicSplineInterpolator<Real> interp;
    interp.setData(x, y);
    const Real oob = oob::default_value<Real>();

    SECTION("derivative above the upper bound is 0") {
        CHECK_OOB(interp.derivative(oob::next_above(Real(4))), oob);
        CHECK_OOB(interp.derivative(Real(5)), oob);
        CHECK_OOB(interp.derivative(std::numeric_limits<Real>::infinity()),
                  oob);
    }
    SECTION("derivative at the boundary knots is nonzero") {
        CHECK(interp.derivative(Real(1)) != 0);
        CHECK(interp.derivative(Real(4)) != 0);
    }
    SECTION("integral limits are clamped to the domain") {
        const Real full = interp.integral(1, 4);
        CHECK(full != 0);
        CHECK(interp.integral(0, 5) == Approx(full));
        CHECK(interp.integral(-std::numeric_limits<Real>::infinity(),
                              std::numeric_limits<Real>::infinity()) ==
              Approx(full));
        CHECK(interp.integral(5, 0) == Approx(-full));
        CHECK(interp.integral(Real(3.5), 10) ==
              Approx(interp.integral(Real(3.5), 4)));
    }
    SECTION("integral over a range that only touches the domain is 0") {
        CHECK_OOB(interp.integral(0, 1), oob);
        CHECK_OOB(interp.integral(4, 5), oob);
    }
}

// Known bugs. [!shouldfail] passes while the bug is present and fails once it
// is fixed, which flags that the tag should be removed.
TEMPLATE_TEST_CASE("KNOWN BUG: CubicSpline derivative extrapolates below x[0]",
                   "[out-of-bounds][1D][spline][!shouldfail]", float, double) {
    std::vector<TestType> x{1, 2, 3, 4}, y{1, 4, 9, 16};
    _1D::CubicSplineInterpolator<TestType> interp(x, y);
    CHECK_OOB(interp.derivative(TestType(0)), oob::default_value<TestType>());
    CHECK_OOB(interp.derivative(oob::next_below(TestType(1))),
              oob::default_value<TestType>());
}

TEMPLATE_TEST_CASE(
    "KNOWN BUG: CubicSpline integral entirely below x[0] is nonzero",
    "[out-of-bounds][1D][spline][!shouldfail]", float, double) {
    std::vector<TestType> x{1, 2, 3, 4}, y{1, 4, 9, 16};
    _1D::CubicSplineInterpolator<TestType> interp(x, y);
    CHECK_OOB(interp.integral(-5, -2), oob::default_value<TestType>());
    CHECK_OOB(interp.integral(TestType(-2), TestType(0.5)),
              oob::default_value<TestType>());
}

// Hidden ([.]) because it reads past the end of the coefficient arrays: an
// Eigen assertion in Debug builds, a heap overflow in Release. Run it on
// purpose with `"[known-crash]"`.
TEST_CASE("KNOWN CRASH: CubicSpline integral entirely above x[N-1]",
          "[.][known-crash][out-of-bounds][1D][spline]") {
    std::vector<double> x{1, 2, 3, 4}, y{1, 4, 9, 16};
    _1D::CubicSplineInterpolator<double> interp(x, y);
    CHECK_OOB(interp.integral(5, 9), 0.0);
}
