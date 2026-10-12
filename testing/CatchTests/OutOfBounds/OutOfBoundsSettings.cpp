#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <libInterpolate/AnyInterpolator.hpp>
#include <libInterpolate/Interpolate.hpp>
#include <type_traits>
#include <vector>

#include "OutOfBoundsTestUtils.hpp"
using namespace Catch;

// Tests for setOutOfBoundsValue(), setOutOfBoundsNaN(), and
// enableExtrapolation(). The default behavior (return 0) is tested in
// OutOfBounds1D.cpp and OutOfBounds2D.cpp.

namespace {
// z = x + 2y + 100 on x in {1,2,3}, y in {10,20,30,40}. Never 0 or -1 in or
// near the domain, so the sentinel values used below are unambiguous. Every 2D
// method except NearestNeighbor reproduces a plane exactly, including when
// extrapolating.
template <typename Real>
struct Plane {
    std::vector<Real> x, y, z;
    static Real f(Real x, Real y) { return x + 2 * y + 100; }
    Plane() {
        for (Real a : {Real(1), Real(2), Real(3)})
            for (Real b : {Real(10), Real(20), Real(30), Real(40)}) {
                x.push_back(a);
                y.push_back(b);
                z.push_back(f(a, b));
            }
    }
};

// Relative tolerance loose enough for float.
template <typename Real>
Approx approx(Real v) {
    return Approx(v).epsilon(1e-4);
}
}  // namespace

using libInterpolate::Utils::supports_extrapolation;
static_assert(supports_extrapolation<_1D::LinearInterpolator<double>>::value);
static_assert(
    supports_extrapolation<_1D::CubicSplineInterpolator<double>>::value);
static_assert(supports_extrapolation<_1D::MonotonicInterpolator<double>>::value);
static_assert(supports_extrapolation<_2D::BilinearInterpolator<double>>::value);
static_assert(supports_extrapolation<_2D::BicubicInterpolator<double>>::value);
static_assert(
    supports_extrapolation<_2D::NearestNeighborInterpolator<double>>::value);
static_assert(
    supports_extrapolation<_2D::ThinPlateSplineInterpolator<double>>::value);
// LinearDelaunay cannot extrapolate, so enableExtrapolation() does not compile
// for it. Its base class is the "Derived" type seen by the settings.
static_assert(!supports_extrapolation<_2D::DelaunayTriangulationInterpolatorBase<
                  _2D::LinearDelaunayTriangleInterpolator<double>>>::value);

TEMPLATE_PRODUCT_TEST_CASE("1D out-of-bounds settings", "[out-of-bounds][1D]",
                           (_1D::LinearInterpolator,
                            _1D::CubicSplineInterpolator,
                            _1D::MonotonicInterpolator),
                           (float, double)) {
    using Real = typename TestType::VectorType::Scalar;
    using L = std::numeric_limits<Real>;
    // y = 2x + 1. All three methods reproduce a straight line exactly, and so
    // do their straight-line extrapolations.
    auto f = [](Real x) { return 2 * x + 1; };
    std::vector<Real> x{1, 2, 3, 4}, y{3, 5, 7, 9};
    TestType interp;
    interp.setData(x, y);

    SECTION("defaults") {
        CHECK(interp.getOutOfBoundsValue() == 0);
        CHECK(!interp.extrapolationEnabled());
    }

    SECTION("custom value") {
        interp.setOutOfBoundsValue(-1);
        CHECK(interp.getOutOfBoundsValue() == -1);
        CHECK_OOB(interp(0), Real(-1));
        CHECK_OOB(interp(oob::next_below(Real(1))), Real(-1));
        CHECK_OOB(interp(oob::next_above(Real(4))), Real(-1));
        CHECK_OOB(interp(L::infinity()), Real(-1));
        CHECK_OOB(interp(-L::infinity()), Real(-1));
        // inside and on the boundary are unaffected
        CHECK(interp(1) == approx(f(1)));
        CHECK(interp(2.5) == approx(f(2.5)));
        CHECK(interp(4) == approx(f(4)));
    }

    SECTION("NaN") {
        interp.setOutOfBoundsNaN();
        CHECK_THAT(interp.getOutOfBoundsValue(), Matchers::IsNaN());
        CHECK_THAT(interp(0), Matchers::IsNaN());
        CHECK_THAT(interp(5), Matchers::IsNaN());
        CHECK(interp(2.5) == approx(f(2.5)));
    }

    SECTION("extrapolation") {
        interp.enableExtrapolation();
        CHECK(interp.extrapolationEnabled());
        CHECK(interp(0) == approx(f(0)));
        CHECK(interp(-10) == approx(f(-10)));
        CHECK(interp(5) == approx(f(5)));
        CHECK(interp(100) == approx(f(100)));
        CHECK(interp(2.5) == approx(f(2.5)));

        SECTION("setOutOfBoundsValue turns it off") {
            interp.setOutOfBoundsValue(-1);
            CHECK(!interp.extrapolationEnabled());
            CHECK_OOB(interp(0), Real(-1));
        }
    }

    SECTION("settings are kept by copies, AnyInterpolator, and batch") {
        interp.setOutOfBoundsValue(-1);
        TestType copy(interp);
        TestType assigned;
        assigned = interp;
        _1D::AnyInterpolator<Real> any(interp);
        CHECK_OOB(copy(0), Real(-1));
        CHECK_OOB(assigned(5), Real(-1));
        CHECK_OOB(any(0), Real(-1));
        auto ys = interp.batch(std::vector<Real>{0, 2.5, 5});
        CHECK_OOB(ys[0], Real(-1));
        CHECK(ys[1] == approx(f(2.5)));
        CHECK_OOB(ys[2], Real(-1));

        interp.enableExtrapolation();
        TestType ecopy(interp);
        _1D::AnyInterpolator<Real> eany(interp);
        CHECK(ecopy(5) == approx(f(5)));
        CHECK(eany(0) == approx(f(0)));
        CHECK(interp.batch(std::vector<Real>{0})[0] == approx(f(0)));
    }
}

TEST_CASE("1D extrapolation continues in a straight line from the ends",
          "[out-of-bounds][1D]") {
    // curved, increasing data
    std::vector<double> x{1, 2, 3, 4}, y{1, 4, 9, 16};

    SECTION("CubicSpline uses the slope at the end points") {
        _1D::CubicSplineInterpolator<double> interp(x, y);
        interp.enableExtrapolation();
        const double s0 = interp.derivative(1), s1 = interp.derivative(4);
        CHECK(interp(0) == Approx(1 - s0));
        CHECK(interp(-9) == Approx(1 - 10 * s0));
        CHECK(interp(5) == Approx(16 + s1));
        CHECK(interp(14) == Approx(16 + 10 * s1));
        // continuous at the boundary
        CHECK(interp(oob::next_below(1.)) == Approx(interp(1)));
        CHECK(interp(oob::next_above(4.)) == Approx(interp(4)));
    }

    SECTION("Monotonic stays monotonic") {
        _1D::MonotonicInterpolator<double> interp(x, y);
        interp.enableExtrapolation();
        // the end slopes are the slopes of the end segments
        CHECK(interp(0) == Approx(1 - 3));
        CHECK(interp(5) == Approx(16 + 7));
        CHECK(interp(oob::next_below(1.)) == Approx(interp(1)));
        CHECK(interp(oob::next_above(4.)) == Approx(interp(4)));
        double last = interp(-10);
        for (double q = -9.5; q <= 15; q += 0.5) {
            CHECK(interp(q) > last);
            last = interp(q);
        }
    }

    SECTION("Monotonic with a flat end segment extrapolates flat") {
        std::vector<double> xf{1, 2, 3, 4}, yf{1, 2, 5, 5};
        _1D::MonotonicInterpolator<double> interp(xf, yf);
        interp.enableExtrapolation();
        CHECK(interp(10) == Approx(5));
    }
}

TEMPLATE_PRODUCT_TEST_CASE(
    "2D out-of-bounds settings", "[out-of-bounds][2D]",
    (_2D::BilinearInterpolator, _2D::BicubicInterpolator,
     _2D::NearestNeighborInterpolator, _2D::ThinPlateSplineInterpolator,
     _2D::LinearDelaunayTriangleInterpolator),
    (float, double)) {
    using Real = typename TestType::VectorType::Scalar;
    using L = std::numeric_limits<Real>;
    Plane<Real> g;
    TestType interp;
    interp.setData(g.x, g.y, g.z);

    SECTION("defaults") {
        CHECK(interp.getOutOfBoundsValue() == 0);
        CHECK(!interp.extrapolationEnabled());
    }

    SECTION("custom value") {
        interp.setOutOfBoundsValue(-1);
        CHECK_OOB(interp(0, 20), Real(-1));
        CHECK_OOB(interp(4, 20), Real(-1));
        CHECK_OOB(interp(2, 5), Real(-1));
        CHECK_OOB(interp(2, 50), Real(-1));
        CHECK_OOB(interp(0, 0), Real(-1));
        CHECK_OOB(interp(L::infinity(), 20), Real(-1));
        // inside and on the boundary (data points) are unaffected
        CHECK(interp(2, 20) == approx(g.f(2, 20)));
        CHECK(interp(1, 10) == approx(g.f(1, 10)));
        CHECK(interp(3, 40) == approx(g.f(3, 40)));
    }

    SECTION("NaN") {
        interp.setOutOfBoundsNaN();
        CHECK_THAT(interp(0, 20), Matchers::IsNaN());
        CHECK_THAT(interp(2, 50), Matchers::IsNaN());
        CHECK(interp(2, 20) == approx(g.f(2, 20)));
    }

    SECTION("settings are kept by copies and AnyInterpolator") {
        interp.setOutOfBoundsValue(-1);
        TestType copy(interp);
        TestType assigned;
        assigned = interp;
        _2D::AnyInterpolator<Real> any(interp);
        CHECK_OOB(copy(0, 20), Real(-1));
        CHECK_OOB(assigned(2, 50), Real(-1));
        CHECK_OOB(any(0, 0), Real(-1));
    }
}

TEMPLATE_PRODUCT_TEST_CASE("2D extrapolation reproduces a plane",
                           "[out-of-bounds][2D]",
                           (_2D::BilinearInterpolator, _2D::BicubicInterpolator,
                            _2D::ThinPlateSplineInterpolator),
                           (float, double)) {
    using Real = typename TestType::VectorType::Scalar;
    Plane<Real> g;
    TestType interp;
    interp.setData(g.x, g.y, g.z);
    interp.enableExtrapolation();
    CHECK(interp.extrapolationEnabled());

    // In float, the thin plate spline's radial weights are not exactly 0 for
    // a plane (about 2e-3 error inside the data), and their terms grow with
    // distance, so it needs a looser tolerance.
    const double eps = std::is_same<Real, float>::value ? 2e-3 : 1e-4;

    // x only, y only, and the corner regions
    for (Real x : {Real(-2), Real(0), Real(2), Real(4), Real(10)})
        for (Real y : {Real(0), Real(5), Real(25), Real(45), Real(100)})
            CHECK(interp(x, y) == Approx(g.f(x, y)).epsilon(eps));

    TestType copy(interp);
    _2D::AnyInterpolator<Real> any(interp);
    CHECK(copy(10, 100) == Approx(g.f(10, 100)).epsilon(eps));
    CHECK(any(0, 0) == Approx(g.f(0, 0)).epsilon(eps));
}

TEMPLATE_TEST_CASE("2D NearestNeighbor extrapolation returns the edge value",
                   "[out-of-bounds][2D][nearest]", float, double) {
    Plane<TestType> g;
    _2D::NearestNeighborInterpolator<TestType> interp;
    interp.setData(g.x, g.y, g.z);
    interp.enableExtrapolation();
    CHECK(interp(0, 20) == approx(g.f(1, 20)));
    CHECK(interp(10, 20) == approx(g.f(3, 20)));
    CHECK(interp(2, 0) == approx(g.f(2, 10)));
    CHECK(interp(2, 100) == approx(g.f(2, 40)));
    CHECK(interp(-5, -5) == approx(g.f(1, 10)));
    CHECK(interp(10, 100) == approx(g.f(3, 40)));
    CHECK(interp(10, 24) == approx(g.f(3, 20)));
}
