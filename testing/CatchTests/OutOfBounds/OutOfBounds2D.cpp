#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <functional>
#include <libInterpolate/AnyInterpolator.hpp>
#include <libInterpolate/Interpolators/_2D/BicubicInterpolator.hpp>
#include <libInterpolate/Interpolators/_2D/BilinearInterpolator.hpp>
#include <libInterpolate/Interpolators/_2D/LinearDelaunayTriangleInterpolator.hpp>
#include <libInterpolate/Interpolators/_2D/NearestNeighborInterpolator.hpp>
#include <libInterpolate/Interpolators/_2D/ThinPlateSplineInterpolator.hpp>
#include <vector>

#include "OutOfBoundsTestUtils.hpp"
using namespace Catch;

namespace {
// Grid: x in {1,2,3}, y in {10,20,30,40}, z = x + 2y + 100.
// z is never 0 in the domain, so a 0 result means out of bounds.
// Every method here reproduces z exactly at the knots. Only knots are used
// as expected values, so NearestNeighbor shares the same assertions.
template <typename Real>
struct Grid {
    std::vector<Real> x, y, z;
    Real xlo = 1, xhi = 3, ylo = 10, yhi = 40;
    static Real f(Real x, Real y) { return x + 2 * y + 100; }
    Grid() {
        for (Real a : {Real(1), Real(2), Real(3)})
            for (Real b : {Real(10), Real(20), Real(30), Real(40)}) {
                x.push_back(a);
                y.push_back(b);
                z.push_back(f(a, b));
            }
    }
};
}  // namespace

TEMPLATE_PRODUCT_TEST_CASE("2D out-of-bounds queries return 0 by default",
                           "[out-of-bounds][2D]",
                           (_2D::BilinearInterpolator, _2D::BicubicInterpolator,
                            _2D::NearestNeighborInterpolator,
                            _2D::ThinPlateSplineInterpolator,
                            _2D::LinearDelaunayTriangleInterpolator),
                           (float, double)) {
    using Real = typename TestType::VectorType::Scalar;
    using L = std::numeric_limits<Real>;
    Grid<Real> g;
    const Real oob = oob::default_value<Real>();
    const Real xm = 2, ym = 20;  // an interior knot
    TestType interp;
    interp.setData(g.x, g.y, g.z);

    SECTION("only x out of bounds") {
        CHECK_OOB(interp(g.xlo - 1, ym), oob);
        CHECK_OOB(interp(g.xhi + 1, ym), oob);
        CHECK_OOB(interp(g.xlo - 1, g.ylo), oob);  // y exactly on boundary
        CHECK_OOB(interp(g.xhi + 1, g.yhi), oob);
    }
    SECTION("only y out of bounds") {
        CHECK_OOB(interp(xm, g.ylo - 1), oob);
        CHECK_OOB(interp(xm, g.yhi + 1), oob);
        CHECK_OOB(interp(g.xlo, g.ylo - 1), oob);  // x exactly on boundary
        CHECK_OOB(interp(g.xhi, g.yhi + 1), oob);
    }
    SECTION("both out of bounds (corner regions)") {
        CHECK_OOB(interp(g.xlo - 1, g.ylo - 1), oob);
        CHECK_OOB(interp(g.xlo - 1, g.yhi + 1), oob);
        CHECK_OOB(interp(g.xhi + 1, g.ylo - 1), oob);
        CHECK_OOB(interp(g.xhi + 1, g.yhi + 1), oob);
    }
    SECTION("one ulp outside each edge") {
        CHECK_OOB(interp(oob::next_below(g.xlo), ym), oob);
        CHECK_OOB(interp(oob::next_above(g.xhi), ym), oob);
        CHECK_OOB(interp(xm, oob::next_below(g.ylo)), oob);
        CHECK_OOB(interp(xm, oob::next_above(g.yhi)), oob);
    }
    SECTION("far away and infinite") {
        CHECK_OOB(interp(Real(1e30), Real(1e30)), oob);
        CHECK_OOB(interp(Real(-1e30), Real(-1e30)), oob);
        CHECK_OOB(interp(Real(1e30), ym), oob);
        CHECK_OOB(interp(xm, Real(-1e30)), oob);
        CHECK_OOB(interp(L::infinity(), ym), oob);
        CHECK_OOB(interp(-L::infinity(), ym), oob);
        CHECK_OOB(interp(xm, L::infinity()), oob);
        CHECK_OOB(interp(xm, -L::infinity()), oob);
    }
    SECTION("boundary knots (edges and corners) are inside the domain") {
        for (Real yy : {g.ylo, ym, g.yhi}) {
            CHECK(interp(g.xlo, yy) == Approx(g.f(g.xlo, yy)));
            CHECK(interp(g.xhi, yy) == Approx(g.f(g.xhi, yy)));
        }
        CHECK(interp(xm, g.ylo) == Approx(g.f(xm, g.ylo)));
        CHECK(interp(xm, g.yhi) == Approx(g.f(xm, g.yhi)));
    }
    SECTION("one ulp inside each edge is inside the domain") {
        CHECK(interp(oob::next_above(g.xlo), ym) == Approx(g.f(g.xlo, ym)));
        CHECK(interp(oob::next_below(g.xhi), ym) == Approx(g.f(g.xhi, ym)));
        CHECK(interp(xm, oob::next_above(g.ylo)) == Approx(g.f(xm, g.ylo)));
        CHECK(interp(xm, oob::next_below(g.yhi)) == Approx(g.f(xm, g.yhi)));
    }
    SECTION("through std::function and _2D::AnyInterpolator") {
        std::function<Real(Real, Real)> fn = interp;
        CHECK_OOB(fn(g.xlo - 1, ym), oob);
        CHECK_OOB(fn(xm, g.yhi + 1), oob);
        _2D::AnyInterpolator<Real> any = TestType();
        any.setData(static_cast<int>(g.x.size()), g.x.data(), g.y.data(),
                    g.z.data());
        CHECK_OOB(any(g.xhi + 1, g.yhi + 1), oob);
        CHECK(any(g.xhi, g.yhi) == Approx(g.f(g.xhi, g.yhi)));
    }
    SECTION("copies keep the out-of-bounds behavior") {
        TestType copy(interp);
        CHECK_OOB(copy(g.xlo - 1, ym), oob);
        CHECK(copy(g.xlo, ym) == Approx(g.f(g.xlo, ym)));
    }
}

// NaN results are unspecified (see README) and differ between methods. These
// pin each method's current result so a change is noticed. Update them if the
// behavior changes on purpose. NearestNeighbor is not pinned (see below).
TEMPLATE_PRODUCT_TEST_CASE("2D NaN query propagates NaN",
                           "[out-of-bounds][2D][nan]",
                           (_2D::BilinearInterpolator, _2D::BicubicInterpolator,
                            _2D::ThinPlateSplineInterpolator),
                           (float, double)) {
    using Real = typename TestType::VectorType::Scalar;
    Grid<Real> g;
    TestType interp;
    interp.setData(g.x, g.y, g.z);
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    CHECK_THAT(interp(nan, 20), Matchers::IsNaN());
    CHECK_THAT(interp(2, nan), Matchers::IsNaN());
    CHECK_THAT(interp(nan, nan), Matchers::IsNaN());
}

TEMPLATE_TEST_CASE("2D LinearDelaunay NaN query returns 0",
                   "[out-of-bounds][2D][nan]", float, double) {
    Grid<TestType> g;
    _2D::LinearDelaunayTriangleInterpolator<TestType> interp;
    interp.setData(g.x, g.y, g.z);
    const TestType nan = std::numeric_limits<TestType>::quiet_NaN();
    CHECK_OOB(interp(nan, 20), TestType(0));
    CHECK_OOB(interp(2, nan), TestType(0));
}

// Scattered data. The two scattered-data methods define "the domain"
// differently:
//   ThinPlateSpline: the axis-aligned bounding box of the data (min/max of x,y)
//   LinearDelaunay:  the convex hull of the data (edges and vertices included)
TEMPLATE_TEST_CASE("2D scattered data: domain definition",
                   "[out-of-bounds][2D][scattered]", float, double) {
    using Real = TestType;
    auto f = [](Real x, Real y) { return x + 2 * y + 100; };
    // Convex hull is the triangle (0,0),(4,0),(0,4). Bounding box is [0,4]^2.
    std::vector<Real> x{0, 4, 0, 1, 2, 1}, y{0, 0, 4, 1, 1, 2}, z;
    for (size_t i = 0; i < x.size(); ++i) z.push_back(f(x[i], y[i]));

    SECTION("ThinPlateSpline uses the bounding box") {
        _2D::ThinPlateSplineInterpolator<Real> interp;
        interp.setData(x, y, z);
        CHECK(interp(3, 3) == Approx(f(3, 3)));  // outside hull, inside box
        CHECK(interp(4, 4) == Approx(f(4, 4)));  // box corner, not a data point
        CHECK_OOB(interp(oob::next_below(Real(0)), 1), Real(0));
        CHECK_OOB(interp(5, 1), Real(0));
        CHECK_OOB(interp(1, 5), Real(0));
    }
    SECTION("LinearDelaunay uses the convex hull") {
        _2D::LinearDelaunayTriangleInterpolator<Real> interp;
        interp.setData(x, y, z);
        CHECK(interp(0, 0) == Approx(f(0, 0)));  // hull vertex
        CHECK(interp(2, 2) == Approx(f(2, 2)));  // on hull edge
        CHECK_OOB(interp(3, 3), Real(0));        // in box, outside hull
        CHECK_OOB(interp(4, 4), Real(0));
        CHECK_OOB(interp(oob::next_below(Real(0)), 1), Real(0));
        CHECK_OOB(interp(5, 1), Real(0));
    }
}

// NearestNeighbor used to read one element before the start of the X/Y
// arrays for queries exactly on the lower x or y edge, or NaN. Run under
// AddressSanitizer to catch a regression. What a NaN query should return is
// not decided yet, so it is not checked here.
TEST_CASE("NearestNeighbor lower-edge queries",
          "[out-of-bounds][2D][nearest]") {
    Grid<double> g;
    _2D::NearestNeighborInterpolator<double> interp;
    interp.setData(g.x, g.y, g.z);
    CHECK(interp(g.xlo, 20) == Approx(g.f(g.xlo, 20)));
    CHECK(interp(2, g.ylo) == Approx(g.f(2, g.ylo)));
    (void)interp(std::nan(""), 20);  // memory safety only, value not checked
}
