#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#include <libInterpolate/Interpolators/_2D/BicubicInterpolator.hpp>
#include <libInterpolate/Interpolators/_2D/BilinearInterpolator.hpp>
#include <libInterpolate/Interpolators/_2D/LinearDelaunayTriangleInterpolator.hpp>
#include <libInterpolate/Interpolators/_2D/ThinPlateSplineInterpolator.hpp>

TEST_CASE("2D - Default construction with setData", "[construction]")
{
}

TEMPLATE_TEST_CASE("2D - Construction, Assignment, etc.","", _2D::BilinearInterpolator<double>, _2D::BicubicInterpolator<double>, _2D::ThinPlateSplineInterpolator<double>, _2D::LinearDelaunayTriangleInterpolator<double>)
{
  int    nx, ny;
  double xmin, xmax, dx;
  double ymin, ymax, dy;

  nx = 10;
  ny = 5;

  xmin = -1;
  xmax = 8;

  ymin = -1;
  ymax = 3;

  dx = (xmax - xmin) / (nx - 1);
  dy = (ymax - ymin) / (ny - 1);

  std::vector<double> xx(nx * ny), yy(nx * ny), zz(nx * ny);

  auto f = [](double x, double y) { return x * y + 2 * x + 3 * y; };

  for(int i = 0; i < nx * ny; i++) {
    // gnuplot format is essentially row-major
    xx[i] = xmin + dx * (i / ny);
    yy[i] = ymin + dy * (i % ny);
    zz[i] = f(xx[i], yy[i]);
  }

  SECTION("Copy construct")
  {
    TestType interp;
    interp.setData(xx, yy, zz);
    TestType interp2(interp);

    REQUIRE(interp2(0, 0) == Approx(interp(0, 0)));
    REQUIRE(interp2(1, 2) == Approx(interp(1, 2)));
    REQUIRE(interp2(2, 1) == Approx(interp(2, 1)));
    REQUIRE(interp2(2, -1) == Approx(interp(2, -1)));
    REQUIRE(interp2(8, 3) == Approx(interp(8, 3)));

    REQUIRE(interp2(-2, -1) == Approx(0).scale(1));
    REQUIRE(interp2(10, 3) == Approx(0).scale(1));
  }
}
