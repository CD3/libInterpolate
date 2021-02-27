#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#include <libInterpolate/Interpolate.hpp>

TEMPLATE_TEST_CASE("2D Interplation Benchmarks", "[.][benchmarks]", _2D::BilinearInterpolator<double>, _2D::BicubicInterpolator<double>)
{
  int                 Nx = 2000;
  int                 Ny = 2000;
  std::vector<double> x(Nx*Ny), y(Nx*Ny), z(Nx*Ny);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i] = (i/Ny) * 0.01;
    y[i] = (i%Ny) * 0.01;
    z[i] = x[i]*y[i];
  }

  TestType interp;
  interp.setData(x,y,z);

  BENCHMARK("Interpolate near a corner") { return interp(0.02,0.02); };
  BENCHMARK("Interpolate near the middle") { return interp(10,10); };
}

