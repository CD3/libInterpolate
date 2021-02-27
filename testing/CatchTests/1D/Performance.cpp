#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#include <libInterpolate/Interpolate.hpp>

TEMPLATE_TEST_CASE("1D Interplation Benchmarks", "[.][benchmarks]", _1D::LinearInterpolator<double>, _1D::MonotonicInterpolator<double>, _1D::CubicSplineInterpolator<double>)
{
  int                 N = 200000;
  std::vector<double> x(N), y(N);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i] = i * 0.01;
    y[i] = x[i] * x[i];
  }

  TestType interp;
  interp.setData(x, y);

  BENCHMARK("Find left element near 'front'")
  {
    return interp.get_index_to_left_of(1);
  };
  BENCHMARK("Find left element near 'middle'")
  {
    return interp.get_index_to_left_of(1000);
  };
  BENCHMARK("Find left element near 'end'")
  {
    return interp.get_index_to_left_of(1999);
  };
  BENCHMARK("Find left element near the 1/4 mark") 
  {
    return interp.get_index_to_left_of(498);
  };
  BENCHMARK("Interpolate near 'front'") { return interp(1); };
  BENCHMARK("Interpolate near 'middle'") { return interp(1000); };
  BENCHMARK("Interpolate near 'end'") { return interp(1999); };
  BENCHMARK("Interpolate near the 1/4 mark") { return interp(498); };
}

