#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#include <libInterpolate/Interpolate.hpp>

TEMPLATE_TEST_CASE("1D Interpolation Benchmarks", "[.][benchmarks]", _1D::LinearInterpolator<double>, _1D::MonotonicInterpolator<double>, _1D::CubicSplineInterpolator<double>)
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

  BENCHMARK("setData()")
  {
    TestType interp2;
    interp2.setData(x,y);
  };
  BENCHMARK("setUnsafeDataReference()")
  {
    TestType interp2;
    interp2.setUnsafeDataReference(x.size(),x.data(),y.data());
  };
}



#ifdef HAVE_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

TEST_CASE("GSL Comparisons","[.][benchmarks]")
{
  int                 N = 200000;
  std::vector<double> x(N), y(N);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i] = i * (2000. - 0.)/(N-1);
    y[i] = x[i] * x[i];
  }

  SECTION("Cubic Spline Interpolation")
  {
  _1D::CubicSplineInterpolator<double> interp;
  interp.setData(x, y);

  BENCHMARK("libInterpolate at 'front'")
  {
    return interp(0.1);
  };

  BENCHMARK("libInterpolate at 'front'")
  {
    return interp(1000);
  };
  BENCHMARK("libInterpolate at 'front' and 'back' consecutive")
  {
    return interp(1) + interp(1999);
  };

  SECTION("GSL Standard")
  {
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, N);
    gsl_spline_init (spline, x.data(), y.data(), N);

    BENCHMARK("Interpolation at 'front'")
    {
      return gsl_spline_eval(spline, 0.1, NULL );
    };
    BENCHMARK("Interpolation at 'middle'")
    {
      return gsl_spline_eval(spline, 1000, NULL );
    };
    BENCHMARK("Interpolation at 'front' and 'back' together")
    {
      return gsl_spline_eval (spline, 1, NULL)+gsl_spline_eval(spline,1999, NULL);
    };


    gsl_spline_free(spline);
  }

  SECTION("GSL Accelerated")
  {
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, N);
    gsl_spline_init (spline, x.data(), y.data(), N);

    BENCHMARK("Interpolation at 'front'")
    {
      return gsl_spline_eval (spline, 0.1, acc);
    };
    BENCHMARK("Interpolation at 'middle'")
    {
      return gsl_spline_eval (spline, 1000, acc);
    };
    BENCHMARK("Interpolation at 'front' and 'back' together")
    {
      return gsl_spline_eval (spline, 1, acc)+gsl_spline_eval(spline,1999,acc);
    };


    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }

  }

}

#endif
