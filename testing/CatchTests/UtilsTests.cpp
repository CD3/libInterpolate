#include "catch.hpp"


#include <vector>
#include <libInterpolate/Interpolators/_1D/InterpolatorBase.hpp>
#include <libInterpolate/Interpolators/_1D/CubicSplineInterpolator.hpp>
#include <type_traits>
TEST_CASE("RealTypeOf tests")
{
  CHECK( std::is_same<_1D::RealTypeOf<_1D::CubicSplineInterpolator<double>>::type, double>::value );
  CHECK( std::is_same<_1D::RealTypeOf<_1D::CubicSplineInterpolator<float>>::type, float>::value );
  CHECK( std::is_same<_1D::RealTypeOf<_1D::CubicSplineInterpolator<int>>::type, int>::value );
}
