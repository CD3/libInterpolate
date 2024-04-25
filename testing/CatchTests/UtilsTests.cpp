#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
using namespace Catch;

#include <libInterpolate/Interpolators/_1D/CubicSplineInterpolator.hpp>
#include <libInterpolate/Interpolators/_1D/InterpolatorBase.hpp>
#include <type_traits>
#include <vector>
TEST_CASE("RealTypeOf tests") {
    CHECK(std::is_same<
          _1D::RealTypeOf<_1D::CubicSplineInterpolator<double>>::type,
          double>::value);
    CHECK(
        std::is_same<_1D::RealTypeOf<_1D::CubicSplineInterpolator<float>>::type,
                     float>::value);
    CHECK(std::is_same<_1D::RealTypeOf<_1D::CubicSplineInterpolator<int>>::type,
                       int>::value);
}
