#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <libInterpolate/Interpolators/_1D/MonotonicInterpolator.hpp>
using namespace Catch;

TEST_CASE("MonotonicInterpolator<float> does not compile", "[bugs]") {
    std::vector<float> x(10), y(10);

    _1D::MonotonicInterpolator<float> interp;

    interp.setData(x, y);
}
