#include "catch.hpp"

#include <libInterpolate/Interpolators/_1D/MonotonicInterpolator.hpp>


TEST_CASE( "MonotonicInterpolator<float> does not compile", "[bugs]" ) {

  std::vector<float> x(10),y(10);

  _1D::MonotonicInterpolator<float> interp;

  interp.setData(x,y);




}

