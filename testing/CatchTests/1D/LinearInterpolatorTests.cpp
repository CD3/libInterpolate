#include "catch.hpp"
#include "fakeit.hpp"

#include "Interpolators/_1D/LinearInterpolator.hpp"


TEST_CASE( "LinearInterpolator Tests", "[spline]" ) {

  _1D::LinearInterpolator<double> interp;

  std::vector<double> x,y;

  x.push_back(0); y.push_back(-1);
  x.push_back(5); y.push_back(9);  // slope of 2
  x.push_back(10); y.push_back(14);  // slope of 1 

  interp.setData(x,y);

  SECTION("Interpolation")
  {
    //REQUIRE( interp( -1 ) == Approx( -3 ) );
    REQUIRE( interp(  0 ) == Approx( -1 ) );
    REQUIRE( interp(  1 ) == Approx(  1 ) );
    REQUIRE( interp(  2 ) == Approx(  3 ) );
    REQUIRE( interp(  3 ) == Approx(  5 ) );
    REQUIRE( interp(  4 ) == Approx(  7 ) );
    REQUIRE( interp(  5 ) == Approx(  9 ) );
    REQUIRE( interp(  6 ) == Approx( 10 ) );
    REQUIRE( interp(  7 ) == Approx( 11 ) );
    REQUIRE( interp(  8 ) == Approx( 12 ) );
    REQUIRE( interp(  9 ) == Approx( 13 ) );
    REQUIRE( interp( 10 ) == Approx( 14 ) );
    //REQUIRE( interp( 11 ) == Approx(  0 ) );
  }

  SECTION("Multi-Interpolation")
  {
    std::vector<double> xvals = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    std::vector<double> yvals(9);

    interp( 9, xvals.data(), yvals.data() );

    REQUIRE( yvals[0] == Approx( 1) );
    REQUIRE( yvals[1] == Approx( 3) );
    REQUIRE( yvals[2] == Approx( 5) );
    REQUIRE( yvals[3] == Approx( 7) );
    REQUIRE( yvals[4] == Approx( 9) );
    REQUIRE( yvals[5] == Approx(10) );
    REQUIRE( yvals[6] == Approx(11) );
    REQUIRE( yvals[7] == Approx(12) );
    REQUIRE( yvals[8] == Approx(13) );

  }


  SECTION("Derivative")
  {
    //REQUIRE( interp.derivative( -1 ) == Approx( 0 ) );
    REQUIRE( interp.derivative(  0 ) == Approx( 2 ) );
    REQUIRE( interp.derivative(  1 ) == Approx( 2 ) );
    REQUIRE( interp.derivative(  2 ) == Approx( 2 ) );
    REQUIRE( interp.derivative(  3 ) == Approx( 2 ) );
    REQUIRE( interp.derivative(  4 ) == Approx( 2 ) );
    REQUIRE( interp.derivative(  5 ) == Approx( 2 ) );
    REQUIRE( interp.derivative(  6 ) == Approx( 1 ) );
    REQUIRE( interp.derivative(  7 ) == Approx( 1 ) );
    REQUIRE( interp.derivative(  8 ) == Approx( 1 ) );
    REQUIRE( interp.derivative(  9 ) == Approx( 1 ) );
    REQUIRE( interp.derivative( 10 ) == Approx( 1 ) );
    //REQUIRE( interp.derivative( 11 ) == Approx( 0 ) );
  }

}

