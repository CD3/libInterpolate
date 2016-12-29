#include "catch.hpp"
#include "fakeit.hpp"

#include "Interpolators/_1D/SplineInterpolator.hpp"


TEST_CASE( "SplineInterpolator Tests", "[spline]" ) {

  _1D::SplineInterpolator<double> interp;

  int N = 100;
  double xmin = 0, xmax = 2*M_PI;
  double dx = (xmax - xmin)/(N-1);

  _1D::SplineInterpolator<double>::VectorType xx(N), yy(N);

  for( int i = 0; i < N; i++)
  {
    xx(i) = dx*i;
    yy(i) = sin(xx(i));
  }
  interp.setData( xx, yy );


  SECTION("Interpolation")
  {
    REQUIRE( interp( dx/2          ) == Approx( sin( dx/2          ) ) );
    REQUIRE( interp( M_PI/2 - dx/2 ) == Approx( sin( M_PI/2 - dx/2 ) ) );
    REQUIRE( interp( M_PI/4 - dx/4 ) == Approx( sin( M_PI/4 - dx/4 ) ) );
  }


  SECTION("Derivative")
  {
    REQUIRE( interp.derivative( dx/2          ) == Approx( cos( dx/2          ) ) );
    REQUIRE( interp.derivative( M_PI/2 - dx/2 ) == Approx( cos( M_PI/2 - dx/2 ) ) );
    REQUIRE( interp.derivative( M_PI/4 - dx/4 ) == Approx( cos( M_PI/4 - dx/4 ) ) );
  }


  SECTION("Integral")
  {
    REQUIRE( interp.integral( 0, dx/2          ) == Approx( -cos( dx/2          ) + cos(0) ) );
    REQUIRE( interp.integral( 0, M_PI/2 - dx/2 ) == Approx( -cos( M_PI/2 - dx/2 ) + cos(0) ) );
    REQUIRE( interp.integral( 0, M_PI/4 - dx/4 ) == Approx( -cos( M_PI/4 - dx/4 ) + cos(0) ) );
  }


}

TEST_CASE("SplineInterpolator Edge Case Tests")
{
  _1D::SplineInterpolator<double> interp;

  std::vector<double> x, y;
  x.push_back( 1.0 ); y.push_back( 2.0 );
  x.push_back( 2.0 ); y.push_back( 8.0 );
  x.push_back( 3.0 ); y.push_back(16.0 );
  x.push_back( 4.0 ); y.push_back(64.0 );

  interp.setData(x,y);

  REQUIRE( interp(1.0) == Approx( 2.0));
  REQUIRE( interp(2.0) == Approx( 8.0));
  REQUIRE( interp(3.0) == Approx(16.0));
  REQUIRE( interp(4.0) == Approx(64.0));

}
