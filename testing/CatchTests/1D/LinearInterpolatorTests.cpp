#include "catch.hpp"

#include <libInterpolate/Interpolators/_1D/LinearInterpolator.hpp>


TEST_CASE( "LinearInterpolator Tests", "[linear]" ) {

  SECTION("Double Precision")
  {
  _1D::LinearInterpolator<double> interp;

  std::vector<double> x,y;

  x.push_back(0); y.push_back(-1);
  x.push_back(5); y.push_back(9);  // slope of 2
  x.push_back(10); y.push_back(14);  // slope of 1 


  SECTION("Interpolation w/ Copied Data")
  {
    interp.setData(x,y);
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
    
    std::vector<double> xs = {0,1,2,3,4,5,6,7,8,9,10};
    auto ys = interp.batch(xs);
    REQUIRE( ys[  0 ] == Approx( -1 ) );
    REQUIRE( ys[  1 ] == Approx(  1 ) );
    REQUIRE( ys[  2 ] == Approx(  3 ) );
    REQUIRE( ys[  3 ] == Approx(  5 ) );
    REQUIRE( ys[  4 ] == Approx(  7 ) );
    REQUIRE( ys[  5 ] == Approx(  9 ) );
    REQUIRE( ys[  6 ] == Approx( 10 ) );
    REQUIRE( ys[  7 ] == Approx( 11 ) );
    REQUIRE( ys[  8 ] == Approx( 12 ) );
    REQUIRE( ys[  9 ] == Approx( 13 ) );
    REQUIRE( ys[ 10 ] == Approx( 14 ) );
  }

  SECTION("Interpolation w/ ReferencedData")
  {
    interp.setUnsafeDataReference(x.size(),x.data(),y.data());
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


  SECTION("From single-precision data")
  {
    std::vector<float> x,y;

    x.push_back(0); y.push_back(-1);
    x.push_back(5); y.push_back(9);  // slope of 2
    x.push_back(10); y.push_back(14);  // slope of 1 


    SECTION("Interpolation w/ Copied Data")
    {
      interp.setData(x,y);
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
  }
  }

  SECTION("Single Precision")
  {
  _1D::LinearInterpolator<float> interp;

  std::vector<float> x,y;

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
  }

}

