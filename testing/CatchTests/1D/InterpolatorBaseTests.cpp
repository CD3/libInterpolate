#include "catch.hpp"
#include "fakeit.hpp"

#include "Interpolators/_1D/InterpolatorBase.hpp"


namespace _1D {

class TestInterp : public InterpolatorBase<double>
{
  public:
    virtual double operator()( double x ) const
    {
      return x + 10;
    }

    VectorType getX() { return *(this->xv); }
    VectorType getY() { return *(this->yv); }

};

}


TEST_CASE( "InterpolatorBase Setup Tests", "[plumbing]" ) {

  _1D::TestInterp interp;

  // check that the interpolator works as expected
  REQUIRE( interp(1) == Approx(11) );
  REQUIRE( interp(10) == Approx(20) );

  int N = 10;

  SECTION("Eigen Vector Initialization")
  {
    _1D::InterpolatorBase<double>::VectorType xx(N), yy(N);

    for( int i = 0; i < N; i++)
    {
      xx(i) = 0.1*i;
      yy(i) = xx(i)*xx(i);
    }

    SECTION("Deep Copy")
    {
      interp.setData( xx, yy );

      // clear the original data to make sure deep copy worked
      for(int i = 0; i < N; i++)
        xx(i) = yy(i) = 0;

      auto x = interp.getX();
      auto y = interp.getY();

      for( int i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.1*i ) );
        REQUIRE( y(i) == Approx( x(i)*x(i) ) );
      }
    }

    SECTION("Shallow Copy")
    {
      interp.setData( xx, yy, false );

      // clear the original data to check that shallow copy worked
      for(int i = 0; i < N; i++)
        xx(i) = yy(i) = 0;

      auto x = interp.getX();
      auto y = interp.getY();

      for( int i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.0 ) );
        REQUIRE( y(i) == Approx( 0.0 ) );
      }
    }

  }

  SECTION("Std Vector Initialization")
  {
    std::vector<double> xx(N), yy(N);

    for( int i = 0; i < N; i++)
    {
      xx[i] = 0.1*i;
      yy[i] = xx[i]*xx[i];
    }

    SECTION("Deep Copy")
    {
      interp.setData( xx, yy );

      // clear the original data to make sure deep copy worked
      for(int i = 0; i < N; i++)
        xx[i] = yy[i] = 0;

      auto x = interp.getX();
      auto y = interp.getY();

      for( int i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.1*i ) );
        REQUIRE( y(i) == Approx( x(i)*x(i) ) );
      }
    }

    SECTION("Shallow Copy")
    {
      interp.setData( xx, yy, false );

      // clear the original data to check that shallow copy worked
      for(int i = 0; i < N; i++)
        xx[i] = yy[i] = 0;

      auto x = interp.getX();
      auto y = interp.getY();

      for( int i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.0 ) );
        REQUIRE( y(i) == Approx( 0.0 ) );
      }
    }

  }

  SECTION("Raw Pointer Initialization")
  {
    double *xx, *yy;
    xx = new double[N];
    yy = new double[N];


    for( int i = 0; i < N; i++)
    {
      xx[i] = 0.1*i;
      yy[i] = xx[i]*xx[i];
    }

    SECTION("Deep Copy")
    {
      interp.setData( N, xx, yy );

      // clear the original data to make sure deep copy worked
      for(int i = 0; i < N; i++)
        xx[i] = yy[i] = 0;

      auto x = interp.getX();
      auto y = interp.getY();

      for( int i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.1*i ) );
        REQUIRE( y(i) == Approx( x(i)*x(i) ) );
      }
    }

    SECTION("Shallow Copy")
    {
      interp.setData( N, xx, yy, false );

      // clear the original data to check that shallow copy worked
      for(int i = 0; i < N; i++)
      for(int i = 0; i < N; i++)
        xx[i] = yy[i] = 0;

      auto x = interp.getX();
      auto y = interp.getY();

      for( int i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.0 ) );
        REQUIRE( y(i) == Approx( 0.0 ) );
      }
    }

    delete[] xx;
    delete[] yy;
  }

  SECTION("Data Get Functions")
  {
    _1D::InterpolatorBase<double>::VectorType xx(N), yy(N);

    for( int i = 0; i < N; i++)
    {
      xx(i) = 0.1*i;
      yy(i) = xx(i)*xx(i);
    }

    interp.setData( xx, yy );

    auto x = interp.getXData();
    auto y = interp.getYData();

    REQUIRE( x.size() > 0 );
    CHECK( x[0] == Approx(0) );
    REQUIRE( x.size() == N );
    CHECK( x[N-1] == Approx( xx(N-1) ) );

    REQUIRE( y.size() > 0 );
    CHECK( y[0] == Approx(0) );
    REQUIRE( y.size() == N );
    CHECK( y[N-1] == Approx( pow(xx(N-1),2) ) );

  }

}


TEST_CASE( "InterpolatorBase Derivative/Integration Tests", "[plumbing]" ) {

  _1D::TestInterp interp;

  REQUIRE( interp(1) == Approx(11) );
  REQUIRE( interp(10) == Approx(20) );

  int N = 10;
  _1D::InterpolatorBase<double>::VectorType xx(N), yy(N);

  for( int i = 0; i < N; i++)
  {
    xx(i) = 0.1*i;
    yy(i) = 1;
  }
  interp.setData(xx,yy);

  REQUIRE( interp.derivative(0.0) == Approx(1) );
  REQUIRE( interp.derivative(0.1) == Approx(1) );
  REQUIRE( interp.derivative(0.2) == Approx(1) );

  REQUIRE( interp.integral(0.0,0.1) == Approx(10.1*10.1/2 - 10*10/2) );
  REQUIRE( interp.integral(0.1,0.2) == Approx(10.2*10.2/2 - 10.1*10.1/2) );
  REQUIRE( interp.integral(0.2,0.3) == Approx(10.3*10.3/2 - 10.2*10.2/2) );
  // switching limits should pick up a negative sign
  REQUIRE( interp.integral(0.1,0.0) == Approx(-10.1*10.1/2 + 10*10/2) );
  REQUIRE( interp.integral(0.2,0.1) == Approx(-10.2*10.2/2 + 10.1*10.1/2) );
  REQUIRE( interp.integral(0.3,0.2) == Approx(-10.3*10.3/2 + 10.2*10.2/2) );


}
