#include "catch.hpp"
#include "fakeit.hpp"

#include "Interpolators/_2D/InterpolatorBase.hpp"


namespace _2D {

class TestInterp : public InterpolatorBase<double>
{
  public:
    virtual double operator()( double x, double y )
    {
      return x + 2*y + 10;
    }

    VectorType getX() { return *(this->xv); }
    VectorType getY() { return *(this->yv); }
    VectorType getZ() { return *(this->zv); }

};

}


TEST_CASE( "2D InterpolatorBase Setup Tests", "[plumbing]" ) {

  _2D::TestInterp interp;

  // make sure interpolator works the way we expect
  REQUIRE( interp(1,1) == Approx(13) );
  REQUIRE( interp(10,20) == Approx(60) );

  int N = 10;

  SECTION("Eigen Vector Initialization")
  {
    _2D::InterpolatorBase<double>::VectorType xx(N), yy(N), zz(N);

    for( int i = 0; i < N; i++)
    {
      xx(i) = 0.1*i;
      yy(i) = 0.2*i;
      zz(i) = xx(i) + yy(i);
    }

    SECTION("Deep Copy")
    {
      interp.setData( xx, yy, zz );

      // clear the original data to make sure deep copy worked
      for(int i = 0; i < N; i++)
        xx(i) = yy(i) = 0;

      auto x = interp.getX();
      auto y = interp.getY();
      auto z = interp.getZ();

      for( int i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.1*i ) );
        REQUIRE( y(i) == Approx( 0.2*i ) );
        REQUIRE( z(i) == Approx( x(i)+y(i) ) );
      }
    }

    SECTION("Shallow Copy")
    {
      interp.setData( xx, yy, zz, false );

      // clear the original data to check that shallow copy worked
      for(int i = 0; i < N; i++)
        xx(i) = yy(i) = zz(i) = 0;

      auto x = interp.getX();
      auto y = interp.getY();
      auto z = interp.getZ();

      for( int i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.0 ) );
        REQUIRE( y(i) == Approx( 0.0 ) );
        REQUIRE( z(i) == Approx( 0.0 ) );
      }
    }

  }

  SECTION("Std Vector Initialization")
  {
    std::vector<double> xx(N), yy(N), zz(N);

    for( int i = 0; i < N; i++)
    {
      xx[i] = 0.1*i;
      yy[i] = 0.2*i;
      zz[i] = xx[i] + yy[i];
    }

    SECTION("Deep Copy")
    {
      interp.setData( xx, yy, zz );

      // clear the original data to make sure deep copy worked
      for(int i = 0; i < N; i++)
        xx[i] = yy[i] = zz[i] =  0;

      auto x = interp.getX();
      auto y = interp.getY();
      auto z = interp.getZ();

      for( int i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.1*i ) );
        REQUIRE( y(i) == Approx( 0.2*i ) );
        REQUIRE( z(i) == Approx( x(i) + y(i) ) );
      }
    }

    SECTION("Shallow Copy")
    {
      interp.setData( xx, yy, zz, false );

      // clear the original data to check that shallow copy worked
      for(int i = 0; i < N; i++)
        xx[i] = yy[i] = zz[i] = 0;

      auto x = interp.getX();
      auto y = interp.getY();
      auto z = interp.getZ();

      for( int i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.0 ) );
        REQUIRE( y(i) == Approx( 0.0 ) );
        REQUIRE( z(i) == Approx( 0.0 ) );
      }
    }

  }



  SECTION("Raw Pointer Initialization")
  {
    double *xx, *yy, *zz;
    xx = new double[N];
    yy = new double[N];
    zz = new double[N];


    for( int i = 0; i < N; i++)
    {
      xx[i] = 0.1*i;
      yy[i] = 0.2*i;
      zz[i] = xx[i] + yy[i];
    }

    SECTION("Deep Copy")
    {
      interp.setData( N, xx, yy, zz );

      // clear the original data to make sure deep copy worked
      for(int i = 0; i < N; i++)
        xx[i] = yy[i] = zz[i] = 0;

      auto x = interp.getX();
      auto y = interp.getY();
      auto z = interp.getZ();

      for( int i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.1*i ) );
        REQUIRE( y(i) == Approx( 0.2*i ) );
        REQUIRE( z(i) == Approx( x(i)+y(i) ) );
      }
    }

    SECTION("Shallow Copy")
    {
      interp.setData( N, xx, yy, zz, false );

      // clear the original data to check that shallow copy worked
      for(int i = 0; i < N; i++)
      for(int i = 0; i < N; i++)
        xx[i] = yy[i] = zz[i] = 0;

      auto x = interp.getX();
      auto y = interp.getY();
      auto z = interp.getZ();

      for( int i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.0 ) );
        REQUIRE( y(i) == Approx( 0.0 ) );
        REQUIRE( z(i) == Approx( 0.0 ) );
      }
    }

    delete[] xx;
    delete[] yy;
    delete[] zz;
  }


}


TEST_CASE( "2D InterpolatorBase Derivative/Integration Tests", "[plumbing]" ) {

  _2D::TestInterp interp;

  REQUIRE( interp(1,2) == Approx(15) );
  REQUIRE( interp(10,20) == Approx(60) );

  int N = 10;
  _2D::InterpolatorBase<double>::VectorType xx(N), yy(N), zz(N);

  for( int i = 0; i < N; i++)
  {
    xx(i) = 0.1*i;
    yy(i) = 0.2*i;
    zz(i) = 0; // this won't be used
  }
  interp.setData(xx,yy,zz);

  REQUIRE( interp.gradient(0.0,0.0)(0) == Approx(1) );
  REQUIRE( interp.gradient(0.0,0.0)(1) == Approx(2) );
  REQUIRE( interp.gradient(0.1,0.2)(0) == Approx(1) );
  REQUIRE( interp.gradient(0.1,0.2)(1) == Approx(2) );
  REQUIRE( interp.gradient(0.2,0.2)(0) == Approx(1) );
  REQUIRE( interp.gradient(0.2,0.2)(1) == Approx(2) );

  auto F = [](double x, double y){ return x*x*y/2 + x*y*y + 10*x*y; };

  REQUIRE( interp.integral(0.0,0.1,0.0,0.1) == Approx( (F(0.1,0.1) - F(0.0,0.1)) - (F(0.1,0.0) - F(0.1,0.0)) ) );
  REQUIRE( interp.integral(0.0,0.1,0.2,0.3) == Approx( (F(0.1,0.3) - F(0.0,0.3)) - (F(0.1,0.2) - F(0.0,0.2)) ) );
  REQUIRE( interp.integral(1.0,1.1,2.0,2.1) == Approx( (F(1.1,2.1) - F(1.0,2.1)) - (F(1.1,2.0) - F(1.0,2.0)) ) );

  REQUIRE( interp.integral(1.1,1.0,2.0,2.1) == Approx(-(F(1.1,2.1) - F(1.0,2.1)) + (F(1.1,2.0) - F(1.0,2.0)) ) );
  REQUIRE( interp.integral(1.0,1.1,2.1,2.0) == Approx(-(F(1.1,2.1) - F(1.0,2.1)) + (F(1.1,2.0) - F(1.0,2.0)) ) );
  REQUIRE( interp.integral(1.1,1.0,2.1,2.0) == Approx( (F(1.1,2.1) - F(1.0,2.1)) - (F(1.1,2.0) - F(1.0,2.0)) ) );


}
