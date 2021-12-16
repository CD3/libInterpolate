#include "catch.hpp"

#include <libInterpolate/Interpolators/_2D/TriangleInterpolatorBase.hpp>


namespace _2D {

class TestTriangleInterp : public TriangleInterpolatorBase<TestTriangleInterp>
{
  protected:
    double c = 0;
  public:
    virtual double operator()( double x, double y ) const
    {
      return x + 2*y + this->c;
    }

    VectorType getX() { return *(this->xView); }
    VectorType getY() { return *(this->yView); }
    VectorType getZ() { return *(this->zView); }

    void setupInterpolator()
    {
      this->c = 10;
    }

};

}


TEST_CASE( "2D TriangleInterpolatorBase Setup Tests", "[plumbing]" ) {

  _2D::TestTriangleInterp interp;

  // make sure interpolator works the way we expect
  REQUIRE( interp(1,1) == Approx(3) );
  REQUIRE( interp(10,20) == Approx(50) );

  size_t N = 10;

  SECTION("Eigen Vector Initialization")
  {
    _2D::InterpolatorBase<double>::VectorType xx(N), yy(N), zz(N);

    for( size_t i = 0; i < N; i++)
    {
      xx(i) = 0.1*i;
      yy(i) = 0.2*i;
      zz(i) = xx(i) + yy(i);
    }

    SECTION("Deep Copy")
    {
      interp.setData( xx, yy, zz );

      REQUIRE( interp(1,1) == Approx(13) );
      REQUIRE( interp(10,20) == Approx(60) );

      // clear the original data to make sure deep copy worked
      for(size_t i = 0; i < N; i++)
        xx(i) = yy(i) = 0;

      auto x = interp.getX();
      auto y = interp.getY();
      auto z = interp.getZ();

      for( size_t i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.1*i ) );
        REQUIRE( y(i) == Approx( 0.2*i ) );
        REQUIRE( z(i) == Approx( x(i)+y(i) ) );
      }
    }

    //SECTION("Shallow Copy")
    //{
      //interp.setData( xx, yy, zz, false );

      //// clear the original data to check that shallow copy worked
      //for(int i = 0; i < N; i++)
        //xx(i) = yy(i) = zz(i) = 0;

      //auto x = interp.getX();
      //auto y = interp.getY();
      //auto z = interp.getZ();

      //for( int i = 0; i < N; i++)
      //{
        //REQUIRE( x(i) == Approx( 0.0 ) );
        //REQUIRE( y(i) == Approx( 0.0 ) );
        //REQUIRE( z(i) == Approx( 0.0 ) );
      //}
    //}

  }

  SECTION("Std Vector Initialization")
  {
    std::vector<double> xx(N), yy(N), zz(N);

    for( size_t i = 0; i < N; i++)
    {
      xx[i] = 0.1*i;
      yy[i] = 0.2*i;
      zz[i] = xx[i] + yy[i];
    }

    SECTION("Deep Copy")
    {
      interp.setData( xx, yy, zz );

      // clear the original data to make sure deep copy worked
      for(size_t i = 0; i < N; i++)
        xx[i] = yy[i] = zz[i] =  0;

      auto x = interp.getX();
      auto y = interp.getY();
      auto z = interp.getZ();

      for( size_t i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.1*i ) );
        REQUIRE( y(i) == Approx( 0.2*i ) );
        REQUIRE( z(i) == Approx( x(i) + y(i) ) );
      }
    }

    //SECTION("Shallow Copy")
    //{
      //interp.setData( xx, yy, zz, false );

      //// clear the original data to check that shallow copy worked
      //for(int i = 0; i < N; i++)
        //xx[i] = yy[i] = zz[i] = 0;

      //auto x = interp.getX();
      //auto y = interp.getY();
      //auto z = interp.getZ();

      //for( int i = 0; i < N; i++)
      //{
        //REQUIRE( x(i) == Approx( 0.0 ) );
        //REQUIRE( y(i) == Approx( 0.0 ) );
        //REQUIRE( z(i) == Approx( 0.0 ) );
      //}
    //}

  }

  SECTION("Raw Pointer Initialization")
  {
    double *xx, *yy, *zz;
    xx = new double[N];
    yy = new double[N];
    zz = new double[N];


    for( size_t i = 0; i < N; i++)
    {
      xx[i] = 0.1*i;
      yy[i] = 0.2*i;
      zz[i] = xx[i] + yy[i];
    }

    SECTION("Deep Copy")
    {
      interp.setData( N, xx, yy, zz );

      // clear the original data to make sure deep copy worked
      for(size_t i = 0; i < N; i++)
        xx[i] = yy[i] = zz[i] = 0;

      auto x = interp.getX();
      auto y = interp.getY();
      auto z = interp.getZ();

      for( size_t i = 0; i < N; i++)
      {
        REQUIRE( x(i) == Approx( 0.1*i ) );
        REQUIRE( y(i) == Approx( 0.2*i ) );
        REQUIRE( z(i) == Approx( x(i)+y(i) ) );
      }
    }

    //SECTION("Shallow Copy")
    //{
      //interp.setData( N, xx, yy, zz, false );

      //// clear the original data to check that shallow copy worked
      //for(int i = 0; i < N; i++)
      //for(int i = 0; i < N; i++)
        //xx[i] = yy[i] = zz[i] = 0;

      //auto x = interp.getX();
      //auto y = interp.getY();
      //auto z = interp.getZ();

      //for( int i = 0; i < N; i++)
      //{
        //REQUIRE( x(i) == Approx( 0.0 ) );
        //REQUIRE( y(i) == Approx( 0.0 ) );
        //REQUIRE( z(i) == Approx( 0.0 ) );
      //}
    //}

    delete[] xx;
    delete[] yy;
    delete[] zz;
  }

  SECTION("Eigen Vector Initialization")
  {
    _2D::InterpolatorBase<double>::VectorType xx(N), yy(N), zz(N);

    for( size_t i = 0; i < N; i++)
    {
      xx(i) = 0.1*i;
      yy(i) = 0.2*i;
      zz(i) = xx(i) + yy(i);
    }

    interp.setData( xx, yy, zz );

    auto x = interp.getXData();
    auto y = interp.getYData();
    auto z = interp.getZData();

    REQUIRE( x.size() > 0 );
    CHECK( x[0] == Approx(0) );
    REQUIRE( x.size() == N );
    CHECK( x[N-1] == Approx( xx(N-1) ) );

    REQUIRE( y.size() > 0 );
    CHECK( y[0] == Approx(0) );
    REQUIRE( y.size() == N );
    CHECK( y[N-1] == Approx( yy(N-1) ) );

    REQUIRE( z.size() > 0 );
    CHECK( z[0] == Approx(0) );
    REQUIRE( z.size() == N );
    CHECK( z[N-1] == Approx( xx(N-1) + yy(N-1) ) );



  }

  SECTION("Grid-formatted Data")
  {
    std::vector<int> x(3), y(4), z(12);
    x[0] = 0; x[1] = 1; x[2] = 2;
    y[0] = 2; y[1] = 3; y[2] = 4; y[3] = 5;

    SECTION("Iterators")
    {
      interp.setData(x.begin(), x.end(), y.begin(), y.end(), z.begin(), z.end() );
      auto X = interp.getX();
      auto Y = interp.getY();
      auto Z = interp.getZ();

      CHECK( X.size() == 12 );
      CHECK( X[0] == 0 );
      CHECK( X[1] == 0 );
      CHECK( X[2] == 0 );
      CHECK( X[3] == 0 );

      CHECK( X[4] == 1 );
      CHECK( X[5] == 1 );
      CHECK( X[6] == 1 );
      CHECK( X[7] == 1 );

      CHECK( X[8] == 2 );
      CHECK( X[9] == 2 );
      CHECK( X[10] == 2 );
      CHECK( X[11] == 2 );

      CHECK( Y[0] == 2 );
      CHECK( Y[1] == 3 );
      CHECK( Y[2] == 4 );
      CHECK( Y[3] == 5 );
      CHECK( Y[4] == 2 );
      CHECK( Y[5] == 3 );
      CHECK( Y[6] == 4 );
      CHECK( Y[7] == 5 );
      CHECK( Y[8] == 2 );
      CHECK( Y[9] == 3 );
      CHECK( Y[10] == 4 );
      CHECK( Y[11] == 5 );
    }
    SECTION("std::vector")
    {
      interp.setData(x, y, z);
      auto X = interp.getX();
      auto Y = interp.getY();
      auto Z = interp.getZ();

      CHECK( X.size() == 12 );
      CHECK( X[0] == 0 );
      CHECK( X[1] == 0 );
      CHECK( X[2] == 0 );
      CHECK( X[3] == 0 );

      CHECK( X[4] == 1 );
      CHECK( X[5] == 1 );
      CHECK( X[6] == 1 );
      CHECK( X[7] == 1 );

      CHECK( X[8] == 2 );
      CHECK( X[9] == 2 );
      CHECK( X[10] == 2 );
      CHECK( X[11] == 2 );

      CHECK( Y[0] == 2 );
      CHECK( Y[1] == 3 );
      CHECK( Y[2] == 4 );
      CHECK( Y[3] == 5 );
      CHECK( Y[4] == 2 );
      CHECK( Y[5] == 3 );
      CHECK( Y[6] == 4 );
      CHECK( Y[7] == 5 );
      CHECK( Y[8] == 2 );
      CHECK( Y[9] == 3 );
      CHECK( Y[10] == 4 );
      CHECK( Y[11] == 5 );
    }



  }

}


