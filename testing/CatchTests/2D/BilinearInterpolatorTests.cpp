#include "catch.hpp"
#include "fakeit.hpp"

#include "Interpolators/_2D/BilinearInterpolator.hpp"

namespace _2D {
  class TestBilinearInterp : BilinearInterpolator<double>
  {
    public:
      VectorType getX() { return *(this->X); }
      VectorType getY() { return *(this->Y); }
      MatrixType getZ() { return *(this->Z); }

      using BilinearInterpolator<double>::setData;
      using BilinearInterpolator<double>::operator();
      using BilinearInterpolator<double>::integral;
      using BilinearInterpolator<double>::gradient;
  };
}


TEST_CASE( "BilinearInterpolator Tests", "[bilinear]" ) {

  _2D::TestBilinearInterp interp;

  int nx, ny;
  double xmin, xmax, dx, x;
  double ymin, ymax, dy, y;
  double z;

  nx = 10;
  ny = 5;

  xmin = -1;
  xmax = 8;

  ymin = -1;
  ymax = 3;

  dx = (xmax - xmin)/(nx - 1);
  dy = (ymax - ymin)/(ny - 1);

  _2D::BilinearInterpolator<double>::VectorType xx(nx*ny), yy(nx*ny), zz(nx*ny);

  auto f  = [](double x, double y){return x*y + 2*x + 3*y;};

  for( int i = 0; i < nx*ny; i++)
  {
    // gnuplot format is essentially row-major
    xx(i) = xmin+dx*(i/ny);
    yy(i) = ymin+dy*(i%ny);
    zz(i) = f(xx(i),yy(i));
  }
  interp.setData( xx, yy, zz );

  SECTION("1D -> 2D Mappings")
  {
    auto X = interp.getX();
    auto Y = interp.getY();
    auto Z = interp.getZ();

    REQUIRE( X.size() == nx );
    REQUIRE( X(0) == Approx( xmin+0    ) );
    REQUIRE( X(2) == Approx( xmin+2*dx ) );
    REQUIRE( X(4) == Approx( xmin+4*dx ) );

    REQUIRE( Y.size() == ny );
    REQUIRE( Y(0) == Approx( ymin+0    ) );
    REQUIRE( Y(2) == Approx( ymin+2*dy ) );
    REQUIRE( Y(4) == Approx( ymin+4*dy ) );

    REQUIRE( Z.size() == nx*ny );
    REQUIRE( Z(0,0) == Approx( f(xmin+0*dx,ymin+0*dy) ) );
    REQUIRE( Z(0,2) == Approx( f(xmin+0*dx,ymin+2*dy) ) );
    REQUIRE( Z(2,0) == Approx( f(xmin+2*dx,ymin+0*dy) ) );
    REQUIRE( Z(2,2) == Approx( f(xmin+2*dx,ymin+2*dy) ) );
  }


  SECTION("Interpolation")
  {
    REQUIRE( interp(0,0)   == Approx(f(0,0)).epsilon(0.00002 )) ;
    REQUIRE( interp(1,2)   == Approx(f(1,2)).epsilon(0.00002 )) ;
    REQUIRE( interp(2,1)   == Approx(f(2,1)).epsilon(0.00002 )) ;
    REQUIRE( interp(2,-1)  == Approx(f(2,-1)).epsilon(0.00002 )) ;
    REQUIRE( interp(8,3)   == Approx(f(8,3)).epsilon(0.00002 )) ;

    REQUIRE( interp(-2,-1) == Approx(0).epsilon(0.00002 )) ;
    REQUIRE( interp(10,3)  == Approx(0).epsilon(0.00002 )) ;
  }


  SECTION("Gradient")
  {
    auto dfdx = [](double x, double y){return y + 2;};
    auto dfdy = [](double x, double y){return x + 3;};
    double x,y;
    x = xmin + 1.5*dx; y = ymin + 1.5*dy;
    REQUIRE( interp.gradient( x,y )(0) == Approx( dfdx(x,y) ).epsilon(0.0003));
    REQUIRE( interp.gradient( x,y )(1) == Approx( dfdy(x,y) ).epsilon(0.0003));
    x = xmin + 3*dx; y = ymin + 1*dy;                                               
    REQUIRE( interp.gradient( x,y )(0) == Approx( dfdx(x,y) ).epsilon(0.0003));
    REQUIRE( interp.gradient( x,y )(1) == Approx( dfdy(x,y) ).epsilon(0.0003));
    x = xmin + 1.5*dx; y = ymin + 3.25*dy;                                               
    REQUIRE( interp.gradient( x,y )(0) == Approx( dfdx(x,y) ).epsilon(0.0003));
    REQUIRE( interp.gradient( x,y )(1) == Approx( dfdy(x,y) ).epsilon(0.0003));

    x = -2; y = -1;
    REQUIRE( interp.gradient( x,y )(0) == Approx( 0 ).epsilon(0.0003));
    REQUIRE( interp.gradient( x,y )(1) == Approx( 0 ).epsilon(0.0003));
    x = 10; y = 3;
    REQUIRE( interp.gradient( x,y )(0) == Approx( 0 ).epsilon(0.0003));
    REQUIRE( interp.gradient( x,y )(1) == Approx( 0 ).epsilon(0.0003));
  }



  SECTION("Integral")
  {


    // \int f(x,y) dx dy = x^2 y^2 / 4 + x^2 y + 3 x y^2 / 2 = g(x,y) |_xa^xb | _ya^yb
    // 
    // = ( g(xb,yb) - g(xb,ya) ) - ( g(xa,yb) - g(xa,ya) ) = g(xb,yb) - g(xb,ya) - g(xa,yb) + g(xa,ya)
    //
    
    double sum, xa, xb, ya, yb;

    auto integral = [](double xa, double xb, double ya, double yb) {
      double sum = 0;
      double x,y;
      x = xb; y = yb;
      sum += x*x*y*y/4 + x*x*y + 3*x*y*y/2;
      x = xb; y = ya;
      sum -= x*x*y*y/4 + x*x*y + 3*x*y*y/2;
      x = xa; y = yb;
      sum -= x*x*y*y/4 + x*x*y + 3*x*y*y/2;
      x = xa; y = ya;
      sum += x*x*y*y/4 + x*x*y + 3*x*y*y/2;
      return sum;
    };


  // function we are interpolating is defined over
  // x -> [-1,8] (increments of 1)
  // y -> [-1,3] (increments of 1)

  // integral over domain containing only whole element
  REQUIRE( interp.integral( 1, 3, 0, 2 ) == integral( 1, 3, 0, 2 ) );

  // integral over domain containing some partial elements, but each corner is in a different element
  REQUIRE( interp.integral( 1.5, 3.5, 0.5, 2.5 ) == integral( 1.5, 3.5, 0.5, 2.5 ) );

  // integral over domain with the x limits contained in the same element
  REQUIRE( interp.integral( 1.5, 1.75, 0, 2 ) == integral( 1.5, 1.75, 0, 2 ) );

  // same thing, but offset y limits
  REQUIRE( interp.integral( 1.5, 1.75, 0.5, 2.5 ) == integral( 1.5, 1.75, 0.5, 2.5 ) );

  // integral over domain with the y limits contained in the same element
  REQUIRE( interp.integral( 1, 3, 0.5, 0.75 ) == integral( 1, 3, 0.5, 0.75 ) );

  // same thing, but offset x limits
  REQUIRE( interp.integral( 1.5, 3.5, 0.5, 0.75 ) == integral( 1.5, 3.5, 0.5, 0.75 ) );

  // integral over domain with both x and y limits in the same element
  REQUIRE( interp.integral( 1.5, 1.75, 0.5, 0.75 ) == integral( 1.5, 1.75, 0.5, 0.75 ) );



  }



}


