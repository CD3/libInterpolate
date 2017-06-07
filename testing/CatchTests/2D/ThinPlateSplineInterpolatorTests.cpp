#include "catch.hpp"
#include "fakeit.hpp"

#include <fstream>

#include "Interpolators/_2D/ThinPlateSplineInterpolator.hpp"

namespace _2D {
  class TestThinPlateSplineInterp : ThinPlateSplineInterpolator<double>
  {
    public:
      using ThinPlateSplineInterpolator<double>::setData;
      using ThinPlateSplineInterpolator<double>::operator();
      using ThinPlateSplineInterpolator<double>::integral;
      using ThinPlateSplineInterpolator<double>::gradient;
  };
}


TEST_CASE( "ThinPlateSplineInterpolator Tests - Monotonic Data", "[thin plate spline]" ) {

  _2D::TestThinPlateSplineInterp interp;

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

  _2D::ThinPlateSplineInterpolator<double>::VectorType xx(nx*ny), yy(nx*ny), zz(nx*ny);

  auto f  = [](double x, double y){return x*y + 2*x + 3*y;};

  for( int i = 0; i < nx*ny; i++)
  {
    // gnuplot format is essentially row-major
    xx(i) = xmin+dx*(i/ny);
    yy(i) = ymin+dy*(i%ny);
    zz(i) = f(xx(i),yy(i));
  }
  interp.setData( xx, yy, zz );

  SECTION("Write Data")
  {
    std::ofstream out;

    out.open("TPS-input.txt");
    for(int i = 0; i < nx*ny; i++)
      out << xx(i) << " " << yy(i) << " " << zz(i) << "\n";
    out.close();

    out.open("TPS-output.txt");
    for(int i = 0; i < 5*nx; i++)
    {
      for(int j = 0; j < 5*ny; j++)
      {
        out << xmin + dx*i/5. << " " << ymin + dy*j/5. << " " << interp(xmin+dx*i/5., ymin+dy*j/5.) << "\n";
      }
      out << "\n";
    }
    out.close();

  }



  SECTION("Interpolation")
  {
    CHECK( interp(0,0)   == Approx(f(0,0)).epsilon(0.00002 )) ;
    CHECK( interp(1,2)   == Approx(f(1,2)).epsilon(0.00002 )) ;
    CHECK( interp(2,1)   == Approx(f(2,1)).epsilon(0.00002 )) ;
    CHECK( interp(2,-1)  == Approx(f(2,-1)).epsilon(0.00002 )) ;
    CHECK( interp(8,3)   == Approx(f(8,3)).epsilon(0.00002 )) ;

    CHECK( interp(-2,-1) == Approx(0).epsilon(0.00002 )) ;
    CHECK( interp(10,3)  == Approx(0).epsilon(0.00002 )) ;
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
  CHECK( interp.integral( 1, 3, 0, 2 ) == Approx( integral( 1, 3, 0, 2 ) ).epsilon(0.001) );

  // integral over domain containing some partial elements, but each corner is in a different element
  CHECK( interp.integral( 1.5, 3.5, 0.5, 2.5 ) == Approx( integral( 1.5, 3.5, 0.5, 2.5 ) ).epsilon(0.001) );

  // integral over domain with the x limits contained in the same element
  CHECK( interp.integral( 1.5, 1.75, 0, 2 ) == Approx( integral( 1.5, 1.75, 0, 2 ) ).epsilon(0.001) );

  // same thing, but offset y limits
  CHECK( interp.integral( 1.5, 1.75, 0.5, 2.5 ) == Approx( integral( 1.5, 1.75, 0.5, 2.5 ) ).epsilon(0.001) );

  // integral over domain with the y limits contained in the same element
  CHECK( interp.integral( 1, 3, 0.5, 0.75 ) == Approx( integral( 1, 3, 0.5, 0.75 ) ).epsilon(0.001) );

  // same thing, but offset x limits
  CHECK( interp.integral( 1.5, 3.5, 0.5, 0.75 ) == Approx( integral( 1.5, 3.5, 0.5, 0.75 ) ).epsilon(0.001) );

  // integral over domain with both x and y limits in the same element
  CHECK( interp.integral( 1.5, 1.75, 0.5, 0.75 ) == Approx( integral( 1.5, 1.75, 0.5, 0.75 ) ).epsilon(0.001) );



  }



}

TEST_CASE( "ThinPlateSplineInterpolator Tests - Oscillating Data", "[thin plate spline]" ) {

  _2D::TestThinPlateSplineInterp interp;

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

  _2D::ThinPlateSplineInterpolator<double>::VectorType xx(nx*ny), yy(nx*ny), zz(nx*ny);

  auto f  = [](double x, double y){return sin(x)*sin(y);};

  for( int i = 0; i < nx*ny; i++)
  {
    // gnuplot format is essentially row-major
    xx(i) = xmin+dx*(i/ny);
    yy(i) = ymin+dy*(i%ny);
    zz(i) = f(xx(i),yy(i));
  }
  interp.setData( xx, yy, zz );

  SECTION("Write Data")
  {
    std::ofstream out;

    out.open("TPS-oscillations-input.txt");
    for(int i = 0; i < nx*ny; i++)
      out << xx(i) << " " << yy(i) << " " << zz(i) << "\n";
    out.close();

    out.open("TPS-oscillations-output.txt");
    for(int i = 0; i < 5*nx; i++)
    {
      for(int j = 0; j < 5*ny; j++)
      {
        out << xmin + dx*i/5. << " " << ymin + dy*j/5. << " " << interp(xmin+dx*i/5., ymin+dy*j/5.) << "\n";
      }
      out << "\n";
    }
    out.close();

  }

  SECTION("Interpolation")
  {
    CHECK( interp(0,0)   == Approx(f(0,0)).epsilon(0.00002 )) ;
    CHECK( interp(1,2)   == Approx(f(1,2)).epsilon(0.00002 )) ;
    CHECK( interp(2,1)   == Approx(f(2,1)).epsilon(0.00002 )) ;
    CHECK( interp(2,-1)  == Approx(f(2,-1)).epsilon(0.00002 )) ;
    CHECK( interp(8,3)   == Approx(f(8,3)).epsilon(0.00002 )) ;

    CHECK( interp(-2,-1) == Approx(0).epsilon(0.00002 )) ;
    CHECK( interp(10,3)  == Approx(0).epsilon(0.00002 )) ;

    CHECK( interp(M_PI/2,M_PI/2)  == Approx(sin(M_PI/2)*sin(M_PI/2)).epsilon(0.01)) ;
  }






}
