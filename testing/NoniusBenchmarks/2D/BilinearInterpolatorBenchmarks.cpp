#include "nonius/nonius.h++"


#include "Interpolators/_2D/BilinearInterpolator.hpp"


NONIUS_BENCHMARK("BilinearInterpolator Small Data Set",
//[](nonius::chronometer meter)
[]()
{
  _2D::BilinearInterpolator<double> interp;
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

  for(int i = 0; i < 100; i++)
    interp(2,2);

})












