#include <iostream>
#include <fstream>
#include <vector>

#define BOOST_TEST_MODULE BilinearInterp2DUT
#include <boost/test/included/unit_test.hpp>

#include "interpLib.hpp"

BOOST_AUTO_TEST_SUITE(BilinearInterp2DUT)

BOOST_AUTO_TEST_CASE(Validation)
{

  int nx, ny;
  double xmin, xmax, dx, x;
  double ymin, ymax, dy, y;
  double z;

  nx = 10;
  ny = 6;

  xmin = 0;
  xmax = 4;

  ymin = -1;
  ymax = 5;

  dx = (xmax - xmin)/(nx - 1);
  dy = (ymax - ymin)/(ny - 1);

  double *X, *Y, *Z;
  X = new double[nx*ny];
  Y = new double[nx*ny];
  Z = new double[nx*ny];
  for(int i = 0; i < nx; i++)
  {
    x = xmin + i*dx;
    for(int j = 0; j < ny; j++)
    {
      y = ymin + j*dy;
      z = x*y + 2*x + 3*y;
      X[i*ny + j] = x;
      Y[i*ny + j] = y;
      Z[i*ny + j] = z;
    }
  }

  BilinearInterp2D<double> interp;
  interp.setData( nx*ny, X, Y, Z );

  int res = 10;
  nx *= res;
  dx /= res;
  ny *= res;
  dy /= res;
  for(int i = 0; i < nx; i++)
  {
    x = xmin + i*dx;
    for(int j = 0; j < ny; j++)
    {
      y = ymin + j*dy;
      z = x*y + 2*x + 3*y;
      BOOST_CHECK_CLOSE( z, interp(x,y), 0.1 );
    }
  }

  // \int f(x,y) dx dy = x^2 y^2 / 4 + x^2 y + 3 x y^2 / 2 = g(x,y) |_xa^xb | _ya^yb
  // 
  // = ( g(xb,yb) - g(xb,ya) ) - ( g(xa,yb) - g(xa,ya) ) = g(xb,yb) - g(xb,ya) - g(xa,yb) + g(xa,ya)
  //
  double sum = 0;
  double xa = 1;
  double xb = 3;
  double ya = 0;
  double yb = 4;

  x = xb; y = yb;
  sum += x*x*y*y/4 + x*x*y + 3*x*y*y/2;
  x = xb; y = ya;
  sum -= x*x*y*y/4 + x*x*y + 3*x*y*y/2;
  x = xa; y = yb;
  sum -= x*x*y*y/4 + x*x*y + 3*x*y*y/2;
  x = xa; y = ya;
  sum += x*x*y*y/4 + x*x*y + 3*x*y*y/2;

  BOOST_CHECK_CLOSE( sum, interp.integral( xa, xb, ya, yb ), 0.1 );


  delete[] X;
  delete[] Y;
  delete[] Z;



}

BOOST_AUTO_TEST_SUITE_END()
