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
  ny = 5;

  xmin = -1;
  xmax = 8;

  ymin = -1;
  ymax = 3;

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
      if( abs(z) > 1e-10 )
        BOOST_CHECK_CLOSE( z, interp(x,y), 0.1 );
    }
  }

  // \int f(x,y) dx dy = x^2 y^2 / 4 + x^2 y + 3 x y^2 / 2 = g(x,y) |_xa^xb | _ya^yb
  // 
  // = ( g(xb,yb) - g(xb,ya) ) - ( g(xa,yb) - g(xa,ya) ) = g(xb,yb) - g(xb,ya) - g(xa,yb) + g(xa,ya)
  //
  
  double sum, xa, xb, ya, yb;

#define integrate(_xa,_xb,_ya,_yb) \
  sum = 0; \
  xa = _xa; \
  xb = _xb; \
  ya = _ya; \
  yb = _yb; \
  x = xb; y = yb; \
  sum += x*x*y*y/4 + x*x*y + 3*x*y*y/2; \
  x = xb; y = ya; \
  sum -= x*x*y*y/4 + x*x*y + 3*x*y*y/2; \
  x = xa; y = yb; \
  sum -= x*x*y*y/4 + x*x*y + 3*x*y*y/2; \
  x = xa; y = ya; \
  sum += x*x*y*y/4 + x*x*y + 3*x*y*y/2;


  // function we are interpolating is defined over
  // x -> [-1,8] (increments of 1)
  // y -> [-1,3] (increments of 1)

  // integral over domain containing only whole element
  integrate(                               1, 3, 0, 2 );
  BOOST_CHECK_CLOSE( sum, interp.integral( 1, 3, 0, 2 ), 0.1 );

  // integral over domain containing some partial elements, but each corner is in a different element
  integrate(                               1.5, 3.5, 0.5, 2.5 );
  BOOST_CHECK_CLOSE( sum, interp.integral( 1.5, 3.5, 0.5, 2.5 ), 0.1 );

  // integral over domain with the x limits contained in the same element
  integrate(                               1.5, 1.75, 0, 2 );
  BOOST_CHECK_CLOSE( sum, interp.integral( 1.5, 1.75, 0, 2 ), 0.1 );

  // same thing, but offset y limits
  integrate(                               1.5, 1.75, 0.5, 2.5 );
  BOOST_CHECK_CLOSE( sum, interp.integral( 1.5, 1.75, 0.5, 2.5 ), 0.1 );

  // integral over domain with the y limits contained in the same element
  integrate(                               1, 3, 0.5, 0.75 );
  BOOST_CHECK_CLOSE( sum, interp.integral( 1, 3, 0.5, 0.75 ), 0.1 );

  // same thing, but offset x limits
  integrate(                               1.5, 3.5, 0.5, 0.75 );
  BOOST_CHECK_CLOSE( sum, interp.integral( 1.5, 3.5, 0.5, 0.75 ), 0.1 );

  // integral over domain with both x and y limits in the same element
  integrate(                               1.5, 1.75, 0.5, 0.75 );
  BOOST_CHECK_CLOSE( sum, interp.integral( 1.5, 1.75, 0.5, 0.75 ), 0.1 );


  delete[] X;
  delete[] Y;
  delete[] Z;



}

BOOST_AUTO_TEST_SUITE_END()
