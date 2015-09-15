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


  delete[] X;
  delete[] Y;
  delete[] Z;



}

BOOST_AUTO_TEST_SUITE_END()
