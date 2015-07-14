#include <iostream>
#include <fstream>
#include <vector>

#define BOOST_TEST_MODULE LinearInterpUT
#include <boost/test/included/unit_test.hpp>

#include "interpLib.h"

BOOST_AUTO_TEST_SUITE(CubicSplineInterpUT)

BOOST_AUTO_TEST_CASE(testInterpVecs)
{
    
    //set up the class
    SplineInterp testInterp;
    int n = 10;
    std::vector<double> x, y;

    //set up test data
    for (int i = 0; i < n; ++i)
    {
        x.push_back(i);
        y.push_back(i*i);
    }

    //put data into the class
    testInterp.setData(x,y);

    //test case
    BOOST_CHECK_CLOSE(testInterp(2.5), 6.2566, 1e-3);
}

BOOST_AUTO_TEST_CASE(testSinInterp)
{
    int n = 30;
    double xmin,xmax,dx;
    std::vector<double> x(n), y(n);
    xmin = 0;
    xmax = 10;
    dx = (xmax - xmin) / (n - 1);

    //set up test data
    for (int i = 0; i < n; ++i)
    {
        x[i] = xmin+i*dx;
        y[i] = sin(x[i]);
    }

    // get an interpolator
    SplineInterp testInterp;
    testInterp.setData(x,y);

    ////test case
    int res = 10;
    for(int i = 0; i < res*n; i++)
    {
      double xx = xmin -.1 + i*dx/res;
      if( sin(xx) > 0.01 ) // only check values much greater than zero
        BOOST_CHECK_CLOSE( (xx > 0 && xx < 10) ? sin(xx) : 0, testInterp(xx), 1);
    }
}

BOOST_AUTO_TEST_CASE(testSinDerivative)
{
    int n = 30;
    double xmin,xmax,dx;
    std::vector<double> x(n), y(n);
    xmin = 0;
    xmax = 10;
    dx = (xmax - xmin) / (n - 1);

    //set up test data
    for (int i = 0; i < n; ++i)
    {
        x[i] = xmin+i*dx;
        y[i] = sin(x[i]);
    }

    // get an interpolator
    SplineInterp testInterp;
    testInterp.setData(x,y);

    ////test case
    int res = 10;
    for(int i = 0; i < res*n; i++)
    {
      double xx = xmin -.1 + i*dx/res;
      if( sin(xx) > 0.01 ) // only check values much greater than zero
        BOOST_CHECK_CLOSE( (xx > 0 && xx < 10) ? cos(xx) : 0, testInterp.derivative(xx), 2.9);
    }
}

BOOST_AUTO_TEST_CASE(testSinIntegral)
{
    int n = 30;
    double xmin,xmax,dx;
    std::vector<double> x(n), y(n);
    xmin = 0;
    xmax = 10;
    dx = (xmax - xmin) / (n - 1);

    //set up test data
    for (int i = 0; i < n; ++i)
    {
        x[i] = xmin+i*dx;
        y[i] = sin(x[i]);
    }

    // get an interpolator
    SplineInterp testInterp;
    testInterp.setData(x,y);

    ////test case
    int res = 10;
    for(int i = 0; i < res*n; i++)
    {
      double xx = xmin -.1 + i*dx/res;
      if( abs(cos(xx)) > 0.01 ) // only check values much greater than zero
        BOOST_CHECK_CLOSE( -cos(xx), testInterp.integral(0,xx)-1, 1);
    }

    for(int i = 0; i < res*n; i++)
    {
      double xx = xmin -.1 + i*dx/res;
      if( abs(cos(xx)) > 0.01 ) // only check values much greater than zero
        BOOST_CHECK_CLOSE( cos(xx), testInterp.integral(xx,0)-1, 1);
    }
}

BOOST_AUTO_TEST_CASE(testDefinateIntegrals)
{
    int n = 30;
    double xmin,xmax,dx;
    std::vector<double> x(n), y(n);
    xmin = 0;
    xmax = M_PI;
    dx = (xmax - xmin) / (n - 1);

    //set up test data
    for (int i = 0; i < n; ++i)
    {
        x[i] = xmin+i*dx;
        y[i] = sin(x[i]);
    }

    // get an interpolator
    SplineInterp testInterp;
    testInterp.setData(x,y);


    BOOST_CHECK_CLOSE( 2, testInterp.integral(), 1);
}

BOOST_AUTO_TEST_SUITE_END()
