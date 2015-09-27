#include <iostream>
#include <fstream>
#include <vector>

#define BOOST_TEST_MODULE LinearInterpUT
#include <boost/test/included/unit_test.hpp>

#include "interpLib.hpp"

BOOST_AUTO_TEST_SUITE(CubicSplineInterp2DUT)

BOOST_AUTO_TEST_CASE(testInterpVecs)
{
    
    //set up the class
    SplineInterp2D <double>testInterp;
    int n = 10;
    std::vector<double> x, y, z;

    //set up test data
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            x.push_back(i);
            y.push_back(j);
            z.push_back(i*j);
        }
    }

    //put data into the class
    testInterp.setData(x,y,z);

    //test case
    BOOST_CHECK_CLOSE(testInterp(2.5,2.5), 6.25, 1e-4);
}

BOOST_AUTO_TEST_CASE(testSinCosInterp)
{
    int n = 30;
    std::vector<double> x,y,z;
    double xmin, xmax, dx;
    xmin = 0;
    xmax = 10;
    dx = (xmax-xmin)/(n-1);

    double a,b;
    for (int i = 0; i < n; ++i)
    {
        a = xmin+i*dx;
        for (int j = 0; j < n; ++j)
        {
            b = xmin+j*dx;
            x.push_back(a);
            y.push_back(b);
            z.push_back(sin(a)*cos(b));
        }
    }

    SplineInterp2D<double> testInterp;
    testInterp.setData(x,y,z);

    int res = 10;
    for (int i = 0; i < res*n; ++i)
    {
        double xx = xmin - .1 + i*dx/res;
        for (int j = 0; j < res*n; ++j)
        {
            double yy = xmin -.1 + j*dx/res;

            if( sin(xx) > 0.01 && cos(yy) > 0.01 )
            {
                BOOST_CHECK_CLOSE( (xx > 0 && xx < 10 && yy > 0 && yy < 10 ) ? sin(xx)*cos(yy) : 0, testInterp(xx,yy), 1.5);
            }
        }
    }
}


BOOST_AUTO_TEST_SUITE_END()
