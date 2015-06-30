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

BOOST_AUTO_TEST_SUITE_END()
