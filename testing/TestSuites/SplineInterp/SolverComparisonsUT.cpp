#include <iostream>
#include <fstream>
#include <vector>
#include <boost/progress.hpp>

#define BOOST_TEST_MODULE SolverComparisonsUT
#include <boost/test/included/unit_test.hpp>

#include "interpLib.hpp"

BOOST_AUTO_TEST_SUITE(SolverComparisonsUT)

BOOST_AUTO_TEST_CASE(UnifromGrid)
{
    int n = 3000;
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



    SplineInterp<double> customInterp;
    customInterp.setData(x,y);

}

BOOST_AUTO_TEST_SUITE_END()
