#include <iostream>
#include <fstream>
#include <vector>
#include <boost/progress.hpp>

#define BOOST_TEST_MODULE SolverComparisonsUT
#include <boost/test/included/unit_test.hpp>

namespace custom {
#include "interpLib.hpp"
}

namespace eigen {
#undef interplib_hpp
#define USE_EIGEN
#include "interpLib.hpp"
#undef USE_EIGEN
}

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



    eigen::SplineInterp<double> eigenInterp;
    eigenInterp.setData(x,y);

    custom::SplineInterp<double> customInterp;
    customInterp.setData(x,y);

    n = 10;
    dx = (xmax - xmin) / (n - 1);
    for(int i = 0; i < n; i++)
    {
      BOOST_CHECK_CLOSE( eigenInterp(xmin * dx*i), customInterp(xmin * dx*i), 0.01 );
    }

}

BOOST_AUTO_TEST_SUITE_END()