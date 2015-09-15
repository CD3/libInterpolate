#include <iostream>
#include <fstream>
#include <vector>
#include <boost/progress.hpp>

#define BOOST_TEST_MODULE PerformanceUT
#include <boost/test/included/unit_test.hpp>

#include "interpLib.hpp"

BOOST_AUTO_TEST_SUITE(PerformanceUT)

BOOST_AUTO_TEST_CASE(InitBigData)
{
    int n = 30000;
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


    boost::progress_timer* timer;

    std::cout<<"Custom solver"<<std::endl;
    timer = new boost::progress_timer();
    SplineInterp<double> customInterp;
    customInterp.setData(x,y);
    delete timer;

}

BOOST_AUTO_TEST_SUITE_END()
