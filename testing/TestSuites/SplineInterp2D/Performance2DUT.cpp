#include <iostream>
#include <fstream>
#include <vector>
#include <boost/progress.hpp>

#define BOOST_TEST_MODULE PerformanceUT
#include <boost/test/included/unit_test.hpp>

#include "interpLib.hpp"

BOOST_AUTO_TEST_SUITE(Performance2DUT)

BOOST_AUTO_TEST_CASE(InitData)
{
    int n = 10000;
    double xmin,xmax,dx;
    std::vector<double> x, y, z;
    xmin = 0;
    xmax = 10;
    dx = (xmax - xmin) / (n - 1);

    //set up test data
    for (int i = 0; i < n; ++i)
    {
        double a = xmin+i*dx;
        for (int j = 0; j < n; ++j)
        {
            double b = xmin+j*dx;
            
            x.push_back(a);
            y.push_back(b);
            z.push_back(a*b);
        }
    }

    boost::progress_timer* timer;

    timer = new boost::progress_timer();
    std::cout<<"2D Data Init:"<<std::endl;
    SplineInterp2D<double> customInterp;
    customInterp.setData(x,y,z);
    delete timer;

}

BOOST_AUTO_TEST_SUITE_END()
