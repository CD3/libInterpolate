#include <iostream>
#include <fstream>
#include <vector>
#include <boost/progress.hpp>

#define BOOST_TEST_MODULE PerformanceUT
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

    //std::cout<<"Eigen solver"<<std::endl;
    //timer = new boost::progress_timer();
    //eigen::SplineInterp<double> eigenInterp;
    //eigenInterp.setData(x,y);
    //delete timer;

    timer = new boost::progress_timer();
    std::cout<<"2D Data Init:"<<std::endl;
    eigen::SplineInterp2D<double> customInterp;
    customInterp.setData(x,y,z);
    delete timer;

    //timer = new boost::progress_timer();
    //delete timer;

}

BOOST_AUTO_TEST_SUITE_END()
