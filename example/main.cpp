#include <iostream>
#include <vector>

#include "boost/lexical_cast.hpp"

#include "interpLib.h"

int main(int argc, char *argv[])
{
    //set up some example data:
    double xmin = 0;
    double xmax = 10;
    const int N = 8;
    double dx = (xmax - xmin)/(N-1);

    std::vector<double> X, Y;

    X.resize(N);
    Y.resize(N);

    for (int i = 0; i < N; ++i)
    {
        double xtmp = xmin+i*dx;
        X.at(i) = xtmp;
        Y.at(i) = xtmp*xtmp;
    } 

    //create the interpolation object
    SplineInterp exampleInterp;

    //give data to the interpolation object
    exampleInterp.setData( X, Y);

    //interpolate an example value from user input
    if( argc > 1 )
    {
        std::cout <<  exampleInterp(boost::lexical_cast<double>(argv[1])) << std::endl;
        std::cout << "Derivative: " << exampleInterp.derivative(boost::lexical_cast<double>(argv[1])) << std::endl;
    }

    return 0;
}
