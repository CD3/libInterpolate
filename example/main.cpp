#include <iostream>
#include <vector>

#include "boost/lexical_cast.hpp"

#include "interpLib.hpp"

int main(int argc, char *argv[])
{
    //set up some example data:
    double xmin = 0;
    double xmax = 10;
    const int N = 3;
    double dx = (xmax - xmin)/(N-1);

    std::vector<double> X, Y;

    //X.resize(N);
    //Y.resize(N);

    //for (int i = 0; i < N; ++i)
    //{
        //double xtmp = xmin+i*dx;
        //X.at(i) = xtmp;
        //Y.at(i) = xtmp*xtmp;
    //} 

    //using the example points from the wikipedia page.
    //
    X.push_back(-1.0);
    Y.push_back(0.5);

    X.push_back(0);
    Y.push_back(0);

    X.push_back(3);
    Y.push_back(3);

    //create the interpolation object
    SplineInterp <double>exampleInterp;

    //give data to the interpolation object
    exampleInterp.setData( X, Y);

    //interpolate an example value from user input
    //if( argc > 1 )
    //{
        //std::cout <<  exampleInterp(boost::lexical_cast<double>(argv[1])) << std::endl;
        //std::cout << "Derivative: " << exampleInterp.derivative(boost::lexical_cast<double>(argv[1])) << std::endl;
    //}

    return 0;
}
