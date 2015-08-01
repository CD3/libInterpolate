#include <iostream>
#include <vector>

#include "boost/lexical_cast.hpp"

#include "interpLib.hpp"

int main(int argc, char *argv[])
{
    //set up some example data:
    double xmin = 0;
    double xmax = 10;
    int n = 5;
    double dx = (xmax - xmin)/(n-1);

    double ymin = 0;
    double ymax = 5;
    int m = 6;
    double dy = (ymax-ymin)/(m-1);

    std::vector<double> x,y,z;

    double a,b;
    for (int i = 0; i < n; ++i)
    {
        a = xmin+i*dx;
        for (int j = 0; j < m; ++j)
        {
            b = ymin+j*dy;

            x.push_back(a);
            y.push_back(b);
            z.push_back(a*b);
        }
    }

    //double xlast = x[0];
    //for (size_t i = 0; i < x.size(); ++i)
    //{
        //if( fabs(xlast-x[i]) > 1e-4 )
            //std::cout << std::endl;
        //std::cout << x[i] << "    " << y[i] << "    " << z[i] << std::endl;
        //xlast = x[i];
    //}


    SplineInterp2D<double> testInterp;
    testInterp.setData(x,y,z);

    //std::cout << testInterp( boost::lexical_cast<double>(argv[1]), boost::lexical_cast<double>(argv[2]) ) << std::endl;
    return 0;
}
