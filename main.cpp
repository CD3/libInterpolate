#include <iostream>
#include "splineInterp/interpLib.h"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/lexical_cast.hpp"

namespace ublas=boost::numeric::ublas;

int main(int argc, char *argv[])
{

    std::string path = "test.txt";
    SplineInterp testInterp(path);

    if( argc > 1 )
    {
        std::cout << testInterp.interpAt(boost::lexical_cast<double>(argv[1]) ) << std::endl;;
    }
    return 0;
}
