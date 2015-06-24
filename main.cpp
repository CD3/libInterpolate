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

    //std::cout << testInterp.interpAt( 3.00 ) << std::endl;
    //std::cout << testInterp.interpAt( 3.10 ) << std::endl;
    //std::cout << testInterp.interpAt( 3.20 ) << std::endl;
    //std::cout << testInterp.interpAt( 3.30 ) << std::endl;
    //std::cout << testInterp.interpAt( 3.40 ) << std::endl;
    //std::cout << testInterp.interpAt( 3.50 ) << std::endl;
    //std::cout << testInterp.interpAt( 3.60 ) << std::endl;

    //std::cout << argc << std::endl;

    //if( argc > 1 )
    //{
        //std::cout << testInterp.interpAt(boost::lexical_cast<double>(argv[1]) ) << std::endl;;
    //}
    return 0;
}
