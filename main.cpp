#include <iostream>
#include <fstream>
#include <iomanip>

#include "interpLib/interpLib.h"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/lexical_cast.hpp"

namespace ublas=boost::numeric::ublas;

int main(int argc, char *argv[])
{

    std::string path = "test.txt";
    SplineInterp testInterp;
    testInterp.setData(path);

    std::cout << testInterp(9.2) << std::endl;

    //if( argc > 1 )
    //{
        //std::cout << std::setprecision(12) <<  testInterp(boost::lexical_cast<double>(argv[1]) ) << std::endl;;
    //}

    return 0;
}
