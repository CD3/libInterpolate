#ifndef INTERPLIB_H
#define INTERPLIB_H value

#include <fstream>
#include <string>
#include <vector>
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/lu.hpp"
#include "boost/filesystem.hpp"

namespace fs=boost::filesystem;
namespace ublas=boost::numeric::ublas;

class SplineInterp
{

    private:
        //this is the length of the input file
        int n;

        //These are the vectors that the data will be stored int
        ublas::vector<double> X;
        ublas::vector<double> Y;

        //these will hold intermediate data for the interpolation method
        ublas::matrix<double> A;
        ublas::vector<double> B;

        //these are the "end goal" from the interpolation method. When these coefficients are determined,
        //they can be used to generate the splines that the interpolations are calculated at
        ublas::vector<double> a;
        ublas::vector<double> b;

        void initialSetup();

    protected:

    public:
        SplineInterp(std::string dataFile);
        SplineInterp(std::vector<double> x, std::vector<double> y);
        ~SplineInterp();

        //int getN();
        double interpAt(double targetValue);

};

#endif
