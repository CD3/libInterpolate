#ifndef interplib_hpp
#define interplib_hpp

#include <fstream>
#include <string>
#include <vector>

#include "matrixBuild.hpp"

#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/lu.hpp"
#include "boost/filesystem.hpp"

namespace ublas=boost::numeric::ublas;
namespace fs=boost::filesystem;

template<class Real>
class SplineInterp
{

    private:
        //this is the length of the input file
        int n;

        //These are the vectors that the data will be stored in
        ublas::vector<Real> X;
        ublas::vector<Real> Y;

        //When these coefficients are determined,
        //they can be used to generate the splines that are evaluated for the interpolation
        ublas::vector<Real> a;
        ublas::vector<Real> b;

        void initCoefficients();

    protected:

    public:
        SplineInterp() {};
        ~SplineInterp() {};

        //overload the () operator to return an interpolated value
        Real operator()( Real x );

        //the square bracket one is probably unnecessary
        Real operator[]( Real x );

        Real derivative( Real x );
  
        
        Real integral( Real _a, Real _b);

        Real integral( );

        void setData( size_t _n, Real *_x, Real *_y );

        void setData( std::vector<Real> x, std::vector<Real> y );
};

#include "interpLib_imp.hpp"

#endif
