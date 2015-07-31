#ifndef interplib_hpp
#define interplib_hpp

#include <fstream>
#include <string>
#include <vector>

#ifdef USE_EIGEN
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#endif

template<class Real>
class SplineInterp
{

    private:
        //this is the length of the data vector
        int n;

        //These are the vectors that the data will be stored in
        std::vector<Real> X;
        std::vector<Real> Y;

        //When these coefficients are determined,
        //they can be used to generate the splines that are evaluated for the interpolation
        std::vector<Real> a;
        std::vector<Real> b;

        void initCoefficients();
        void solveForCoefficients();

    protected:

    public:
        SplineInterp() {};
        ~SplineInterp() {};

        //overload the () operator to return an interpolated value
        Real operator()( Real x );

        Real derivative( Real x );
  
        
        Real integral( Real _a, Real _b);

        Real integral( );

        void setData( size_t _n, Real *_x, Real *_y );

        void setData( std::vector<Real> &x, std::vector<Real> &y );
};

#include "interpLib_imp.hpp"

#endif
