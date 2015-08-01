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

template<class Real>
class SplineInterp2D
{
    /*
     * INARTDI is Not A Real Two Dimensional Interpolator
     *
     * I really would like to work on implementing a better 2D interpolation method that works on scattered data.
     * As this works by rearranging blocks of the input data and calling the 1D Spline interpolator multiple times,
     * it is probably not going to be nearly as fast as a real 2D interpolation method would be.
     *
     *
     * Usage:
     *  The setData function expects that the x, y, and z vectors are straight from a gnuplot style 2D data file.
     *
     */

    private:
        std::vector<Real> d2y;
        std::vector< std::vector<Real> > d2x, d2z;

    public:
        void setData( std::vector<Real> &x, std::vector<Real> &y, std::vector<Real> &z );
        void setData( size_t n, Real *_x, Real *_y, Real *_z );

        Real operator() (Real _x, Real _y);
};

#include "interpLib_imp.hpp"
#include "interpLib2D_imp.hpp"

#endif
