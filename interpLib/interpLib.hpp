#ifndef interplib_hpp
#define interplib_hpp

#include <fstream>
#include <string>
#include <vector>

#include <boost/shared_array.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

using namespace Eigen;

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

/**
 * A 2D interpolation class that does bilinear interpoaltion
 */
template<class Real>
class BilinearInterp2D
{

    private:
      //some typedefs for using Eigen vectors and matrixes
      typedef Matrix<Real,Dynamic,1               > VectorType;
      typedef Matrix<Real,Dynamic,Dynamic,RowMajor> MatrixType;
      typedef InnerStride<> InnerStrideType;
      typedef Stride<Dynamic,Dynamic> StrideType;
      typedef Map< VectorType,0,InnerStride<> > VectorMap;
      typedef Map< MatrixType,0,StrideType > MatrixMap;

      typedef Matrix<Real,2,2> Matrix22;
      typedef Array<Matrix22, Dynamic, Dynamic> Matrix22Array;

      typedef Matrix<Real,2,1 > ColVector2;
      typedef Matrix<Real,1,2 > RowVector2;

      // data we are interpolating from
      VectorType X;
      VectorType Y;
      MatrixType Z;


      // coeficient matrixes
      Matrix22Array C;


    public:
      void setData( size_t n, Real* x, Real* y, Real* z );
      void setData( std::vector<Real> &x, std::vector<Real> &y, std::vector<Real> &z );
      Real operator() (Real _x, Real _y);
      Real integral( Real xa, Real xb, Real ya, Real yb );
      //Real integral( );

};

#include "interpLib_imp.hpp"
#include "interpLib2D_imp.hpp"

#endif
