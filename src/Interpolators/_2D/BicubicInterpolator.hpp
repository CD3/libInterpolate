#ifndef Interpolators__2D_BicubicInterpolator_hpp
#define Interpolators__2D_BicubicInterpolator_hpp

/** @file BicubicInterpolator.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 12/27/16
  */

#include "InterpolatorBase.hpp"

/** @class 
  * @brief Cubic spline interpolation for for 2D functions.
  * @author C.D. Clark III, Aaron Hoffman
  *
  * This class implements the bicubic spline interpolation method.
  * It is essentially the 2D equivalent of cubic splines for 1D. 
  *
  */

namespace _2D {

template<class Real>
class BicubicInterpolator : public InterpolatorBase<Real>
{
  public:
    // typedefs
    typedef typename InterpolatorBase<Real>::VectorType VectorType;
    typedef typename InterpolatorBase<Real>::MapType MapType;
    typedef typename InterpolatorBase<Real>::GradientType GradientType;

    // types used to view data as 2D coordinates
    typedef typename Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> MatrixType;
    typedef Eigen::Map<VectorType,Eigen::Unaligned,Eigen::InnerStride<Eigen::Dynamic               >> _2DVectorView;
    typedef Eigen::Map<MatrixType,Eigen::Unaligned,     Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>> _2DMatrixView;

    typedef Eigen::Matrix<Real,4,4 > Matrix44;
    typedef Eigen::Array< Matrix44, Eigen::Dynamic, Eigen::Dynamic > Matrix44Array;
    typedef Eigen::Matrix<Real,4,1 > ColVector4;
    typedef Eigen::Matrix<Real,1,4 > RowVector4;

    // methods required by the interface
    virtual Real operator()( Real x, Real y ) const;
    virtual void setData( size_t _n, Real *x, Real *y, Real *z, bool deep_copy = true );
    using InterpolatorBase<Real>::setData;

    // additional methods
    virtual void setData( std::vector<Real> &x, std::vector<Real> &y, std::vector<Real> &z, bool deep_copy = true );
    virtual void setData( VectorType  &x, VectorType &y, VectorType &z, bool deep_copy = true );

  protected:
    using InterpolatorBase<Real>::xv;
    using InterpolatorBase<Real>::yv;
    using InterpolatorBase<Real>::zv;

    // these maps are used to view the x,y,z data as two coordinate vectors and a function matrix, instead of three vectors.
    std::shared_ptr<_2DVectorView> X,Y;
    std::shared_ptr<_2DMatrixView> Z;
    
    void calcCoefficients();
    Matrix44Array a; // naming convention used by wikipedia article (see Wikipedia https://en.wikipedia.org/wiki/Bicubic_interpolation)



};

template<class Real>
void
BicubicInterpolator<Real>::setData( size_t n, Real *x, Real *y, Real *z, bool deep_copy )
{
  InterpolatorBase<Real>::setData(n,x,y,z,deep_copy);
  calcCoefficients();
}

template<class Real>
void
BicubicInterpolator<Real>::setData( std::vector<Real> &x, std::vector<Real> &y, std::vector<Real> &z, bool deep_copy )
{
  InterpolatorBase<Real>::setData(x,y,z,deep_copy);
  calcCoefficients();
}

template<class Real>
void
BicubicInterpolator<Real>::setData( VectorType  &x, VectorType &y, VectorType &z, bool deep_copy )
{
  InterpolatorBase<Real>::setData(x,y,z,deep_copy);
  calcCoefficients();
}

template<class Real>
void
BicubicInterpolator<Real>::calcCoefficients()
{
  // setup 2D view of the data
  // We need to figure out what the x and y dimensions are.
  int Nx = 0, Ny = 0;
  int N = xv->size();
  // Ny will be the number of elements that have the same x coordinate
  Real xlast = (*xv)(0);
  while( Ny < N-1 && fabs((*xv)(Ny)-xlast) < 1e-40 )
    Ny++;
  Nx = N/Ny;

  // consecutive values in the x data are separated by Ny, so this is the inner stride for X
  X.reset( new _2DVectorView( &(*xv)(0), Nx, Eigen::InnerStride<Eigen::Dynamic>(Ny) ) );

  // consecutive values in the y data are next to each other, so the stride is just 1
  Y.reset( new _2DVectorView( &(*yv)(0), Ny, Eigen::InnerStride<Eigen::Dynamic>(1) ) );

  // Eigen defaults to COLUMN MAJOR
  // consecutive elements in a column are separated by Ny (this is the inner stride)
  // consecutive elements in a row are located next to each other (this is the outer stride)
  // Stride object takes outer,inner as aruments.
  Z.reset( new _2DMatrixView( &(*zv)(0), Nx, Ny, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(1,Ny) ) );

  // Interpolation will be done by multiplying the coordinates by coefficients.

  a = Matrix44Array( X->size()-1, Y->size()-1 );

  // We are going to precompute the interpolation coefficients so
  // that we can interpolate quickly. This requires a 4x4 matrix for each "patch".
  Matrix44 Left, Right;

  Left <<  1,  0,  0,  0,
           0,  0,  1,  0,
          -3,  3, -2, -1,
           2, -2,  1,  1;

  Right<<  1,  0, -3,  2,
           0,  0,  3, -2,
           0,  1, -2,  1,
           0,  0, -1,  1;

  for(int i = 0; i < X->size() - 1; i++)
  {
    for( int j = 0; j < Y->size() - 1; j++)
    {
      Matrix44 F;

      Real f00,   f01,   f10,   f11;
      Real fx00,  fx01,  fx10,  fx11;
      Real fy00,  fy01,  fy10,  fy11;
      Real fxy00, fxy01, fxy10, fxy11;

      Real fm, fp;
      int im, ip, jm, jp;

      int iN = X->size();
      int jN = Y->size();

      // function values
      f00 = (*Z)(i    ,j    ); // <<<<<<
      f01 = (*Z)(i    ,j + 1); // <<<<<<
      f10 = (*Z)(i + 1,j    ); // <<<<<<
      f11 = (*Z)(i + 1,j + 1); // <<<<<<

      // need to calculate function values and derivatives
      // at each corner.
      //
      // note: interpolation algorithm is derived for the unit square.
      // so we need to take the derivatives assuming X(i+1) - X(i) = Y(j+1) - Y(j) = 1
      
      Real xL = (*X)(i+1) - (*X)(i);
      Real yL = (*Y)(j+1) - (*Y)(j);
      Real dx, dy;

      // x derivatives

      im = std::max(i-1,0);
      ip = std::min(i+1,iN-1);

      dx = ((*X)(ip) - (*X)(im))/xL;

      fp = (*Z)(ip,j);
      fm = (*Z)(im,j);
      fx00 = (fp - fm) / dx; // <<<<<<

      fp = (*Z)(ip,j+1);
      fm = (*Z)(im,j+1);
      fx01 = (fp - fm) / dx; // <<<<<<


      im = std::max(i,0);
      ip = std::min(i+2,iN-1);

      dx = ((*X)(ip) - (*X)(im))/xL;

      fp = (*Z)(ip,j);
      fm = (*Z)(im,j);
      fx10 = (fp - fm) / dx; // <<<<<<

      fp = (*Z)(ip,j+1);
      fm = (*Z)(im,j+1);
      fx11 = (fp - fm) / dx; // <<<<<<


      // y derivatives

      jm = std::max(j-1,0);
      jp = std::min(j+1,jN-1);

      dy = ((*Y)(jp) - (*Y)(jm))/yL;

      fp = (*Z)(i,jp);
      fm = (*Z)(i,jm);
      fy00 = (fp - fm) / dy; // <<<<<<

      fp = (*Z)(i+1,jp);
      fm = (*Z)(i+1,jm);
      fy10 = (fp - fm) / dy; // <<<<<<


      jm = std::max(j,0);
      jp = std::min(j+2,jN-1);

      (*Y)(jp) = (*Y)(jp);
      (*Y)(jm) = (*Y)(jm);
      dy = ((*Y)(jp) - (*Y)(jm))/yL;

      fp = (*Z)(i,jp);
      fm = (*Z)(i,jm);
      fy01 = (fp - fm) / yL; // <<<<<<

      fp = (*Z)(i+1,jp);
      fm = (*Z)(i+1,jm);
      fy11 = (fp - fm) / yL; // <<<<<<

      // xy derivatives

      im = std::max(i-1,0);
      ip = std::min(i+1,iN-1);
      jm = std::max(j-1,0);
      jp = std::min(j+1,jN-1);

      dx = ((*X)(ip) - (*X)(im)) / xL;

      dy = ((*Y)(jp) - (*Y)(jm)) / yL;

      fp = ((*Z)(ip,jp) - (*Z)(im,jp))/dx;
      fm = ((*Z)(ip,jm) - (*Z)(im,jm))/dx;
      fxy00 = (fp - fm) / dy; // <<<<<<


      jm = std::max(j,0);
      jp = std::min(j+2,jN-1);

      dy = ((*Y)(jp) - (*Y)(jm))/yL;

      fp = ((*Z)(ip,jp) - (*Z)(im,jp))/dx;
      fm = ((*Z)(ip,jm) - (*Z)(im,jm))/dx;
      fxy01 = (fp - fm) / dy; // <<<<<<


      im = std::max(i,0);
      ip = std::min(i+2,iN-1);
      jm = std::max(j-1,0);
      jp = std::min(j+1,jN-1);

      dx = ((*X)(ip) - (*X)(im)) / xL;

      dy = ((*Y)(jp) - (*Y)(jm)) / yL;

      fp = ((*Z)(ip,jp) - (*Z)(im,jp))/dx;
      fm = ((*Z)(ip,jm) - (*Z)(im,jm))/dx;
      fxy10 = (fp - fm) / dy; // <<<<<<

      jm = std::max(j,0);
      jp = std::min(j+2,jN-1);

      dy = ((*Y)(jp) - (*Y)(jm)) / yL;

      fp = ((*Z)(ip,jp) - (*Z)(im,jp))/dx;
      fm = ((*Z)(ip,jm) - (*Z)(im,jm))/dx;
      fxy11 = (fp - fm) / dy; // <<<<<<


      F <<  f00,  f01,  fy00,  fy01,
            f10,  f11,  fy10,  fy11,
           fx00, fx01, fxy00, fxy01,
           fx10, fx11, fxy10, fxy11;

      a(i,j) = Left * F * Right;
    }
  }
}

template<class Real>
Real
BicubicInterpolator<Real>::operator()( Real x, Real y ) const
{
  InterpolatorBase<Real>::checkData();
  
  // no extrapolation...
  if( x < this->xv->minCoeff()
   || x > this->xv->maxCoeff()
   || y < this->yv->minCoeff()
   || y > this->yv->maxCoeff() )
  {
    return 0;
  }


  // find the x index that is just to the LEFT of x
  int i  = Utils::index_last_lt( x, *X );
  if(i < 0)
    i = 0;

  // find the y index that is just BELOW y
  int j  = Utils::index_last_lt( y, *Y );
  if(j < 0)
    j = 0;
  

  Real xL = (*X)(i+1) - (*X)(i);
  Real yL = (*Y)(j+1) - (*Y)(j);

  // now, create the coordinate vectors (see Wikipedia https://en.wikipedia.org/wiki/Bicubic_interpolation)
  RowVector4 vx;
  vx[0] = 1;                                   // x^0
  vx[1] = (x - (*X)(i))/xL;                    // x^1
  vx[2] = vx[1] * vx[1];                       // x^2
  vx[3] = vx[2] * vx[1];                       // x^3

  ColVector4 vy;
  vy[0] = 1;                                   // y^0
  vy[1] = (y - (*Y)(j))/yL;                    // y^1
  vy[2] = vy[1] * vy[1];                       // y^2
  vy[3] = vy[2] * vy[1];                       // y^3


  // interpolation is just x*a*y

  return vx*a(i,j)*vy;

}





}

#endif // include protector
