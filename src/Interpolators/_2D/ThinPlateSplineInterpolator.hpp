#ifndef Interpolators__2D_ThinPlateSplineInterpolator_hpp
#define Interpolators__2D_ThinPlateSplineInterpolator_hpp

/** @file ThinPlateSplineInterpolator.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 12/27/16
  */

#include "InterpolatorBase.hpp"

/** @class 
  * @brief Cubic spline interpolation for for 2D functions.
  * @author C.D. Clark III, Aaron Hoffman
  *
  * This class implements the "Thin Plate Spline" method as derived by David Eberly (https://www.geometrictools.com/Documentation/ThinPlateSplines.pdf)
  * It is essentially the 2D equivalent of cubic splines for 1D. 
  *
  * This method has a performance cost over other methods. Specifically, for N interpolation points, an NxN matrix
  * must be inverted. The matrix inversion is only required during setup, but interpolation will be based on all
  * N points as well, rather than the nearest neighbors. So, it is possible that this method will be slow for
  * large data sets.
  * 
  */

namespace _2D {

template<class Real>
class ThinPlateSplineInterpolator : public InterpolatorBase<ThinPlateSplineInterpolator<Real>>
{
  public:
    using BASE = InterpolatorBase<ThinPlateSplineInterpolator<Real>>;
    using VectorType = typename BASE::VectorType;
    using MapType = typename BASE::MapType;

    // types used to view data as 2D coordinates
    using MatrixType    = typename BASE::MatrixType;
    using _2DVectorView = typename BASE::_2DVectorView;
    using _2DMatrixView = typename BASE::_2DMatrixView;

    // types used for 2x2 matrix algebra
    using Matrix22 = Eigen::Matrix<Real,2,2 >;
    using Matrix22Array = Eigen::Array< Matrix22, Eigen::Dynamic, Eigen::Dynamic >;
    using ColVector2 = Eigen::Matrix<Real,2,1 >;
    using RowVector2 = Eigen::Matrix<Real,1,2 >;

    ThinPlateSplineInterpolator() = default;
    ThinPlateSplineInterpolator(const ThinPlateSplineInterpolator& interp) = default;

    template<typename I>
    ThinPlateSplineInterpolator( I n, Real *x, Real *y, Real *z, bool deep_copy = true )
    {this->setData(n,x,y,z,deep_copy);}
    template<typename X, typename Y, typename Z>
    ThinPlateSplineInterpolator( X &x, Y &y, Z &z, bool deep_copy = true )
    {this->setData(x,y,z,deep_copy);}

    Real operator()( Real x, Real y ) const;


  protected:
    using BASE::xv;
    using BASE::yv;
    using BASE::zv;
    using BASE::X;
    using BASE::Y;
    using BASE::Z;
    
    MatrixType a, b;

    Real G(Real x, Real y, Real xi, Real yi) const;

    void setupInterpolator();
    friend BASE;




};

template<class Real>
void
ThinPlateSplineInterpolator<Real>::setupInterpolator()
{
  a = MatrixType(this->xv->rows(),1);
  b = MatrixType(this->xv->rows(),3);

  MatrixType M = MatrixType(this->xv->rows(),this->xv->rows());
  MatrixType N = MatrixType(this->xv->rows(),3);

  // rows
  for( int i = 0; i < this->xv->rows(); i++ )
  {
    // N
    N(i,0) = 1;
    N(i,1) = this->xv->coeff(i);
    N(i,2) = this->yv->coeff(i);

    for( int j = 0; j < this->xv->rows(); j++ )
    {
      M(i,j) = G(this->xv->coeff(i), this->yv->coeff(i), this->xv->coeff(j), this->yv->coeff(j) );
    }


  }

  MatrixType Minv = M.inverse();
  MatrixType Ntra = N.transpose();

  b = (Ntra*Minv*N).inverse()*Ntra*Minv*this->zd;
  a = Minv*(this->zd - N*b);

  return;

}

template<class Real>
Real
ThinPlateSplineInterpolator<Real>::G(Real x1, Real y1, Real x2, Real y2) const
{
  if( x1 == x2 && y1 == y2 )
    return 0;

  Real r = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );

  return r*r*log(r);
}


template<class Real>
Real
ThinPlateSplineInterpolator<Real>::operator()( Real x, Real y ) const
{
  BASE::checkData();
  
  // no extrapolation...
  if( x < this->xv->minCoeff()
   || x > this->xv->maxCoeff()
   || y < this->yv->minCoeff()
   || y > this->yv->maxCoeff() )
  {
    return 0;
  }
  
  MatrixType Gx( 1, this->xv->rows() );
  for(int i = 0; i < this->xv->rows(); i++)
    Gx(i) = G( x, y, this->xv->coeff(i), this->yv->coeff(i) );

  Real f = (Gx*a)(0,0) + b(0) + b(1)*x + b(2)*y;

  return f;
}





}

#endif // include protector
