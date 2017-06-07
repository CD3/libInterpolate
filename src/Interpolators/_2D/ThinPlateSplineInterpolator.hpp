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
class ThinPlateSplineInterpolator : public InterpolatorBase<Real>
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
    
    MatrixType a, b;

    Real G(Real x, Real y, Real xi, Real yi) const;
    void calcCoefficients();




};

template<class Real>
void
ThinPlateSplineInterpolator<Real>::setData( size_t n, Real *x, Real *y, Real *z, bool deep_copy )
{
  InterpolatorBase<Real>::setData(n,x,y,z,deep_copy);
  calcCoefficients();
}

template<class Real>
void
ThinPlateSplineInterpolator<Real>::setData( std::vector<Real> &x, std::vector<Real> &y, std::vector<Real> &z, bool deep_copy )
{
  InterpolatorBase<Real>::setData(x,y,z,deep_copy);
  calcCoefficients();
}

template<class Real>
void
ThinPlateSplineInterpolator<Real>::setData( VectorType  &x, VectorType &y, VectorType &z, bool deep_copy )
{
  InterpolatorBase<Real>::setData(x,y,z,deep_copy);
  calcCoefficients();
}

template<class Real>
void
ThinPlateSplineInterpolator<Real>::calcCoefficients()
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
  auto Minv = M.inverse().eval();

  // TODO: could optimize for computations, but would require more memory. i.e. store the
  // inverse of M and transpose of N in temporary matrices.
  b = (N.transpose()*M.inverse()*N).inverse()*N.transpose()*M.inverse()*this->zd;
  a = M.inverse()*(this->zd - N*b);

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
  InterpolatorBase<Real>::checkData();
  
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
