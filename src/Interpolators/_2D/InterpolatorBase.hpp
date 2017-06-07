#ifndef Interpolators__2D_InterpolatorBase_hpp
#define Interpolators__2D_InterpolatorBase_hpp

/** @file InterpolatorBase.hpp
  * @author C.D. Clark III
  * @date 12/24/16
  */

#include <memory>
#include "InterpolatorInterface.hpp"
#include "Utils/Indexing.hpp"


namespace _2D {

/** @class 
  * @brief A base class for 2D interpolators that provides default implementations of the interface.
  * @author C.D. Clark III
  *
  * This class provides an implementation for the setData method as well as adding a few additional
  * useful methods, including derivative and integral methods.
  */
template<class Real>
class InterpolatorBase : public InterpolatorInterface<Real>
{
  public:
    typedef Eigen::Matrix<Real,Eigen::Dynamic,1> VectorType;
    typedef Eigen::Map<VectorType> MapType;
    typedef Eigen::Matrix<Real,2,1> GradientType;

    // methods required by interface
    virtual void setData( size_t _n, Real *x, Real *y, Real *z, bool deep_copy = true );

    // methods to get interpolated data
    virtual std::vector<Real> getXData() const { return std::vector<Real>(&xd(0),&xd(0)+xd.size()); }
    virtual std::vector<Real> getYData() const { return std::vector<Real>(&yd(0),&yd(0)+yd.size()); }
    virtual std::vector<Real> getZData() const { return std::vector<Real>(&zd(0),&zd(0)+zd.size()); }
    

    // additional methods
    virtual GradientType gradient( Real x, Real y ) const;
    virtual Real integral(   Real a, Real b, Real c, Real d ) const;
    virtual void setData( std::vector<Real> &x, std::vector<Real> &y, std::vector<Real> &z, bool deep_copy = true );
    virtual void setData( VectorType  &x, VectorType &y, VectorType &z, bool deep_copy = true );

  protected:
    VectorType xd, yd, zd;               // data
    std::shared_ptr<MapType> xv, yv, zv; // map view of the data
    void checkData() const; ///< Check that data has been initialized and throw exception if not.



};

template<class Real>
void
InterpolatorBase<Real>::checkData() const
{
  if(!this->xv || !this->yv || !this->zv)
    throw std::logic_error("Interpolator data is not initialized. Did you call setData()?");
}

/**
 * Reads data from x, y, and z arrays. Arrays should all be the same length with each element corresponding to a data point.
 * Basically, x[i], y[i], and z[i] should correspond to the first, second, and third columns of a gnuplot file.
 */
template<class Real>
void
InterpolatorBase<Real>::setData( size_t n, Real *x, Real *y, Real *z, bool deep_copy )
{
  Real *xp, *yp, *zp;
  if( deep_copy )
  {
    xd = Eigen::Map<VectorType>( x, n );
    yd = Eigen::Map<VectorType>( y, n );
    zd = Eigen::Map<VectorType>( z, n );
    xp = &xd(0);
    yp = &yd(0);
    zp = &zd(0);
  }
  else
  {
    xp = x;
    yp = y;
    zp = z;
  }
  this->xv.reset( new Eigen::Map<VectorType>( xp, n ) );
  this->yv.reset( new Eigen::Map<VectorType>( yp, n ) );
  this->zv.reset( new Eigen::Map<VectorType>( zp, n ) );
}

template<class Real>
void
InterpolatorBase<Real>::setData( std::vector<Real> &x, std::vector<Real> &y, std::vector<Real> &z, bool deep_copy )
{
  this->setData( x.size(), x.data(), y.data(), z.data(), deep_copy );
}

template<class Real>
void
InterpolatorBase<Real>::setData( VectorType  &x, VectorType &y, VectorType &z, bool deep_copy )
{
  this->setData( x.size(), &x(0), &y(0), &z(0), deep_copy );
}


template<class Real>
auto
InterpolatorBase<Real>::gradient( Real x, Real y ) const -> GradientType
{
  if( xv->size() < 1 )
    return {0,0};

  // simple Finite-Difference approximation
  // NOTE: consecutive values in xv and yv may be the same.
  //       so be careful here
  Real dx = (xv->maxCoeff() - xv->minCoeff())/xv->cols();
  Real dy = (yv->maxCoeff() - yv->minCoeff())/yv->cols();
  Real dzx = this->operator()( x + dx/2, y        ) - this->operator()( x - dx/2, y        );
  Real dzy = this->operator()( x       , y + dy/2 ) - this->operator()( x       , y - dy/2 );

  return {dzx/dx, dzy/dy};
}

template<class Real>
Real
InterpolatorBase<Real>::integral( Real a, Real b, Real c, Real d ) const
{
  if( xv->size() < 1 )
    return 0;

  int sign = 1;
  // support for flipped integration limits
  if( a > b )
  {
    std::swap( a, b );
    sign *= -1;
  }
  if( c > d )
  {
    std::swap( c, d );
    sign *= -1;
  }

  // simple Trapezoid sum
  Real dx = (*xv)(1) - (*xv)(0);
  int Nx = 2*(b-a)/dx;
  Nx = std::max(2,Nx); // make sure N is at least 2
  dx = (b-a)/(Nx-1);

  Real dy = (*yv)(1) - (*yv)(0);
  int Ny = 2*(d-c)/dy;
  Ny = std::max(2,Ny);
  dy = (d-c)/(Ny-1);

  Real sum = 0;
  for(int j = 0; j < Ny-1; j++)
  {
    Real sum1 = 0, sum2 = 0;
    for(int i = 0; i < Nx-1; i++)
      sum1 += this->operator()(a+i*dx,c+j*dy) + this->operator()(a+(i+1)*dx,c+j*dy);
    for(int i = 0; i < Nx-1; i++)
      sum2 += this->operator()(a+i*dx,c+(j+1)*dy) + this->operator()(a+(i+1)*dx,c+(j+1)*dy);

    sum1 *= 0.5*dx; // slight optimization
    sum2 *= 0.5*dx; // only works if dx is constant

    sum += sum1 + sum2;
  }
  sum *= 0.5*dy; // again, requires that dy is the same for all element

  return sign*sum;
}

}

#endif // include protector
