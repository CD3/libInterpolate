#ifndef Interpolators__1D_InterpolatorBase_hpp
#define Interpolators__1D_InterpolatorBase_hpp

/** @file InterpolatorBase.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 12/24/16
  */

#include <memory>
#include <Eigen/Dense>
#include "Utils/Indexing.hpp"


namespace _1D {

/** @class 
  * @brief A base class that provides some useful functions for interpolating data.
  * @author C.D. Clark III
  *
  * This class provides eigen matrices to store the data that is interpolated, a few setData methods
  * to populate the data.
  */
template<class Real>
class InterpolatorBase
{
  public:
    typedef Eigen::Matrix<Real,Eigen::Dynamic,1> VectorType;
    typedef Eigen::Map<VectorType> MapType;

    template<typename I>
    void setData( I n, Real *x, Real *y, bool deep_copy = true );

    template<typename X, typename Y>
    void setData( X &x, Y &y, bool deep_copy = true );

    // methods to get the data
    std::vector<Real> getXData() const { return std::vector<Real>(&xd(0),&xd(0)+xd.size()); }
    std::vector<Real> getYData() const { return std::vector<Real>(&yd(0),&yd(0)+yd.size()); }


  protected:
    VectorType xd, yd;               // data
    std::shared_ptr<MapType> xv, yv; // map view of the data
    void checkData() const; ///< Check that data has been initialized and throw exception if not.


};

template<class Real>
void
InterpolatorBase<Real>::checkData() const
{
  if(!this->xv || !this->yv)
    throw std::logic_error("Interpolator data is not initialized. Did you call setData()?");
}

template<class Real>
template<typename I>
void
InterpolatorBase<Real>::setData( I n, Real *x, Real *y, bool deep_copy )
{
  Real *xp, *yp;
  if( deep_copy )
  {
    xd = Eigen::Map<VectorType>( x, n );
    yd = Eigen::Map<VectorType>( y, n );
    xp = &xd(0);
    yp = &yd(0);
  }
  else
  {
    xp = x;
    yp = y;
  }
  this->xv.reset( new Eigen::Map<VectorType>( xp, n ) );
  this->yv.reset( new Eigen::Map<VectorType>( yp, n ) );
}

template<class Real>
template<typename X, typename Y>
void
InterpolatorBase<Real>::setData( X &x, Y &y, bool deep_copy )
{
  this->setData( x.size(), x.data(), y.data(), deep_copy );
}


}

#endif // include protector
