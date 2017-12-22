#ifndef Interpolators__2D_InterpolatorBase_hpp
#define Interpolators__2D_InterpolatorBase_hpp

/** @file InterpolatorBase.hpp
  * @author C.D. Clark III
  * @date 12/24/16
  */

#include <memory>
#include <Eigen/Dense>
#include "Utils/Indexing.hpp"


namespace _2D {

template<typename T>
struct RealTypeOf { using type = double; };
template<template<typename> typename T,typename R>
struct RealTypeOf<T<R>> { using type = R; };

/** @class 
  * @brief A base class for 2D interpolators that provides default implementations of the interface.
  * @author C.D. Clark III
  *
  * This class provides an implementation for the setData method as well as adding a few additional
  * useful methods, including derivative and integral methods.
  */
template<class Derived, typename Real = typename RealTypeOf<Derived>::type>
class InterpolatorBase
{
  public:
    using VectorType = Eigen::Matrix<Real,Eigen::Dynamic,1>;
    using MapType = Eigen::Map<VectorType>;

    template<typename I>
    void setData( I n, Real *x, Real *y, Real *z, bool deep_copy = true );
    template<typename X, typename Y, typename Z>
    void setData( X &x, Y &y, Z &z, bool deep_copy = true );

    std::vector<Real> getXData() const { return std::vector<Real>(&xd(0),&xd(0)+xd.size()); }
    std::vector<Real> getYData() const { return std::vector<Real>(&yd(0),&yd(0)+yd.size()); }
    std::vector<Real> getZData() const { return std::vector<Real>(&zd(0),&zd(0)+zd.size()); }


  protected:
    VectorType xd, yd, zd;               // data
    std::shared_ptr<MapType> xv, yv, zv; // map view of the data
    void checkData() const; ///< Check that data has been initialized and throw exception if not.

  private:
    // this helps to make sure that the derived class actually
    // passes itself as the template argument.
    InterpolatorBase(){ }
    friend Derived;

    // below is some template magic to detect if the derived class has implemented a
    // setupInterpolator function.

    template<typename T>
    struct has_setupInterpolator
    {
      private:
        typedef std::true_type yes;
        typedef std::false_type no;

        template<typename U> static auto test(int) -> decltype(std::declval<U>().setupInterpolator(),yes());
        template<typename  > static no   test(...);

      public:
        static constexpr bool value = std::is_same<decltype(test<T>(0)),yes>::value;
    };
    
    template<typename T>
    typename std::enable_if<has_setupInterpolator<T>::value>::type
    callSetupInterpolator(){ static_cast<Derived*>(this)->setupInterpolator(); }

    template<typename T>
    typename std::enable_if<!has_setupInterpolator<T>::value>::type
    callSetupInterpolator(){ }


};

template<class Derived, typename Real>
void
InterpolatorBase<Derived,Real>::checkData() const
{
  if(!this->xv || !this->yv || !this->zv)
    throw std::logic_error("Interpolator data is not initialized. Did you call setData()?");
}

/**
 * Reads data from x, y, and z arrays. Arrays should all be the same length with each element corresponding to a data point.
 * Basically, x[i], y[i], and z[i] should correspond to the first, second, and third columns of a gnuplot file.
 */
template<class Derived, typename Real>
template<typename I>
void
InterpolatorBase<Derived,Real>::setData( I n, Real *x, Real *y, Real *z, bool deep_copy )
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

  this->callSetupInterpolator<Derived>();
}

template<class Derived, typename Real>
template<typename X, typename Y, typename Z>
void
InterpolatorBase<Derived,Real>::setData( X &x, Y &y, Z &z, bool deep_copy )
{
  this->setData( x.size(), x.data(), y.data(), z.data(), deep_copy );
}

}

#endif // include protector
