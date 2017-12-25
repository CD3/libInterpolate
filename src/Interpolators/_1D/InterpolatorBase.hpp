#ifndef Interpolators__1D_InterpolatorBase_hpp
#define Interpolators__1D_InterpolatorBase_hpp

/** @file InterpolatorBase.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 12/24/16
  */

#include <memory>
#include <type_traits>
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

template<typename T>
struct RealTypeOf { using type = double; };
template<template<typename> typename T,typename R>
struct RealTypeOf<T<R>> { using type = R; };

template<class Derived, typename Real = typename RealTypeOf<Derived>::type>
class InterpolatorBase
{

  public:
    using VectorType =  Eigen::Matrix<Real,Eigen::Dynamic,1>;
    using MapType =  Eigen::Map<const VectorType>;

    InterpolatorBase(const InterpolatorBase& interp) = default;

    template<typename I>
    void setData( I n, const Real *x, const Real *y, bool deep_copy = true );

    template<typename X, typename Y>
    // this template is ambiguous with the pointer template above,
    // wo we want to disable it for pointers
    typename std::enable_if<!std::is_pointer<Y>::value>::type
    setData( const X &x, const Y &y, bool deep_copy = true );

    // methods to get the data
    std::vector<Real> getXData() const { return std::vector<Real>(&xd(0),&xd(0)+xd.size()); }
    std::vector<Real> getYData() const { return std::vector<Real>(&yd(0),&yd(0)+yd.size()); }


  protected:
    VectorType xd, yd;               // data
    std::shared_ptr<MapType> xv, yv; // map view of the data
    void checkData() const;         ///< Check that data has been initialized and throw exception if not.

  private:
    // this helps to make sure that the derived class actually
    // passes itself as the template argument.
    friend Derived;
    InterpolatorBase(){}






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
  if(!this->xv || !this->yv)
    throw std::logic_error("Interpolator data is not initialized. Did you call setData()?");
}

template<class Derived, typename Real>
template<typename I>
void
InterpolatorBase<Derived,Real>::setData( I n, const Real *x, const Real *y, bool deep_copy )
{
  const Real *xp, *yp;
  if( deep_copy )
  {
    xd = MapType( x, n );
    yd = MapType( y, n );
    xp = &xd(0);
    yp = &yd(0);
  }
  else
  {
    xp = x;
    yp = y;
  }
  this->xv.reset( new MapType( xp, n ) );
  this->yv.reset( new MapType( yp, n ) );

  this->callSetupInterpolator<Derived>();
}

template<class Derived, typename Real>
template<typename X, typename Y>
typename std::enable_if<!std::is_pointer<Y>::value>::type
InterpolatorBase<Derived,Real>::setData( const X &x, const Y &y, bool deep_copy )
{
  this->setData( x.size(), x.data(), y.data(), deep_copy );
}


}

#endif // include protector
