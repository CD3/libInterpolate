#ifndef Interpolators__1D_InterpolatorBase_hpp
#define Interpolators__1D_InterpolatorBase_hpp

/** @file InterpolatorBase.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 12/24/16
  */

#include <memory>
#include <algorithm>
#include <type_traits>
#include <vector>
#include <Eigen/Dense>
#include<boost/range/algorithm.hpp>


namespace _1D {

/** @class 
  * @brief A base class that provides some useful functions for interpolating data.
  * @author C.D. Clark III
  *
  * This class provides eigen matrices to store the data that is interpolated, a few setData methods
  * to populate the data.
  */

template<typename T>
struct RealTypeOf { };
template<template<typename> class T,typename R>
struct RealTypeOf<T<R>> { using type = R; };

template<class Derived, typename Real = typename RealTypeOf<Derived>::type>
class InterpolatorBase
{

  public:
    using VectorType =  Eigen::Matrix<Real,Eigen::Dynamic,1>;
    using MapType =  Eigen::Map<const VectorType>;


  private:
    // we don't want base classes accessing these directly,
    // they should use xView and yView instead.
    std::vector<Real> xData, yData;        ///< storage for interpolated data

  protected:
    std::unique_ptr<MapType> xView, yView; ///< eigen matrix view of the data

  private:
    // making constructors private and the derived class a friend
    // helps to make sure that the derived class actually
    // passes itself as the template argument.
    friend Derived;

    InterpolatorBase(){ }

    // copy constructor
    // we only want to initialize x and y views if the object
    // we are copying from is initialized.
    InterpolatorBase(const InterpolatorBase& rhs)
    :xData(rhs.xData)
    ,yData(rhs.yData)
    ,xView( rhs.xView ? new MapType( xData.data(), xData.size() ) : nullptr )
    ,yView( rhs.xView ? new MapType( yData.data(), yData.size() ) : nullptr )
    { 
    }

    // we use the copy-swap idiom for copy assignment
    friend void swap( InterpolatorBase& lhs, InterpolatorBase& rhs)
    {
      std::swap( lhs.xData, rhs.xData );
      std::swap( lhs.yData, rhs.yData );
      std::swap( lhs.xView, rhs.xView );
      std::swap( lhs.yView, rhs.yView );
    }

    InterpolatorBase& operator=(InterpolatorBase rhs)
    {
      swap(*this,rhs);
      return *this;
    }


  public:
    // methods to get the data
    const std::vector<Real>& getXData() const { return xData; }
    const std::vector<Real>& getYData() const { return yData; }
    std::vector<Real> getXData() { return xData; }
    std::vector<Real> getYData() { return yData; }



    /**
     * Set references to the interpolated data to memory outside the interpolator.
     * This is potentially *unsafe*. The interpolator will not take ownership of the memory,
     * the caller will still be required to free it when it is no longer needed. Freeing memory
     * before calling the interpolator may cause uninitialized memory access and is undefined
     * behavior.
     *
     * Modifying the data after setting the interpolator references is undefined behavior as some
     * interpolators do some setup work that will be negated when the data has changed.
     */
    template<typename I>
    void
    setUnsafeDataReference( I n, const Real *x, const Real *y)
    {
      this->xView.reset( new MapType( x, n ) );
      this->yView.reset( new MapType( y, n ) );

      this->callSetupInterpolator<Derived>();
    }


    /**
     * Set the data that will be interpolated from a pair of std::vectors.
     */
    template<typename X, typename Y>
    auto
    setData( const X &x, const Y &y)
    -> decltype(x.size(),x.data(),y.data(),void())
    {
      this->setData( x.size(), x.data(), y.data() );
    }

    /**
     * Set the data that will be interpolated from raw data pointers.
     *
     * The data pointed at by x and y will be *copied* into the interpolator.
     */
    template<typename I, typename X, typename Y>
    void
    setData( I n, const X *x, const Y *y)
    {
      this->setData(x,x+n,y,y+n);
    }

    template<typename XIter, typename YIter>
    auto
    setData( const XIter &x_begin, const XIter &x_end, const YIter &y_begin, const YIter &y_end)
    -> decltype(std::copy(x_begin,x_end,xData.begin()),std::copy(y_begin,y_end,yData.begin()),void())
    {
      xData.clear();
      yData.clear();
      xData.reserve(x_end-x_begin);
      yData.reserve(y_end-y_begin);
      std::copy( x_begin, x_end, std::back_inserter(xData) );
      std::copy( y_begin, y_end, std::back_inserter(yData) );

      this->setUnsafeDataReference( xData.size(), xData.data(), yData.data() );
    }



    /**
     * Given an x value, returns the index of the stored x data
     * that is just to the "left" of x such that X[i] < x < X[i+1]
     */
    int
    get_index_to_left_of(Real x) const
    {
      auto rng = std::make_pair( xView->data()+1, xView->data()+xView->size() );
      return boost::lower_bound( rng, x) - xView->data() - 1;
    }
    /**
     * Given an x value, returns the index of the stored x data
     * that is just to the "right" of x such that X[i-1] < x < X[i]
     */
    int
    get_index_to_right_of(Real x) const
    {
      auto rng = std::make_pair( xView->data()+1, xView->data()+xView->size() );
      return boost::lower_bound( rng, x) - xView->data();
    }


  protected:

    void checkData() const; ///< Check that data has been initialized and throw exception if not.

  private:
    // callSetupInterpolator will call a function named setupInterpolator in the derived class, if
    // it exists. this is just some template magic to detect if the derived class has implemented a
    // setupInterpolator function, and to call it if it does.

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


// implementations


template<class Derived, typename Real>
void
InterpolatorBase<Derived,Real>::checkData() const
{
  if(!this->xView || !this->yView)
    throw std::logic_error("Interpolator data is not initialized. Did you call setData()?");
  if(this->xView->size() == 0 || this->yView->size() == 0)
    throw std::logic_error("Interpolator data is zero size. Did you call setData() with non-zero sized vectors?");
}


}

#endif // include protector
