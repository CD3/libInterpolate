#ifndef Interpolators__2D_InterpolatorBase_hpp
#define Interpolators__2D_InterpolatorBase_hpp

/** @file InterpolatorBase.hpp
  * @author C.D. Clark III
  * @date 12/24/16
  */

#include <memory>
#include <vector>
#include <Eigen/Dense>
#include<boost/range/algorithm.hpp>
#include <boost/range/adaptor/strided.hpp>


namespace _2D {

template<typename T>
struct RealTypeOf { using type = double; };
template<template<typename> class T,typename R>
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
    using MapType = Eigen::Map<const VectorType>;
    // types used to view data as 2D coordinates
    using MatrixType = typename Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic>;
    using _2DVectorView = Eigen::Map<const VectorType,Eigen::Unaligned,Eigen::InnerStride<Eigen::Dynamic>>;
    using _2DMatrixView = Eigen::Map<const MatrixType,Eigen::Unaligned,Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>>;

  private:
    std::vector<Real> xData, yData, zData;        // data storage

  protected:
    std::unique_ptr<MapType> xView, yView, zView; // map view of the data

    // these maps are used to view the x,y,z data as two coordinate vectors and a function matrix, instead of three vectors.
    std::unique_ptr<_2DVectorView> X,Y;
    std::unique_ptr<_2DMatrixView> Z;

  private:
    // making constructors private and the derived class a friend
    // helps to make sure that the derived class actually
    // passes itself as the template argument.
    friend Derived;

    InterpolatorBase(){ }

    InterpolatorBase( const InterpolatorBase& rhs )
    :xData(rhs.xData)
    ,yData(rhs.yData)
    ,zData(rhs.zData)
    ,xView( rhs.xView ? new MapType( xData.data(), xData.size() ) : nullptr )
    ,yView( rhs.xView ? new MapType( yData.data(), yData.size() ) : nullptr )
    ,zView( rhs.zView ? new MapType( zData.data(), zData.size() ) : nullptr )
    {
      this->setup2DDataViews();
    }

    // we use the copy-swap idiom for copy assignment
    friend void swap( InterpolatorBase& lhs, InterpolatorBase& rhs)
    {
      std::swap( lhs.xData, rhs.xData );
      std::swap( lhs.yData, rhs.yData );
      std::swap( lhs.zData, rhs.zData );
      std::swap( lhs.xView, rhs.xView );
      std::swap( lhs.yView, rhs.yView );
      std::swap( lhs.zView, rhs.zView );
      std::swap( lhs.X, rhs.X);
      std::swap( lhs.Y, rhs.Y);
      std::swap( lhs.Z, rhs.Z);
    }

    InterpolatorBase& operator=(InterpolatorBase rhs)
    {
      swap(*this,rhs);
      return *this;
    }

  public:

    const std::vector<Real>& getXData() const { return xData; }
    const std::vector<Real>& getYData() const { return yData; }
    const std::vector<Real>& getZData() const { return zData; }
    std::vector<Real> getXData() { return xData; }
    std::vector<Real> getYData() { return yData; }
    std::vector<Real> getZData() { return zData; }


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
    setUnsafeDataReference( I n, const Real *x, const Real *y, const Real *z)
    {
      this->xView.reset( new MapType( x, n ) );
      this->yView.reset( new MapType( y, n ) );
      this->zView.reset( new MapType( z, n ) );

      this->setup2DDataViews();
      this->callSetupInterpolator<Derived>();
    }


    /**
     * Set the data taht will be interpolated from a set of iterators.
     */
    template<typename XIter, typename YIter, typename ZIter>
    auto
    setData( const XIter &x_begin, const XIter &x_end, const YIter &y_begin, const YIter &y_end, const ZIter &z_begin, const ZIter &z_end)
    -> decltype(std::copy(x_begin,x_end,xData.begin()),std::copy(y_begin,y_end,yData.begin()),std::copy(z_begin,z_end,zData.begin()),void())
    {
      xData.clear();
      yData.clear();
      zData.clear();
      xData.reserve(x_end-x_begin);
      yData.reserve(y_end-y_begin);
      zData.reserve(z_end-z_begin);
      std::copy( x_begin, x_end, std::back_inserter(xData) );
      std::copy( y_begin, y_end, std::back_inserter(yData) );
      std::copy( z_begin, z_end, std::back_inserter(zData) );

      this->setUnsafeDataReference( xData.size(), xData.data(), yData.data(), zData.data() );
    }


    /**
     * Reads data from x, y, and z arrays. Arrays should all be the same length with each element corresponding to a data point.
     * Basically, x[i], y[i], and z[i] should correspond to the first, second, and third columns of a gnuplot file.
     */
    template<typename I, typename X, typename Y, typename Z>
    void
    setData( I n, const X *x, const Y *y, const Z *z )
    {
      this->setData(x,x+n,y,y+n,z,z+n);
    }

    /**
     * Set the data to be interpolated from vector-like containers that provide a .data() and .size()
     * methods.
     */
    template<typename XT, typename YT, typename ZT>
    auto
    setData( const XT &x, const YT &y, const ZT &z )
    -> decltype( x.size(), x.data(), y.data(), z.data(), void())
    {
      return this->setData( x.size(), x.data(), y.data(), z.data() );
    }






    /**
     * Given an x value, returns the index i of the stored x data X
     * that is just to the "left" of x such that X[i] < x < X[i+1]
     */
    int
    get_x_index_to_left_of(Real x) const
    {
      // NOTE: X data is strided.
      auto xrng = std::make_pair( X->data(), X->data()+X->size()*X->innerStride() ) | boost::adaptors::strided(X->innerStride());
      return boost::lower_bound( xrng, x) - boost::begin(xrng) - 1;
    }
    /**
     * Given an x value, returns the index i of the stored x data X
     * that is just to the "right" of x such that X[i-1] < x < X[i]
     */
    int
    get_x_index_to_right_of(Real x) const
    {
      return this->get_x_index_to_right_of(x)+1;
    }

    /**
     * Given an y value, returns the index i of the stored y data Y
     * that is just "below" y such that Y[i] < y < Y[i+1]
     */
    int
    get_y_index_below(Real y) const
    {
    // NOTE: Y data is NOT strided
      auto yrng = std::make_pair( Y->data(), Y->data()+Y->size() );
      return boost::lower_bound( yrng, y) - boost::begin(yrng) - 1;
    }
    /**
     * Given an y value, returns the index i of the stored y data y
     * that is just "above" y such that Y[i-1] < y < Y[i]
     */
    int
    get_y_index_above(Real y) const
    {
      return this->get_y_index_below(y)+1;
    }

  protected:
    void checkData() const; ///< Check that data has been initialized and throw exception if not.
    void setup2DDataViews(); ///< Setups up 2D views of 1D data arrays

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

template<class Derived, typename Real>
void
InterpolatorBase<Derived,Real>::checkData() const
{
  if(!this->xView || !this->yView || !this->zView)
    throw std::logic_error("Interpolator data is not initialized. Did you call setData()?");
  if(this->xView->size() == 0 || this->yView->size() == 0 || this->zView->size() == 0)
    throw std::logic_error("Interpolator data is zero size. Did you call setData() with non-zero sized vectors?");
}

template<class Derived, typename Real>
void
InterpolatorBase<Derived,Real>::setup2DDataViews()
{
  // make sure data has actually been initialized
  if( !xView || !yView || !zView )
    return;


  // setup 2D view of the data
  // We need to figure out what the x and y dimensions are.
  int N = zView->size();
  int Nx = 0, Ny = 0;
  // Ny will be the number of elements that have the same x coordinate
  Real xlast = (*xView)(0);
  while( Ny < N-1 && fabs((*xView)(Ny)-xlast) < 1e-40 )
    Ny++;
  Nx = N/Ny;

  // consecutive values in the x data are separated by Ny, so this is the inner stride for X
  X.reset( new _2DVectorView( xView->data(), Nx, Eigen::InnerStride<Eigen::Dynamic>(Ny) ) );

  // consecutive values in the y data are next to each other, so the stride is just 1
  Y.reset( new _2DVectorView( yView->data(), Ny, Eigen::InnerStride<Eigen::Dynamic>(1) ) );

  // Eigen defaults to COLUMN MAJOR
  // consecutive elements in a column are separated by Ny (this is the inner stride)
  // consecutive elements in a row are located next to each other (this is the outer stride)
  // Stride object takes outer,inner as arguments.
  Z.reset( new _2DMatrixView( zView->data(), Nx, Ny, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(1,Ny) ) );
}


}

#endif // include protector
