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
    using MapType = Eigen::Map<const VectorType>;
    // types used to view data as 2D coordinates
    using MatrixType = typename Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic>;
    using _2DVectorView = Eigen::Map<const VectorType,Eigen::Unaligned,Eigen::InnerStride<Eigen::Dynamic>>;
    using _2DMatrixView = Eigen::Map<const MatrixType,Eigen::Unaligned,Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>>;

    InterpolatorBase( const InterpolatorBase& interp ) = default;

    template<typename I>
    void setData( I n, const Real *x, const Real *y, const Real *z, bool deep_copy = true );
    template<typename XT, typename YT, typename ZT>
    // this template is ambiguous with the pointer template above,
    // wo we want to disable it for pointers
    typename std::enable_if<!std::is_pointer<YT>::value>::type
    setData( const XT &x, const YT &y, const ZT &z, bool deep_copy = true );

    std::vector<Real> getXData() const { return std::vector<Real>(&xd(0),&xd(0)+xd.size()); }
    std::vector<Real> getYData() const { return std::vector<Real>(&yd(0),&yd(0)+yd.size()); }
    std::vector<Real> getZData() const { return std::vector<Real>(&zd(0),&zd(0)+zd.size()); }


  protected:
    VectorType xd, yd, zd;               // data
    std::shared_ptr<MapType> xv, yv, zv; // map view of the data

    // these maps are used to view the x,y,z data as two coordinate vectors and a function matrix, instead of three vectors.
    std::shared_ptr<_2DVectorView> X,Y;
    std::shared_ptr<_2DMatrixView> Z;

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
InterpolatorBase<Derived,Real>::setData( I n, const Real *x, const Real *y, const Real *z, bool deep_copy )
{
  const Real *xp, *yp, *zp;
  if( deep_copy )
  {
    xd = MapType( x, n );
    yd = MapType( y, n );
    zd = MapType( z, n );
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
  this->xv.reset( new MapType( xp, n ) );
  this->yv.reset( new MapType( yp, n ) );
  this->zv.reset( new MapType( zp, n ) );


  // setup 2D view of the data
  // We need to figure out what the x and y dimensions are.
  int N = zv->size();
  int Nx = 0, Ny = 0;
  // Ny will be the number of elements that have the same x coordinate
  Real xlast = (*xv)(0);
  while( Ny < N-1 && fabs((*xv)(Ny)-xlast) < 1e-40 )
    Ny++;
  Nx = N/Ny;

  // consecutive values in the x data are separated by Ny, so this is the inner stride for X
  X.reset( new _2DVectorView( xv->data(), Nx, Eigen::InnerStride<Eigen::Dynamic>(Ny) ) );

  // consecutive values in the y data are next to each other, so the stride is just 1
  Y.reset( new _2DVectorView( yv->data(), Ny, Eigen::InnerStride<Eigen::Dynamic>(1) ) );

  // Eigen defaults to COLUMN MAJOR
  // consecutive elements in a column are separated by Ny (this is the inner stride)
  // consecutive elements in a row are located next to each other (this is the outer stride)
  // Stride object takes outer,inner as arguments.
  Z.reset( new _2DMatrixView( zv->data(), Nx, Ny, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(1,Ny) ) );





  this->callSetupInterpolator<Derived>();
}

template<class Derived, typename Real>
template<typename XT, typename YT, typename ZT>
typename std::enable_if<!std::is_pointer<YT>::value>::type
InterpolatorBase<Derived,Real>::setData( const XT &x, const YT &y, const ZT &z, bool deep_copy )
{
  this->setData( x.size(), x.data(), y.data(), z.data(), deep_copy );
}

}

#endif // include protector
