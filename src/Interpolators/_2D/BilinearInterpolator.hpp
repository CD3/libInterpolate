#ifndef Interpolators__2D_BilinearInterpolator_hpp
#define Interpolators__2D_BilinearInterpolator_hpp

/** @file BilinearInterpolator.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 12/27/16
  */

#include "InterpolatorBase.hpp"
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/range/adaptor/strided.hpp>

/** @class 
  * @brief Linear interpolation for 2D functions.
  * @author C.D. Clark III, Aaron Hoffman
  */

namespace _2D {

template<class Real>
class BilinearInterpolator : public InterpolatorBase<BilinearInterpolator<Real>>
{
  public:
    using BASE = InterpolatorBase<BilinearInterpolator<Real>>;
    using VectorType = typename BASE::VectorType;
    using MapType = typename BASE::MapType; 
    // types used to view data as 2D coordinates
    using MatrixType = typename Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic>;
    using _2DVectorView = Eigen::Map<VectorType,Eigen::Unaligned,Eigen::InnerStride<Eigen::Dynamic>>;
    using _2DMatrixView = Eigen::Map<MatrixType,Eigen::Unaligned,Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>>;

    // types used for 2x2 matrix algebra
    using Matrix22 = Eigen::Matrix<Real,2,2 >;
    using Matrix22Array = Eigen::Array< Matrix22, Eigen::Dynamic, Eigen::Dynamic >;
    using ColVector2 = Eigen::Matrix<Real,2,1 >;
    using RowVector2 = Eigen::Matrix<Real,1,2 >;

    Real operator()( Real x, Real y ) const;

  protected:
    using BASE::xv;
    using BASE::yv;
    using BASE::zv;
    // these maps are used to view the x,y,z data as two coordinate vectors and a function matrix, instead of three vectors.
    std::shared_ptr<_2DVectorView> X,Y;
    std::shared_ptr<_2DMatrixView> Z;
    
    Matrix22Array Q; // naming convention used by wikipedia article (see Wikipedia https://en.wikipedia.org/wiki/Bilinear_interpolation)

    void setupInterpolator();
    friend BASE;

};

template<class Real>
void
BilinearInterpolator<Real>::setupInterpolator()
{

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
  X.reset( new _2DVectorView( &(*xv)(0), Nx, Eigen::InnerStride<Eigen::Dynamic>(Ny) ) );

  // consecutive values in the y data are next to each other, so the stride is just 1
  Y.reset( new _2DVectorView( &(*yv)(0), Ny, Eigen::InnerStride<Eigen::Dynamic>(1) ) );

  // Eigen defaults to COLUMN MAJOR
  // consecutive elements in a column are separated by Ny (this is the inner stride)
  // consecutive elements in a row are located next to each other (this is the outer stride)
  // Stride object takes outer,inner as arguments.
  Z.reset( new _2DMatrixView( &(*zv)(0), Nx, Ny, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(1,Ny) ) );

  // Interpolation will be done by multiplying the coordinates by coefficients.

  Q = Matrix22Array( X->size()-1, Y->size()-1 );

  // We are going to pre-compute the interpolation coefficients so
  // that we can interpolate quickly
  for(int i = 0; i < X->size() - 1; i++)
  {
    for( int j = 0; j < Y->size() - 1; j++)
    {
      Real tmp = ( ((*X)(i+1) - (*X)(i) )*( (*Y)(j+1) - (*Y)(j) ) );
      Q(i,j) = Z->block(i,j,2,2)/tmp;
    }
  }
}

template<class Real>
Real
BilinearInterpolator<Real>::operator()( Real x, Real y ) const
{
  BASE::checkData();
  
  // no extrapolation...
  if( x < (*X)(0)
   || x > (*X)(X->size()-1)
   || y < (*Y)(0)
   || y > (*Y)(Y->size()-1) )
  {
    return 0;
  }
  // find the x index that is just to the LEFT of x
  //int i  = Utils::index__last_lt( x, *X, X->size() );
  // NOTE: X data is strided.
  auto xrng = std::make_pair( X->data(), X->data()+X->size()*X->innerStride() ) | boost::adaptors::strided(X->innerStride());
  int i = boost::lower_bound( xrng, x) - boost::begin(xrng) - 1;
  if(i < 0)
    i = 0;

  // find the y index that is just BELOW y
  //int j  = Utils::index__last_lt( y, *Y, Y->size() );
  // NOTE: Y data is NOT strided
  auto yrng = std::make_pair( Y->data(), Y->data()+Y->size() );
  int j = boost::lower_bound( yrng, y) - boost::begin(yrng) - 1;
  if(j < 0)
    j = 0;
  


  // now, create the coordinate vectors (see Wikipedia https://en.wikipedia.org/wiki/Bilinear_interpolation)
  RowVector2 vx;
  ColVector2 vy;
  vx << ((*X)(i+1) - x), (x - (*X)(i));
  vy << ((*Y)(j+1) - y), (y - (*Y)(j));

  // interpolation is just x*Q*y

  return vx*Q(i,j)*vy;
}



}

#endif // include protector
