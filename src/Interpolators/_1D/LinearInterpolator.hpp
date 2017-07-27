#ifndef Interpolators__1D_LinearInterpolator_hpp
#define Interpolators__1D_LinearInterpolator_hpp

#include "InterpolatorBase.hpp"

#include<boost/range/algorithm.hpp>

namespace _1D {


/** @class 
  * @brief A simple linear interpolator.
  * @author C.D. Clark III
  *
  * This class does *not* do extrapolation.
  */
template<class Real>
class LinearInterpolator : public InterpolatorBase<Real>
{

  public:
    typedef typename InterpolatorBase<Real>::VectorType VectorType;
    typedef typename InterpolatorBase<Real>::MapType MapType;

    // function required by interface
    virtual Real operator()( Real x ) const {return call(x);}
    virtual Real derivative( Real x ) const {return derivative_imp(x);}
    // virtual Real integral( Real a, Real b) const {return integral_imp(a,b);}
    // use the base class operator to interpolate multiple points at once.
    using InterpolatorBase<Real>::operator();

    // non-virtual implementations
    // we want to implement the interpolation as a non-virtual method
    // in case somebody thinks that will be too slow. the virtual methods
    // will then just call these (they are already slow, right?).
    Real call(Real x ) const; // can't really create a function called oeprator()_imp
    Real derivative_imp( Real x ) const;


};

template<class Real>
Real
LinearInterpolator<Real>::call( Real x ) const
{
  InterpolatorBase<Real>::checkData();

  const VectorType &X = *(this->xv);
  const VectorType &Y = *(this->yv);

  // find the index that is just to the left of the x
  //int i = Utils::index__last_lt( x, X, X.size(), 0);
  auto rng = std::make_pair( X.data()+1, X.data()+X.size() );
  int i = boost::lower_bound( rng, x) - X.data() - 1;

  // don't extrapolate at all
  if( i == -1 || i == X.size() - 1 )
    return 0;
     
  Real b  = Y(i);
  Real m  = Y(i+1)-Y(i);
       m /= X(i+1)-X(i);
  
  return m*(x - X(i)) + b;
}

template<typename Real>
Real
LinearInterpolator<Real>::derivative_imp( Real x ) const
{
  const VectorType &X = *(this->xv);
  const VectorType &Y = *(this->yv);

  // find the index that is just to the left of the x
  //int i = Utils::index__last_lt( x, X, X.size(), 0);
  auto rng = std::make_pair( X.data()+1, X.data()+X.size() );
  int i = boost::lower_bound( rng, x) - X.data() - 1;

  // don't extrapolate at all
  if( i == -1 || i == X.size() - 1 )
    return 0;

  Real m  = Y(i+1)-Y(i);
       m /= X(i+1)-X(i);
  
  return m;

}








}

#endif
