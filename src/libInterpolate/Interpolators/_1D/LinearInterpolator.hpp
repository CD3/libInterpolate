#ifndef Interpolators__1D_LinearInterpolator_hpp
#define Interpolators__1D_LinearInterpolator_hpp

#include "InterpolatorBase.hpp"

namespace _1D {


/** @class 
  * @brief A simple linear interpolator.
  * @author C.D. Clark III
  *
  * This class does *not* do extrapolation.
  */
template<class Real>
class LinearInterpolator : public InterpolatorBase<LinearInterpolator<Real>>
{

  public:
    using BASE = InterpolatorBase<LinearInterpolator<Real>>;
    using VectorType = typename BASE::VectorType;
    using MapType = typename BASE::MapType;

    template<typename I>
    LinearInterpolator( I n, Real *x, Real *y ) { this->setData(n,x,y); }

    template<typename X, typename Y>
    LinearInterpolator( X &x, Y &y ) { this->setData(x,y); }

    LinearInterpolator( ):BASE(){}

    LinearInterpolator(const LinearInterpolator& rhs) = default;
    LinearInterpolator& operator=(const LinearInterpolator& rhs) = default;


    Real operator()( Real x ) const;





};

template<class Real>
Real
LinearInterpolator<Real>::operator()( Real x ) const
{
  BASE::checkData();

  // don't extrapolate at all
  if( x < this->xData[0] || x > this->xData[this->xData.size()-1] )
    return 0;

  const MapType &X = *(this->xView);
  const MapType &Y = *(this->yView);

  int i = this->get_index_to_left_of(x);
     
  Real b  = Y(i);
  Real m  = Y(i+1)-Y(i);
       m /= X(i+1)-X(i);
  
  return m*(x - X(i)) + b;
}





}

#endif
