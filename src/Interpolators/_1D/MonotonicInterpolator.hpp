#ifndef Interpolators__1D_MonotonicInterpolator_hpp
#define Interpolators__1D_MonotonicInterpolator_hpp

/** @file MonotonicInterpolator.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 06/06/17
  */

#include "InterpolatorBase.hpp"
#include <boost/range/algorithm/lower_bound.hpp>

namespace _1D {

/** @class MonotonicInterpolator

  * @brief A simple monotonic interpolator based on Steffen, "A simple method for monotonic interpolation in one dimension", 1990.
  * @author C.D. Clark III
  */
template<class Real>
class MonotonicInterpolator : public InterpolatorBase<MonotonicInterpolator<Real>>
{

  public:
    using BASE = InterpolatorBase<MonotonicInterpolator<Real>>;
    using VectorType = typename BASE::VectorType;
    using MapType = typename BASE::MapType;

    Real operator()( Real x ) const;

    MonotonicInterpolator() = default;
    MonotonicInterpolator( const MonotonicInterpolator& interp ) = default;

    template<typename I>
    MonotonicInterpolator( I n, Real *x, Real *y, bool deep_copy = true )
    { this->setData(n,x,y,deep_copy); }
    template<typename X, typename Y>
    MonotonicInterpolator( X &x, Y &y, bool deep_copy = true )
    { this->setData(x,y,deep_copy); }

  protected:

    void setupInterpolator();
    friend BASE;

  protected:
    VectorType a,b,yplow,yphigh;
};

template<class Real>
void
MonotonicInterpolator<Real>::setupInterpolator()
{
  const VectorType &X = *(this->xv);
  const VectorType &Y = *(this->yv);

  a      = VectorType(X.size()-1);
  b      = VectorType(X.size()-1);
  yplow  = VectorType(X.size()-1);
  yphigh = VectorType(X.size()-1);


  for(int i = 0; i < X.size()-1; i++)
  {
    // i is the *interval* index
    // the interval spans from X(i) to X(i+1)
    Real xlow  = X(i);
    Real xhigh = X(i+1);
    Real ylow  = Y(i);
    Real yhigh = Y(i+1);
    Real h = xhigh - xlow;

    Real slope = (yhigh - ylow)/h;

    // Determine the first derivative values used at the endpoints of the interval
    // These values may be limited to preserve monotonicity
    if (i==0)
    {
      // first interval does not have an interval to its "left", so just
      // use the slope in the interval.
      yplow[i] = slope;
    }
    else
    {
      // Determine the lower slope value
      Real hlow = xlow - X(i-1);
      Real slope_low = 0.0;
      if (hlow > 0.0)
        slope_low = (ylow - Y(i-1))/hlow;
      if (slope_low* slope <= 0.0)
      {
        // Set derivative as zero
        yplow[i] = 0.0;
      }
      else
      {
        yplow[i] = ((slope_low*h) + (slope*hlow))/(hlow + h);
        if (yplow[i] >= 0.0)
        {
          yplow[i] = (std::min)(yplow[i], 2.0*(std::min)(slope_low, slope));
        }
        else
        {
          yplow[i] = (std::max)(yplow[i], 2.0*(std::max)(slope_low, slope));
        }
      }
    }

  
	if(i == (X.size()-2))
	{
    // last interval does not have an interval to its "right", so just
    // use the slope in the interval.
		yphigh[i] = slope;
	}
	else
	{
		// Determine the upper slope value
		Real hhigh = Y(i+2) - xhigh;
		Real slope_high = 0.0;
		if (hhigh > 0.0)
			slope_high = (Y(i+2) - yhigh)/hhigh;
		if (slope*slope_high <= 0.0)
		{
			// Set derivative as zero
			yphigh[i] = 0.0;
		}
		else
		{
			yphigh[i] = ((slope*hhigh) + (slope_high*h))/(h + hhigh);
			if (yphigh[i] >= 0.0)
			{
				yphigh[i] = (std::min)(yphigh[i], 2.0*(std::min)(slope, slope_high));
			}
			else
			{
				yphigh[i] = (std::max)(yphigh[i], 2.0*(std::max)(slope, slope_high));
			}
		}
	}

  a[i] = (yplow[i] + yphigh[i] - (2.0*slope))/(h*h);
	b[i] = ((3.0*slope) + (-2.0*yplow[i]) - yphigh[i])/h;



  }



}

template<class Real>
Real
MonotonicInterpolator<Real>::operator()( Real x ) const
{
  BASE::checkData();

  const VectorType &X = *(this->xv);
  const VectorType &Y = *(this->yv);

  // find the index that is just to the left of the x
  // this will correspond to the "interval index"
  //int i = Utils::index__first_ge( x, X, X.size(), 1);
  auto rng = std::make_pair( X.data()+1, X.data()+X.size() );
  int i = boost::lower_bound( rng, x ) - X.data();

  // don't extrapolate at all
  if( i == 0 || i == X.size())
    return 0;
  i--; // we need the interval index, not the right point index.

	// Deal with the degenerate case of xval = xlow = xhigh
	if (X(i+1) <= X(i))
		return 0.5*(Y(i) + Y(i+1));

	Real xm = x - X(i);
	Real xm2 = xm*xm;
	Real f = (a[i]*xm*xm2) + (b[i]*xm2) + (yplow[i]*xm) + Y(i);

  
  return f;
}

}

#endif // include protector
