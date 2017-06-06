#ifndef Interpolators__1D_MonotonicInterpolator_hpp
#define Interpolators__1D_MonotonicInterpolator_hpp

/** @file MonotonicInterpolator.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 06/06/17
  */

#include "InterpolatorBase.hpp"

namespace _1D {

/** @class MonotonicInterpolator

  * @brief A simple monotonic interpolator based on Steffen, "A simple method for monotonic interpolation in one dimension", 1990.
  * @author C.D. Clark III
  */
template<class Real>
class MonotonicInterpolator : public InterpolatorBase<Real>
{
  public:
    typedef typename InterpolatorBase<Real>::VectorType VectorType;
    typedef typename InterpolatorBase<Real>::MapType MapType;

    using InterpolatorBase<Real>::operator();
    virtual Real operator()( Real x ) const;

    virtual void setData( size_t _n, Real *x, Real *y, bool deep_copy = true );
    virtual void setData( std::vector<Real> &x, std::vector<Real> &y, bool deep_copy = true );
    virtual void setData( VectorType  &x, VectorType &y, bool deep_copy = true );

  protected:

    void calcCoefficients();
};

template<class Real>
void
MonotonicInterpolator<Real>::setData( size_t n, Real *x, Real *y, bool deep_copy )
{
  InterpolatorBase<Real>::setData( n, x, y, deep_copy );
  calcCoefficients();
}

template<class Real>
void
MonotonicInterpolator<Real>::setData( std::vector<Real> &x, std::vector<Real> &y, bool deep_copy )
{
  InterpolatorBase<Real>::setData( x, y, deep_copy );
  calcCoefficients();
}

template<class Real>
void
MonotonicInterpolator<Real>::setData( VectorType  &x, VectorType &y, bool deep_copy )
{
  InterpolatorBase<Real>::setData( x, y, deep_copy );
  calcCoefficients();
}

template<class Real>
void
MonotonicInterpolator<Real>::calcCoefficients()
{
}

template<class Real>
Real
MonotonicInterpolator<Real>::operator()( Real x ) const
{
  InterpolatorBase<Real>::checkData();

  const VectorType &X = *(this->xv);
  const VectorType &Y = *(this->yv);

  // find the index that is just to the right of the x
  int i = Utils::index_first_ge( x, X, 1);

  // don't extrapolate at all
  if( i == 0 || i == X.size())
    return 0;

  Real xlow  = X(i-1);
  Real xhigh = X(i);
  Real ylow  = Y(i-1);
  Real yhigh = Y(i);
  Real h = xhigh - xlow;

	// Deal with the degenerate case of xval = xlow = xhigh
	if (h <= 0.0)
		return (Real(0.5*(ylow + yhigh)));

  Real slope = (yhigh - ylow)/h;

	// Determine the first derivative values used at the endpoints of the interval
	// These values may be limited to preserve monotonicity
	Real yplow, yphigh;
	if (i==1)
	{
		// Extrapolate the slope
		yplow = slope;
	}
	else
	{
		// Determine the lower slope value
		Real hlow = xlow - X(i-2);
		Real slope_low = 0.0;
		if (hlow > 0.0)
			slope_low = (ylow - Y(i-2))/hlow;
		if (slope_low* slope <= 0.0)
		{
			// Set derivative as zero
			yplow = 0.0;
		}
		else
		{
			yplow = ((slope_low*h) + (slope*hlow))/(hlow + h);
			if (yplow >= 0.0)
			{
				yplow = (std::min)(yplow, 2.0*(std::min)(slope_low, slope));
			}
			else
			{
				yplow = (std::max)(yplow, 2.0*(std::max)(slope_low, slope));
			}
		}
	}

	if (i == (X.size()-1))
	{
		// Extrapolate the slope
		yphigh = slope;
	}
	else
	{
		// Determine the upper slope value
		Real hhigh = Y(i+1) - xhigh;
		Real slope_high = 0.0;
		if (hhigh > 0.0)
			slope_high = (Y(i+1) - yhigh)/hhigh;
		if (slope*slope_high <= 0.0)
		{
			// Set derivative as zero
			yphigh = 0.0;
		}
		else
		{
			yphigh = ((slope*hhigh) + (slope_high*h))/(h + hhigh);
			if (yphigh >= 0.0)
			{
				yphigh = (std::min)(yphigh, 2.0*(std::min)(slope, slope_high));
			}
			else
			{
				yphigh = (std::max)(yphigh, 2.0*(std::max)(slope, slope_high));
			}
		}
	}

	//Solve for the interpolated value as a cubic polynomic
	Real a = (yplow + yphigh - (2.0*slope))/(h*h);
	Real b = ((3.0*slope) + (-2.0*yplow) - yphigh)/h;
	Real xm = x - xlow;
	Real xm2 = xm*xm;
	Real f = (a*xm*xm2) + (b*xm2) + (yplow*xm) + ylow;

  
  return f;
}

}

#endif // include protector
