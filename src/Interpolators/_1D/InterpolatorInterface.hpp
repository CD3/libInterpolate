
#ifndef Interpolators__1D_InterpolatorInterface_hpp
#define Interpolators__1D_InterpolatorInterface_hpp

/** @file InterpolatorInterface.hpp
  * @author C.D. Clark III
  * @date 12/24/16
  */

#include <vector>
#include <Eigen/Dense>

namespace _1D {

/** @class 
  * @brief Interface definition for 1D interpolation classes.
  * @author C.D. Clark III
  */
template<class Real>
class InterpolatorInterface
{
  public:
    virtual Real operator()( Real x ) const = 0;
    virtual void setData( size_t n, Real *x, Real *y, bool deep_copy = true ) = 0;
};

}

#endif // include protector
