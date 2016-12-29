#ifndef Interpolators__2D_InterpolatorInterface_hpp
#define Interpolators__2D_InterpolatorInterface_hpp

/** @file InterpolatorInterface.hpp
  * @author C.D. Clark III
  * @date 12/24/16
  */

#include <vector>
#include <eigen3/Eigen/Dense>

namespace _2D {

/** @class 
  * @brief Interface definition for 2D interpolation classes.
  * @author C.D. Clark III
  */
template<class Real>
class InterpolatorInterface
{
  public:
    virtual Real operator()( Real x, Real y ) const = 0;
    virtual void setData( size_t _n, Real *x, Real *y, Real *z, bool deep_copy = true ) = 0;
};

}

#endif // include protector
