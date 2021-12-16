#pragma once

/** @file TriangleInterpolatorBase.hpp
  * @brief A base class for interpolators based on triangles.
  * @author C.D. Clark III
  * @date 2021-12-15
  */

#include <iostream>

#include "InterpolatorBase.hpp"

namespace _2D
{
template<class Derived, typename Real = typename RealTypeOf<Derived>::type>
class TriangleInterpolatorBase : public InterpolatorBase<TriangleInterpolatorBase<Derived, Real>>
{
 public:
  using BASE = InterpolatorBase<TriangleInterpolatorBase<Derived, Real>>;

 private:
  friend BASE;    // this is necessary to allow base class to call setupInterpolator()
  friend Derived; // this is necessary to allow derived class to call constructors

  void setupInterpolator()
  {
    this->template callSetupInterpolator<Derived>();
  }
};
}  // namespace _2D
