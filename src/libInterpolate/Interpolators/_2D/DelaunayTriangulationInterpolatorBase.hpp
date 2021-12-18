#pragma once

/** @file DelaunayTriangulationInterpolatorBase.hpp
  * @brief A base class for interpolators based on triangles.
  * @author C.D. Clark III
  * @date 2021-12-15
  */

#include <iostream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/ring.hpp>

#include "../../Utils/Meshing/delaunator-cpp.hpp"
#include "InterpolatorBase.hpp"

namespace _2D
{
template<class Derived, typename Real = typename RealTypeOf<Derived>::type>
class DelaunayTriangulationInterpolatorBase : public InterpolatorBase<DelaunayTriangulationInterpolatorBase<Derived, Real>>
{
 public:
  using BASE       = InterpolatorBase<DelaunayTriangulationInterpolatorBase<Derived, Real>>;
  using point_t    = boost::geometry::model::d2::point_xy<Real>;
  using triangle_t = boost::geometry::model::ring<point_t>;

  std::vector<triangle_t> getTriangles() const { return triangles; }

 private:
  friend BASE;     // this is necessary to allow base class to call setupInterpolator()
  friend Derived;  // this is necessary to allow derived class to call constructors

  std::vector<triangle_t> triangles;

  void setupInterpolator()
  {
    std::vector<Real> coords;
    coords.reserve(2 * this->xData.size());
    for(size_t i = 0; i < this->xData.size(); ++i) {
      coords.push_back(this->xData[i]);
      coords.push_back(this->yData[i]);
    }

    delaunator::Delaunator triangulation(coords);

    for(size_t i = 0; i < triangulation.triangles.size(); i += 3) {
      point_t p1(triangulation.coords[2 * triangulation.triangles[i]], triangulation.coords[2 * triangulation.triangles[i] + 1]);
      point_t p2(triangulation.coords[2 * triangulation.triangles[i + 1]], triangulation.coords[2 * triangulation.triangles[i + 1] + 1]);
      point_t p3(triangulation.coords[2 * triangulation.triangles[i + 2]], triangulation.coords[2 * triangulation.triangles[i + 2] + 1]);

      triangle_t t{p1, p2, p3, p1};

      triangles.push_back(t);
    }

    this->template callSetupInterpolator<Derived>();
  }
};
}  // namespace _2D
