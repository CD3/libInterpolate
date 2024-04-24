#pragma once

/** @file DelaunayTriangulationInterpolatorBase.hpp
 * @brief A base class for interpolators based on triangles.
 * @author C.D. Clark III
 * @date 2021-12-15
 */

#include <array>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/adapted/std_array.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <iostream>
// #include <boost/geometry/geometries/point_xy.hpp>
BOOST_GEOMETRY_REGISTER_STD_ARRAY_CS(cs::cartesian)

#include "../../Utils/Meshing/delaunator-cpp.hpp"
#include "InterpolatorBase.hpp"

namespace _2D {
template <class Derived, typename Real = typename RealTypeOf<Derived>::type>
class DelaunayTriangulationInterpolatorBase
    : public InterpolatorBase<
          DelaunayTriangulationInterpolatorBase<Derived, Real>> {
   public:
    using BASE =
        InterpolatorBase<DelaunayTriangulationInterpolatorBase<Derived, Real>>;
    using point_t =
        std::array<Real, 2>;  // boost::geometry::model::d2::point_xy<Real>;
    using triangle_t = boost::geometry::model::ring<point_t>;
    using box_t = boost::geometry::model::box<point_t>;
    using rtree_value_t = std::pair<box_t, size_t>;

    std::vector<triangle_t> getTriangles() const { return m_xy_triangles; }

   private:
    friend BASE;     // this is necessary to allow base class to call
                     // setupInterpolator()
    friend Derived;  // this is necessary to allow derived class to call
                     // constructors

    std::vector<std::array<size_t, 3>> m_triangle_datapoints;
    boost::geometry::index::rtree<rtree_value_t,
                                  boost::geometry::index::quadratic<16>>
        m_triangles_index;

   protected:
    std::vector<triangle_t>
        m_xy_triangles;  // WORKAROUND: MSVC does not let derived class access
                         // this if it is private

    void setupInterpolator() {
        clearTriangleBuffers();  // Clear triangle buffers before initialization

        std::vector<Real> coords;
        coords.reserve(2 * this->xData.size());
        for (size_t i = 0; i < this->xData.size(); ++i) {
            coords.push_back(this->xData[i]);
            coords.push_back(this->yData[i]);
        }

        delaunator::Delaunator<Real> triangulation(coords);

        for (size_t i = 0; i < triangulation.triangles.size(); i += 3) {
            point_t p1 = {
                triangulation.coords[2 * triangulation.triangles[i]],
                triangulation.coords[2 * triangulation.triangles[i] + 1]};
            point_t p2 = {
                triangulation.coords[2 * triangulation.triangles[i + 1]],
                triangulation.coords[2 * triangulation.triangles[i + 1] + 1]};
            point_t p3 = {
                triangulation.coords[2 * triangulation.triangles[i + 2]],
                triangulation.coords[2 * triangulation.triangles[i + 2] + 1]};

            triangle_t t{p1, p2, p3, p1};

            m_xy_triangles.push_back(t);
            m_triangle_datapoints.push_back({triangulation.triangles[i],
                                             triangulation.triangles[i + 1],
                                             triangulation.triangles[i + 2]});

            box_t b = boost::geometry::return_envelope<box_t>(t);
            m_triangles_index.insert(
                std::make_pair(b, m_xy_triangles.size() - 1));
        }

        this->template callSetupInterpolator<Derived>();
    }

    void clearTriangleBuffers() {
        m_xy_triangles.clear();
        m_triangle_datapoints.clear();
        m_triangles_index.clear();
    }
};
}  // namespace _2D
