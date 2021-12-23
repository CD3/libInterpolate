#pragma once

#include "DelaunayTriangulationInterpolatorBase.hpp"

namespace _2D {

template<class Real>
class LinearDelaunayTriangleInterpolator : public DelaunayTriangulationInterpolatorBase<LinearDelaunayTriangleInterpolator<Real>>
{
  public:

    using BASE = DelaunayTriangulationInterpolatorBase<LinearDelaunayTriangleInterpolator<Real>>;
    using VectorType = typename BASE::VectorType;
    using MapType = typename BASE::MapType; 

    using point_t = typename BASE::point_t;
    using triangle_t = typename BASE::triangle_t;

    // types used to view data as 2D coordinates
    using MatrixType    = typename BASE::MatrixType;
    using _2DVectorView = typename BASE::_2DVectorView;
    using _2DMatrixView = typename BASE::_2DMatrixView;

    // types used for 2x2 matrix algebra
    using Matrix33 = Eigen::Matrix<Real,3,3 >;
    using ColVector2 = Eigen::Matrix<Real,3,1 >;
    using RowVector2 = Eigen::Matrix<Real,1,3 >;

    Real operator()( Real x, Real y) const
    {

      point_t p{x,y};
      for( auto& t: this->triangles )
      {
        if( boost::geometry::covered_by(p,t) )
        {
          // std::array<std::array<Real,3>,4> points;
          // int i = 0;
          // boost::geometry::for_each_point(t,[&i,&points](const point_t& pp){
          //     points[i++] = pp;
          //     });
          // std::cout << boost::geometry::wkt( points[0] ) << std::endl;
          // std::cout << boost::geometry::wkt( points[1] ) << std::endl;
          // boost::geometry::subtract_point(points[0],points[1]);
          // boost::geometry::subtract_point(points[1],points[2]);
          // std::cout << boost::geometry::wkt( points[0] ) << std::endl;
          // std::cout << boost::geometry::wkt( points[1] ) << std::endl;
          // std::cout << boost::geometry::wkt( boost::geometry::cross_product(points[0],points[1]) ) << std::endl;
          // return 1;
        }
      }

      // do not extrapolate
      return 0;
    }

  protected:
    using BASE::triangles;

  private:
    friend BASE;
    void setupInterpolator()
    {

    }





};

}
