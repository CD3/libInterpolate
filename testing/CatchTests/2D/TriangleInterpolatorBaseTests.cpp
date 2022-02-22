#include "catch.hpp"

#include <libInterpolate/Interpolators/_2D/DelaunayTriangulationInterpolatorBase.hpp>

#include <numeric>

namespace _2D
{
class TestTriangleInterp : public DelaunayTriangulationInterpolatorBase<TestTriangleInterp>
{
 protected:
  double c = 0;

 public:
  virtual double operator()(double x, double y) const
  {
    return x + 2 * y + this->c;
  }

  VectorType getX() { return *(this->xView); }
  VectorType getY() { return *(this->yView); }
  VectorType getZ() { return *(this->zView); }

  void setupInterpolator()
  {
    this->c = 10;
  }
};

}  // namespace _2D

TEST_CASE("2D DelaunayTriangulationInterpolatorBase Setup Tests", "[plumbing]")
{
  _2D::TestTriangleInterp interp;

  // make sure interpolator works the way we expect
  REQUIRE(interp(1, 1) == Approx(3));
  REQUIRE(interp(10, 20) == Approx(50));

  SECTION("Grid-formatted Data")
  {
    std::vector<int> x(3), y(4), z(12);
    x[0] = 0;
    x[1] = 1;
    x[2] = 2;
    y[0] = 2;
    y[1] = 3;
    y[2] = 4;
    y[3] = 5;

    SECTION("Iterators")
    {
      interp.setData(x.begin(), x.end(), y.begin(), y.end(), z.begin(), z.end());
      auto X = interp.getX();
      auto Y = interp.getY();
      auto Z = interp.getZ();

      CHECK(X.size() == 12);
      CHECK(X[0] == 0);
      CHECK(X[1] == 0);
      CHECK(X[2] == 0);
      CHECK(X[3] == 0);

      CHECK(X[4] == 1);
      CHECK(X[5] == 1);
      CHECK(X[6] == 1);
      CHECK(X[7] == 1);

      CHECK(X[8] == 2);
      CHECK(X[9] == 2);
      CHECK(X[10] == 2);
      CHECK(X[11] == 2);

      CHECK(Y[0] == 2);
      CHECK(Y[1] == 3);
      CHECK(Y[2] == 4);
      CHECK(Y[3] == 5);
      CHECK(Y[4] == 2);
      CHECK(Y[5] == 3);
      CHECK(Y[6] == 4);
      CHECK(Y[7] == 5);
      CHECK(Y[8] == 2);
      CHECK(Y[9] == 3);
      CHECK(Y[10] == 4);
      CHECK(Y[11] == 5);
    }
    SECTION("std::vector")
    {
      interp.setData(x, y, z);
      auto X = interp.getX();
      auto Y = interp.getY();
      auto Z = interp.getZ();

      CHECK(X.size() == 12);
      CHECK(X[0] == 0);
      CHECK(X[1] == 0);
      CHECK(X[2] == 0);
      CHECK(X[3] == 0);

      CHECK(X[4] == 1);
      CHECK(X[5] == 1);
      CHECK(X[6] == 1);
      CHECK(X[7] == 1);

      CHECK(X[8] == 2);
      CHECK(X[9] == 2);
      CHECK(X[10] == 2);
      CHECK(X[11] == 2);

      CHECK(Y[0] == 2);
      CHECK(Y[1] == 3);
      CHECK(Y[2] == 4);
      CHECK(Y[3] == 5);
      CHECK(Y[4] == 2);
      CHECK(Y[5] == 3);
      CHECK(Y[6] == 4);
      CHECK(Y[7] == 5);
      CHECK(Y[8] == 2);
      CHECK(Y[9] == 3);
      CHECK(Y[10] == 4);
      CHECK(Y[11] == 5);

      using point_t = _2D::DelaunayTriangulationInterpolatorBase<_2D::TestTriangleInterp>::point_t;
      using triangle_t = _2D::DelaunayTriangulationInterpolatorBase<_2D::TestTriangleInterp>::triangle_t;
      auto                                                                    triangles = interp.getTriangles();
      std::vector<int> covered_counts(triangles.size());
      point_t p{2., 3.};
      std::transform(triangles.begin(),triangles.end(),covered_counts.begin(), [&p](const triangle_t& t){ return boost::geometry::covered_by(p,t); } );
      CHECK( std::accumulate(covered_counts.begin(), covered_counts.end(), 0, std::plus<int>()) == 3 );

      p = point_t{2.-0.1, 3.};
      std::transform(triangles.begin(),triangles.end(),covered_counts.begin(), [&p](const triangle_t& t){ return boost::geometry::covered_by(p,t); } );
      CHECK( std::accumulate(covered_counts.begin(), covered_counts.end(), 0, std::plus<int>()) == 2 );

      p = point_t{2.-0.1, 3.-0.1};
      std::transform(triangles.begin(),triangles.end(),covered_counts.begin(), [&p](const triangle_t& t){ return boost::geometry::covered_by(p,t); } );
      CHECK( std::accumulate(covered_counts.begin(), covered_counts.end(), 0, std::plus<int>()) == 1 );
    }
  }
}
