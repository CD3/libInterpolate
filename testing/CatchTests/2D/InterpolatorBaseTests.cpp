#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <libInterpolate/Interpolators/_2D/InterpolatorBase.hpp>
#include <limits>
using namespace Catch;
namespace _2D {

class TestInterp : public InterpolatorBase<TestInterp> {
   public:
    virtual double operator()(double x, double y) const {
        return x + 2 * y + 10;
    }

    VectorType getX() { return *(this->xView); }
    VectorType getY() { return *(this->yView); }
    VectorType getZ() { return *(this->zView); }
};

}  // namespace _2D

TEST_CASE("2D InterpolatorBase Setup Tests", "[plumbing]") {
    _2D::TestInterp interp;

    // make sure interpolator works the way we expect
    REQUIRE(interp(1, 1) == Approx(13));
    REQUIRE(interp(10, 20) == Approx(60));

    size_t N = 10;

    SECTION("Eigen Vector Initialization") {
        _2D::InterpolatorBase<double>::VectorType xx(N), yy(N), zz(N);

        for (size_t i = 0; i < N; i++) {
            xx(i) = 0.1 * i;
            yy(i) = 0.2 * i;
            zz(i) = xx(i) + yy(i);
        }

        SECTION("Deep Copy") {
            interp.setData(xx, yy, zz);

            // clear the original data to make sure deep copy worked
            for (size_t i = 0; i < N; i++) xx(i) = yy(i) = 0;

            auto x = interp.getX();
            auto y = interp.getY();
            auto z = interp.getZ();

            for (size_t i = 0; i < N; i++) {
                REQUIRE(x(i) == Approx(0.1 * i));
                REQUIRE(y(i) == Approx(0.2 * i));
                REQUIRE(z(i) == Approx(x(i) + y(i)));
            }
        }

        // SECTION("Shallow Copy")
        //{
        // interp.setData( xx, yy, zz, false );

        //// clear the original data to check that shallow copy worked
        // for(int i = 0; i < N; i++)
        // xx(i) = yy(i) = zz(i) = 0;

        // auto x = interp.getX();
        // auto y = interp.getY();
        // auto z = interp.getZ();

        // for( int i = 0; i < N; i++)
        //{
        // REQUIRE( x(i) == Approx( 0.0 ) );
        // REQUIRE( y(i) == Approx( 0.0 ) );
        // REQUIRE( z(i) == Approx( 0.0 ) );
        //}
        //}
    }

    SECTION("Std Vector Initialization") {
        std::vector<double> xx(N), yy(N), zz(N);

        for (size_t i = 0; i < N; i++) {
            xx[i] = 0.1 * i;
            yy[i] = 0.2 * i;
            zz[i] = xx[i] + yy[i];
        }

        SECTION("Deep Copy") {
            interp.setData(xx, yy, zz);

            // clear the original data to make sure deep copy worked
            for (size_t i = 0; i < N; i++) xx[i] = yy[i] = zz[i] = 0;

            auto x = interp.getX();
            auto y = interp.getY();
            auto z = interp.getZ();

            for (size_t i = 0; i < N; i++) {
                REQUIRE(x(i) == Approx(0.1 * i));
                REQUIRE(y(i) == Approx(0.2 * i));
                REQUIRE(z(i) == Approx(x(i) + y(i)));
            }
        }

        // SECTION("Shallow Copy")
        //{
        // interp.setData( xx, yy, zz, false );

        //// clear the original data to check that shallow copy worked
        // for(int i = 0; i < N; i++)
        // xx[i] = yy[i] = zz[i] = 0;

        // auto x = interp.getX();
        // auto y = interp.getY();
        // auto z = interp.getZ();

        // for( int i = 0; i < N; i++)
        //{
        // REQUIRE( x(i) == Approx( 0.0 ) );
        // REQUIRE( y(i) == Approx( 0.0 ) );
        // REQUIRE( z(i) == Approx( 0.0 ) );
        //}
        //}
    }

    SECTION("Raw Pointer Initialization") {
        double *xx, *yy, *zz;
        xx = new double[N];
        yy = new double[N];
        zz = new double[N];

        for (size_t i = 0; i < N; i++) {
            xx[i] = 0.1 * i;
            yy[i] = 0.2 * i;
            zz[i] = xx[i] + yy[i];
        }

        SECTION("Deep Copy") {
            interp.setData(N, xx, yy, zz);

            // clear the original data to make sure deep copy worked
            for (size_t i = 0; i < N; i++) xx[i] = yy[i] = zz[i] = 0;

            auto x = interp.getX();
            auto y = interp.getY();
            auto z = interp.getZ();

            for (size_t i = 0; i < N; i++) {
                REQUIRE(x(i) == Approx(0.1 * i));
                REQUIRE(y(i) == Approx(0.2 * i));
                REQUIRE(z(i) == Approx(x(i) + y(i)));
            }
        }

        // SECTION("Shallow Copy")
        //{
        // interp.setData( N, xx, yy, zz, false );

        //// clear the original data to check that shallow copy worked
        // for(int i = 0; i < N; i++)
        // for(int i = 0; i < N; i++)
        // xx[i] = yy[i] = zz[i] = 0;

        // auto x = interp.getX();
        // auto y = interp.getY();
        // auto z = interp.getZ();

        // for( int i = 0; i < N; i++)
        //{
        // REQUIRE( x(i) == Approx( 0.0 ) );
        // REQUIRE( y(i) == Approx( 0.0 ) );
        // REQUIRE( z(i) == Approx( 0.0 ) );
        //}
        //}

        delete[] xx;
        delete[] yy;
        delete[] zz;
    }

    SECTION("Eigen Vector Initialization") {
        _2D::InterpolatorBase<double>::VectorType xx(N), yy(N), zz(N);

        for (size_t i = 0; i < N; i++) {
            xx(i) = 0.1 * i;
            yy(i) = 0.2 * i;
            zz(i) = xx(i) + yy(i);
        }

        interp.setData(xx, yy, zz);

        auto x = interp.getXData();
        auto y = interp.getYData();
        auto z = interp.getZData();

        REQUIRE(x.size() > 0);
        CHECK(x[0] == Approx(0));
        REQUIRE(x.size() == N);
        CHECK(x[N - 1] == Approx(xx(N - 1)));

        REQUIRE(y.size() > 0);
        CHECK(y[0] == Approx(0));
        REQUIRE(y.size() == N);
        CHECK(y[N - 1] == Approx(yy(N - 1)));

        REQUIRE(z.size() > 0);
        CHECK(z[0] == Approx(0));
        REQUIRE(z.size() == N);
        CHECK(z[N - 1] == Approx(xx(N - 1) + yy(N - 1)));
    }

    SECTION("Grid-formatted Data") {
        std::vector<int> x(3), y(4), z(12);
        x[0] = 0;
        x[1] = 1;
        x[2] = 2;
        y[0] = 2;
        y[1] = 3;
        y[2] = 4;
        y[3] = 5;
        z[0] = 10;
        z[1] = 11;
        z[2] = 12;
        z[3] = 13;
        z[4] = 14;
        z[5] = 15;
        z[6] = 16;
        z[7] = 17;
        z[8] = 18;
        z[9] = 19;
        z[10] = 20;
        z[11] = 21;

        SECTION("Iterators") {
            interp.setData(x.begin(), x.end(), y.begin(), y.end(), z.begin(),
                           z.end());
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

            CHECK(Z[0] == 10);
            CHECK(Z[1] == 11);
            CHECK(Z[2] == 12);
            CHECK(Z[3] == 13);
            CHECK(Z[4] == 14);
            CHECK(Z[5] == 15);
            CHECK(Z[6] == 16);
            CHECK(Z[7] == 17);
            CHECK(Z[8] == 18);
            CHECK(Z[9] == 19);
            CHECK(Z[10] == 20);
            CHECK(Z[11] == 21);
        }
        SECTION("std::vector") {
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

            CHECK(Z[0] == 10);
            CHECK(Z[1] == 11);
            CHECK(Z[2] == 12);
            CHECK(Z[3] == 13);
            CHECK(Z[4] == 14);
            CHECK(Z[5] == 15);
            CHECK(Z[6] == 16);
            CHECK(Z[7] == 17);
            CHECK(Z[8] == 18);
            CHECK(Z[9] == 19);
            CHECK(Z[10] == 20);
            CHECK(Z[11] == 21);
        }
    }
}

TEST_CASE("2D InterpolatorBase closest index helpers", "[plumbing]") {
    _2D::TestInterp interp;

    // unevenly spaced grid so that "closer to the lower knot", "closer to the
    // upper knot" and "exactly halfway" are different cases.
    // x = {0, 1, 3}, y = {2, 3, 5, 10}
    std::vector<double> x, y, z;
    for (double xx : {0., 1., 3.}) {
        for (double yy : {2., 3., 5., 10.}) {
            x.push_back(xx);
            y.push_back(yy);
            z.push_back(xx + yy);
        }
    }
    interp.setData(x, y, z);

    const double inf = std::numeric_limits<double>::infinity();
    auto below = [&](double v) { return std::nextafter(v, -inf); };
    auto above = [&](double v) { return std::nextafter(v, inf); };

    SECTION("x") {
        // knots
        CHECK(interp.get_x_index_closest_to(0) == 0);
        CHECK(interp.get_x_index_closest_to(1) == 1);
        CHECK(interp.get_x_index_closest_to(3) == 2);
        // between knots
        CHECK(interp.get_x_index_closest_to(0.25) == 0);
        CHECK(interp.get_x_index_closest_to(0.75) == 1);
        CHECK(interp.get_x_index_closest_to(1.5) == 1);
        CHECK(interp.get_x_index_closest_to(2.5) == 2);
        // ties go to the upper knot
        CHECK(interp.get_x_index_closest_to(0.5) == 1);
        CHECK(interp.get_x_index_closest_to(2) == 2);
        // just either side of a knot
        CHECK(interp.get_x_index_closest_to(above(0)) == 0);
        CHECK(interp.get_x_index_closest_to(below(1)) == 1);
        CHECK(interp.get_x_index_closest_to(above(1)) == 1);
        CHECK(interp.get_x_index_closest_to(below(3)) == 2);
        // outside the data clamps to the first or last index
        CHECK(interp.get_x_index_closest_to(below(0)) == 0);
        CHECK(interp.get_x_index_closest_to(-100) == 0);
        CHECK(interp.get_x_index_closest_to(-inf) == 0);
        CHECK(interp.get_x_index_closest_to(above(3)) == 2);
        CHECK(interp.get_x_index_closest_to(100) == 2);
        CHECK(interp.get_x_index_closest_to(inf) == 2);
        // NaN must give a valid index. Which one is not part of the contract.
        int i = interp.get_x_index_closest_to(std::nan(""));
        CHECK(i >= 0);
        CHECK(i <= 2);
    }

    SECTION("y") {
        // knots
        CHECK(interp.get_y_index_closest_to(2) == 0);
        CHECK(interp.get_y_index_closest_to(3) == 1);
        CHECK(interp.get_y_index_closest_to(5) == 2);
        CHECK(interp.get_y_index_closest_to(10) == 3);
        // between knots
        CHECK(interp.get_y_index_closest_to(2.25) == 0);
        CHECK(interp.get_y_index_closest_to(2.75) == 1);
        CHECK(interp.get_y_index_closest_to(3.5) == 1);
        CHECK(interp.get_y_index_closest_to(4.5) == 2);
        CHECK(interp.get_y_index_closest_to(7) == 2);
        CHECK(interp.get_y_index_closest_to(8) == 3);
        // ties go to the upper knot
        CHECK(interp.get_y_index_closest_to(2.5) == 1);
        CHECK(interp.get_y_index_closest_to(4) == 2);
        CHECK(interp.get_y_index_closest_to(7.5) == 3);
        // just either side of a knot
        CHECK(interp.get_y_index_closest_to(above(2)) == 0);
        CHECK(interp.get_y_index_closest_to(below(3)) == 1);
        CHECK(interp.get_y_index_closest_to(above(3)) == 1);
        CHECK(interp.get_y_index_closest_to(below(10)) == 3);
        // outside the data clamps to the first or last index
        CHECK(interp.get_y_index_closest_to(below(2)) == 0);
        CHECK(interp.get_y_index_closest_to(-100) == 0);
        CHECK(interp.get_y_index_closest_to(-inf) == 0);
        CHECK(interp.get_y_index_closest_to(above(10)) == 3);
        CHECK(interp.get_y_index_closest_to(100) == 3);
        CHECK(interp.get_y_index_closest_to(inf) == 3);
        // NaN must give a valid index. Which one is not part of the contract.
        int j = interp.get_y_index_closest_to(std::nan(""));
        CHECK(j >= 0);
        CHECK(j <= 3);
    }
}
