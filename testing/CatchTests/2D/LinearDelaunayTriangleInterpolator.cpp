#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <libInterpolate/Interpolators/_2D/BilinearInterpolator.hpp>
#include <libInterpolate/Interpolators/_2D/LinearDelaunayTriangleInterpolator.hpp>
using namespace Catch;
TEMPLATE_TEST_CASE("2D LinearDelaunayTriangleInterpolator Tests", "[current]",
                   double, float) {
    _2D::LinearDelaunayTriangleInterpolator<TestType> interp;

    SECTION("Single Triangle") {
        std::vector<TestType> x(3), y(3), z(3);

        SECTION("Level") {
            // clang-format off
      x[0] = 0; y[0] = 0; z[0] = 1;
      x[1] = 1; y[1] = 0; z[1] = 1;
      x[2] = 0; y[2] = 1; z[2] = 1;
            // clang-format on

            interp.setData(x, y, z);

            CHECK(interp(0, 0) == Approx(1));
            CHECK(interp(0, 1) == Approx(1));
            CHECK(interp(1, 0) == Approx(1));
            CHECK(interp(0.1, 0.1) == Approx(1));

            CHECK(interp(0, -0.1) == Approx(0));
            CHECK(interp(1, 1) == Approx(0));
            CHECK(interp(-0.1, 0) == Approx(0));
        }

        SECTION("Angled") {
            // clang-format off
      x[0] = 0; y[0] = 1; z[0] = 0;
      x[1] = 0; y[1] =-1; z[1] = 0;
      x[2] = 1; y[2] = 0; z[2] = 1;
            // clang-format on

            interp.setData(x, y, z);

            CHECK(interp(0, 0) == Approx(0).scale(1));
            CHECK(interp(0.5, 0) == Approx(0.5));
            CHECK(interp(0.5, 0.1) == Approx(0.5));
            CHECK(interp(0.5, 1) == Approx(0));
            CHECK(interp(1, 0) == Approx(1));
        }
    }
    SECTION("Multiple Triangles") {
        std::vector<TestType> x(5), y(5), z(5);
        SECTION("Table Top") {
            // clang-format off
      x[0] = -1; y[0] = -1; z[0] = 3.14;
      x[1] = -1; y[1] =  1; z[1] = 3.14;
      x[2] =  1; y[2] = -1; z[2] = 3.14;
      x[3] =  1; y[3] =  1; z[3] = 3.14;
      x[4] =  0; y[4] =  0; z[4] = 3.14;
            // clang-format on

            interp.setData(x, y, z);

            CHECK(interp(0, 0) == Approx(3.14));
            CHECK(interp(2, 0) == Approx(0).scale(1));
            CHECK(interp(0, 2) == Approx(0).scale(1));
            CHECK(interp(0.5, 0.5) == Approx(3.14));
        }

        SECTION("Pyramid") {
            // clang-format off
      x[0] = -1; y[0] = -1; z[0] = 0;
      x[1] = -1; y[1] =  1; z[1] = 0;
      x[2] =  1; y[2] = -1; z[2] = 0;
      x[3] =  1; y[3] =  1; z[3] = 0;
      x[4] =  0; y[4] =  0; z[4] = 2;
            // clang-format on

            interp.setData(x, y, z);

            CHECK(interp(0, 0) == Approx(2));
            CHECK(interp(-1, 0) == Approx(0).scale(1));
            CHECK(interp(0, 1) == Approx(0).scale(1));

            CHECK(interp(0, 0.5) == Approx(1));
            CHECK(interp(0.5, 0) == Approx(1));
            CHECK(interp(0, -0.5) == Approx(1));
            CHECK(interp(-0.5, 0) == Approx(1));
            CHECK(interp(0.5, 0.5) == Approx(1));
        }
    }
}

TEMPLATE_TEST_CASE("2D LinearDelaunayTriangleInterpolator Benchmarks",
                   "[.][benchmarks][current]", double, float) {
    _2D::LinearDelaunayTriangleInterpolator<TestType> interp;
    SECTION("Benchmarks", "[.][benchmarks]") {
        _2D::BilinearInterpolator<TestType> bilin_interp;
        SECTION("Setup Interpolator") {
            size_t Nx = 100;
            size_t Ny = 100;
            size_t N = Nx * Ny;
            std::vector<TestType> x(N), y(N), z(N);

            for (size_t i = 0; i < Nx; ++i) {
                for (size_t j = 0; j < Ny; ++j) {
                    x[i * Ny + j] = i * 0.2;
                    y[i * Ny + j] = j * 0.4;
                    z[i * Ny + j] = i * j * 0.4;
                }
            }

            // BENCHMARK("setData(...)")
            // {
            //   interp.setData(x, y, z);
            // };
        }
        SECTION("Interpolating") {
            size_t Nx = 1000;
            size_t Ny = 1000;
            size_t N = Nx * Ny;
            std::vector<TestType> x(N), y(N), z(N);

            for (size_t i = 0; i < Nx; ++i) {
                for (size_t j = 0; j < Ny; ++j) {
                    x[i * Ny + j] = i * 0.2;
                    y[i * Ny + j] = j * 0.4;
                    z[i * Ny + j] = i * j * 0.4;
                }
            }
            interp.setData(x, y, z);
            bilin_interp.setData(x, y, z);

            BENCHMARK(
                "LinearDelaunayTriangleInterpolator::operator()(30.3,55.3)") {
                return interp(30.3, 55.5);
            };
            BENCHMARK("Bilinear::operator()(30.3,55.3)") {
                return bilin_interp(30.3, 55.5);
            };
        }
    }
}

TEMPLATE_TEST_CASE("2D LinearDelaunayTriangleInterpolator data setup",
                   "[delaunay]", double, float) {
    // z = x + 2y + 10 on a square with a center point, so the data is never 0
    // and a planar interpolator reproduces it exactly.
    std::vector<TestType> x{-1, -1, 1, 1, 0}, y{-1, 1, -1, 1, 0}, z(5);
    for (size_t i = 0; i < z.size(); ++i) z[i] = x[i] + 2 * y[i] + 10;

    SECTION("no data throws") {
        _2D::LinearDelaunayTriangleInterpolator<TestType> interp;
        CHECK_THROWS_AS(interp(0, 0), std::logic_error);
    }

    // setUnsafeDataReference used to segfault because the triangulation read
    // the interpolator's own (empty) copy of the data.
    SECTION("setUnsafeDataReference matches setData") {
        _2D::LinearDelaunayTriangleInterpolator<TestType> owned, ref;
        owned.setData(x, y, z);
        ref.setUnsafeDataReference(x.size(), x.data(), y.data(), z.data());

        for (TestType xx : {-1., -0.5, 0., 0.25, 1.}) {
            for (TestType yy : {-1., -0.75, 0., 0.5, 1.}) {
                CHECK(ref(xx, yy) == Approx(xx + 2 * yy + 10));
                CHECK(ref(xx, yy) == owned(xx, yy));
            }
        }
        // outside the convex hull
        CHECK(ref(2, 0) == 0);
        CHECK(ref(0, -2) == 0);
    }
}
