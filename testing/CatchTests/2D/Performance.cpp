#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <libInterpolate/Interpolate.hpp>
using namespace Catch;
TEMPLATE_TEST_CASE("2D Interplation Benchmarks", "[.][benchmarks]",
                   _2D::BilinearInterpolator<double>,
                   _2D::BicubicInterpolator<double>) {
    int Nx = 2000;
    int Ny = 2000;
    std::vector<double> x(Nx * Ny), y(Nx * Ny), z(Nx * Ny);
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = (i / Ny) * 0.01;
        y[i] = (i % Ny) * 0.01;
        z[i] = x[i] * y[i];
    }

    TestType interp;
    interp.setData(x, y, z);

    BENCHMARK("Interpolate near a corner") { return interp(0.02, 0.02); };
    BENCHMARK("Interpolate near the middle") { return interp(10, 10); };
}
