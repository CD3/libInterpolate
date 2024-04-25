#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <fstream>
#include <libInterpolate/AnyInterpolator.hpp>
#include <libInterpolate/Interpolate.hpp>
#include <string>
#include <vector>

using namespace Catch;
unsigned int Factorial(unsigned int number) {
    return number <= 1 ? number : Factorial(number - 1) * number;
}

TEST_CASE("Runtime dispatch with AnyInterpolator") {
    std::string method;
    std::vector<double> x, y;

    x.push_back(0);
    x.push_back(1);
    x.push_back(2);

    y.push_back(10);
    y.push_back(20);
    y.push_back(30);

    method = "linear";

    _1D::AnyInterpolator<double> interp;

    // select the interpolation method on user input
    if (method == "linear") interp = _1D::LinearInterpolator<double>();
    if (method == "cubicspline")
        interp = _1D::CubicSplineInterpolator<double>();
    if (method == "monotonic") interp = _1D::MonotonicInterpolator<double>();

    interp.setData(x.size(), x.data(), y.data());

    // interpolation is done with the operator() method.
    double val = interp(1.5);

    CHECK(val == Approx(25));
}

TEST_CASE("Interpolate the perimeter of a polygon") {
    std::vector<double> x, y;  // x-y coordinates of corners
    x.push_back(-39.363024);
    y.push_back(52.872388);
    // 1
    x.push_back(-37.615071);
    y.push_back(119.453522);
    // 2
    x.push_back(-75.928);
    y.push_back(119.754988);
    // 3
    x.push_back(-75.93703);
    y.push_back(119.689077);
    // 4
    x.push_back(-69.968429);
    y.push_back(76.1259);
    // 5
    x.push_back(-58.196226);
    y.push_back(59.691912);
    // 6
    x.push_back(-44.978186);
    y.push_back(53.709062);
    // 7
    x.push_back(-39.363024);
    y.push_back(52.872388);

    // no create a parametric variable that we can use to interpolate
    std::vector<double> t;

    t.push_back(0);
    auto distance = [](double x1, double y1, double x2, double y2) {
        return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
    };
    for (size_t i = 1; i < x.size(); ++i) {
        t.push_back(t[i - 1] + distance(x[i - 1], y[i - 1], x[i], y[i]));
    }

    _1D::CubicSplineInterpolator<double> interp_x, interp_y;
    interp_x.setData(t, x);
    interp_y.setData(t, y);

    {
        std::ofstream out("polygon-raw.txt");
        for (size_t i = 0; i < x.size(); ++i) {
            out << x[i] << " " << y[i] << "\n";
        }
    }

    {
        std::ofstream out("polygon-interpolated.txt");
        int N = 100;
        double dt = t[t.size() - 1] / (N - 1);
        for (int i = 0; i < N; ++i) {
            out << interp_x(dt * i) << " " << interp_y(dt * i) << "\n";
        }
    }
}
