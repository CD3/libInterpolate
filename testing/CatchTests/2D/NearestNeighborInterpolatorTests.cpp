#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <fstream>
#include <libInterpolate/Interpolators/_2D/NearestNeighborInterpolator.hpp>
using namespace Catch;
namespace _2D {
class TestNearestNeighborInterp : public NearestNeighborInterpolator<double> {
   public:
    VectorType getX() { return *(this->X); }
    VectorType getY() { return *(this->Y); }
    MatrixType getZ() { return *(this->Z); }
};
}  // namespace _2D

TEST_CASE("NearestNeighborInterpolator Tests - Monotonic Data", "[bicubic]") {
    _2D::TestNearestNeighborInterp interp;

    //     0        1         2
    //
    // 0   x        x         x
    //
    // 1   x        x         x
    //
    // 2   x        x         x
    //
    //
    //
    //
    //
    // clang-format off
    _2D::NearestNeighborInterpolator<double>::VectorType xx(9), yy(9), zz(9);
    xx(0) = 0; yy(0) = 0; zz(0) = 0 + 0;
    xx(1) = 0; yy(1) = 1; zz(1) = 0 + 1;
    xx(2) = 0; yy(2) = 2; zz(2) = 0 + 2;
    xx(3) = 1; yy(3) = 0; zz(3) = 1 + 0;
    xx(4) = 1; yy(4) = 1; zz(4) = 1 + 1;
    xx(5) = 1; yy(5) = 2; zz(5) = 1 + 2;
    xx(6) = 2; yy(6) = 0; zz(6) = 2 + 0;
    xx(7) = 2; yy(7) = 1; zz(7) = 2 + 1;
    xx(8) = 2; yy(8) = 2; zz(8) = 2 + 2;
    // clang-format on

    interp.setData(xx, yy, zz);

    std::ofstream out;

    SECTION("Interpolation") {
        CHECK(interp(0, 0) == Approx(0).scale(1));
        CHECK(interp(0, 0.1) == Approx(0).scale(1));
        CHECK(interp(0, 0.49) == Approx(0).scale(1));
        CHECK(interp(0, 0.5) == Approx(1));
        CHECK(interp(0, 0.51) == Approx(1));
        CHECK(interp(0, 0.6) == Approx(1));
        CHECK(interp(0, 0.9) == Approx(1));
        CHECK(interp(0, 0.99) == Approx(1));
        CHECK(interp(0, 1.01) == Approx(1));
        CHECK(interp(0, 1.49) == Approx(1));
        CHECK(interp(0, 1.5) == Approx(2));
        CHECK(interp(0, 1.51) == Approx(2));
        CHECK(interp(0, 1.99) == Approx(2));
        CHECK(interp(0, 2) == Approx(2));
        CHECK(interp(0, 2.01) == Approx(0).scale(1));

        CHECK(interp(0, 0) == Approx(0).scale(1));
        CHECK(interp(0.1, 0) == Approx(0).scale(1));
        CHECK(interp(0.49, 0) == Approx(0).scale(1));
        CHECK(interp(0.5, 0) == Approx(1));
        CHECK(interp(0.51, 0) == Approx(1));
        CHECK(interp(0.6, 0) == Approx(1));
        CHECK(interp(0.9, 0) == Approx(1));
        CHECK(interp(0.99, 0) == Approx(1));
        CHECK(interp(1.01, 0) == Approx(1));
        CHECK(interp(1.49, 0) == Approx(1));
        CHECK(interp(1.5, 0) == Approx(2));
        CHECK(interp(1.51, 0) == Approx(2));
        CHECK(interp(1.99, 0) == Approx(2));
        CHECK(interp(2, 0) == Approx(2));
        CHECK(interp(2.01, 0) == Approx(0).scale(1));

        CHECK(interp(0.5, 0) == Approx(1));
        CHECK(interp(0.5, 0.1) == Approx(1));
        CHECK(interp(0.5, 0.49) == Approx(1));
        CHECK(interp(0.5, 0.5) == Approx(2));
        CHECK(interp(0.5, 0.51) == Approx(2));
        CHECK(interp(0.5, 0.6) == Approx(2));
        CHECK(interp(0.5, 0.9) == Approx(2));
        CHECK(interp(0.5, 0.99) == Approx(2));
        CHECK(interp(0.5, 1.01) == Approx(2));
        CHECK(interp(0.5, 1.49) == Approx(2));
        CHECK(interp(0.5, 1.5) == Approx(3));
        CHECK(interp(0.5, 1.51) == Approx(3));
        CHECK(interp(0.5, 1.99) == Approx(3));
        CHECK(interp(0.5, 2) == Approx(3));
        CHECK(interp(0.5, 2.01) == Approx(0).scale(1));

        CHECK(interp(1.01, 0) == Approx(1));
        CHECK(interp(1.01, 0.1) == Approx(1));
        CHECK(interp(1.01, 0.49) == Approx(1));
        CHECK(interp(1.01, 0.5) == Approx(2));
        CHECK(interp(1.01, 0.51) == Approx(2));
        CHECK(interp(1.01, 0.6) == Approx(2));
        CHECK(interp(1.01, 0.9) == Approx(2));
        CHECK(interp(1.01, 0.99) == Approx(2));
        CHECK(interp(1.01, 1.01) == Approx(2));
        CHECK(interp(1.01, 1.49) == Approx(2));
        CHECK(interp(1.01, 1.5) == Approx(3));
        CHECK(interp(1.01, 1.51) == Approx(3));
        CHECK(interp(1.01, 1.99) == Approx(3));
        CHECK(interp(1.01, 2) == Approx(3));
        CHECK(interp(1.01, 2.01) == Approx(0).scale(1));

        CHECK(interp(1.51, 0) == Approx(2));
        CHECK(interp(1.51, 0.1) == Approx(2));
        CHECK(interp(1.51, 0.49) == Approx(2));
        CHECK(interp(1.51, 0.5) == Approx(3));
        CHECK(interp(1.51, 0.51) == Approx(3));
        CHECK(interp(1.51, 0.6) == Approx(3));
        CHECK(interp(1.51, 0.9) == Approx(3));
        CHECK(interp(1.51, 0.99) == Approx(3));
        CHECK(interp(1.51, 1.01) == Approx(3));
        CHECK(interp(1.51, 1.49) == Approx(3));
        CHECK(interp(1.51, 1.5) == Approx(4));
        CHECK(interp(1.51, 1.51) == Approx(4));
        CHECK(interp(1.51, 1.99) == Approx(4));
        CHECK(interp(1.51, 2) == Approx(4));
        CHECK(interp(1.51, 2.01) == Approx(0).scale(1));
    }
}
