#ifndef Interpolators__2D_BicubicInterpolator_hpp
#define Interpolators__2D_BicubicInterpolator_hpp

/** @file BicubicInterpolator.hpp
 * @brief
 * @author C.D. Clark III
 * @date 12/27/16
 */

#include "InterpolatorBase.hpp"

/** @class
 * @brief Cubic spline interpolation for for 2D functions.
 * @author C.D. Clark III, Aaron Hoffman
 *
 * This class implements the bicubic spline interpolation method.
 * It is essentially the 2D equivalent of cubic splines for 1D.
 *
 */

namespace _2D {

template <class Real>
class BicubicInterpolator : public InterpolatorBase<BicubicInterpolator<Real>> {
   public:
    using BASE = InterpolatorBase<BicubicInterpolator<Real>>;
    using VectorType = typename BASE::VectorType;
    using MapType = typename BASE::MapType;

    // types used to view data as 2D coordinates
    using MatrixType = typename BASE::MatrixType;
    using _2DVectorView = typename BASE::_2DVectorView;
    using _2DMatrixView = typename BASE::_2DMatrixView;

    // types used for 4x4 matrix algebra
    using Matrix44 = Eigen::Matrix<Real, 4, 4>;
    using Matrix44Array =
        Eigen::Array<Matrix44, Eigen::Dynamic, Eigen::Dynamic>;
    using ColVector4 = Eigen::Matrix<Real, 4, 1>;
    using RowVector4 = Eigen::Matrix<Real, 1, 4>;

   protected:
    using BASE::X;
    using BASE::xView;
    using BASE::Y;
    using BASE::yView;
    using BASE::Z;
    using BASE::zView;

    Matrix44Array
        a;  // naming convention used by wikipedia article (see Wikipedia
            // https://en.wikipedia.org/wiki/Bicubic_interpolation)

   public:
    template <typename I>
    BicubicInterpolator(I n, Real* x, Real* y, Real* z) {
        this->setData(n, x, y, z);
    }

    template <typename X, typename Y, typename Z>
    BicubicInterpolator(X& x, Y& y, Z& z) {
        this->setData(x, y, z);
    }

    BicubicInterpolator() : BASE() {}

    BicubicInterpolator(const BicubicInterpolator& rhs) : BASE(rhs), a(rhs.a) {}

    // copy-swap idiom
    friend void swap(BicubicInterpolator& lhs, BicubicInterpolator& rhs) {
        lhs.a.swap(rhs.a);
        swap(static_cast<BASE&>(lhs), static_cast<BASE&>(rhs));
    }

    BicubicInterpolator& operator=(BicubicInterpolator rhs) {
        swap(*this, rhs);
        return *this;
    }

    // methods required by the interface
    Real operator()(Real x, Real y) const;

   protected:
    void setupInterpolator();
    friend BASE;
};

template <class Real>
void BicubicInterpolator<Real>::setupInterpolator() {
    // Interpolation will be done by multiplying the coordinates by
    // coefficients.

    a = Matrix44Array(X->size() - 1, Y->size() - 1);

    // We are going to precompute the interpolation coefficients so
    // that we can interpolate quickly. This requires a 4x4 matrix for each
    // "patch".
    Matrix44 Left, Right;

    Left << 1, 0, 0, 0, 0, 0, 1, 0, -3, 3, -2, -1, 2, -2, 1, 1;

    Right << 1, 0, -3, 2, 0, 0, 3, -2, 0, 1, -2, 1, 0, 0, -1, 1;

    for (int i = 0; i < X->size() - 1; i++) {
        for (int j = 0; j < Y->size() - 1; j++) {
            Matrix44 F;

            Real f00, f01, f10, f11;
            Real fx00, fx01, fx10, fx11;
            Real fy00, fy01, fy10, fy11;
            Real fxy00, fxy01, fxy10, fxy11;

            Real fm, fp;
            int im, ip, jm, jp;

            int iN = X->size();
            int jN = Y->size();

            // function values
            f00 = (*Z)(i, j);          // <<<<<<
            f01 = (*Z)(i, j + 1);      // <<<<<<
            f10 = (*Z)(i + 1, j);      // <<<<<<
            f11 = (*Z)(i + 1, j + 1);  // <<<<<<

            // need to calculate function values and derivatives
            // at each corner.
            //
            // note: interpolation algorithm is derived for the unit square.
            // so we need to take the derivatives assuming X(i+1) - X(i) =
            // Y(j+1) - Y(j) = 1
            //
            // using finite-difference approximation. central-difference on the
            // interior, forwared and backward difference on the ends.

            Real xL =
                (*X)(i + 1) - (*X)(i);  // actual side length in x direction
            Real yL =
                (*Y)(j + 1) - (*Y)(j);  // actual side length in y direction
            Real dx, dy;

            // x derivatives
            //
            //   (i,j+1)     (i+1,j+1)
            //
            //   +        +
            //
            //
            //   +        +
            //
            //  (i,j)     (i+1,j)

            im = std::max(i - 1, 0);  // get index to the "left" of current i
            ip = std::min(i + 1,
                          iN - 1);  // get index to the "right" of current i

            dx = ((*X)(ip) - (*X)(im)) /
                 xL;  // relative distance between "right" and "left" neighbors

            fp = (*Z)(ip, j);
            fm = (*Z)(im, j);
            fx00 = (fp - fm) / dx;  // <<<<<<

            fp = (*Z)(ip, j + 1);
            fm = (*Z)(im, j + 1);
            fx01 = (fp - fm) / dx;  // <<<<<<

            im = std::max(i, 0);
            ip = std::min(i + 2, iN - 1);

            dx = ((*X)(ip) - (*X)(im)) / xL;

            fp = (*Z)(ip, j);
            fm = (*Z)(im, j);
            fx10 = (fp - fm) / dx;  // <<<<<<

            fp = (*Z)(ip, j + 1);
            fm = (*Z)(im, j + 1);
            fx11 = (fp - fm) / dx;  // <<<<<<

            // y derivatives

            jm = std::max(j - 1, 0);
            jp = std::min(j + 1, jN - 1);

            dy = ((*Y)(jp) - (*Y)(jm)) / yL;

            fp = (*Z)(i, jp);
            fm = (*Z)(i, jm);
            fy00 = (fp - fm) / dy;  // <<<<<<

            fp = (*Z)(i + 1, jp);
            fm = (*Z)(i + 1, jm);
            fy10 = (fp - fm) / dy;  // <<<<<<

            jm = std::max(j, 0);
            jp = std::min(j + 2, jN - 1);

            dy = ((*Y)(jp) - (*Y)(jm)) / yL;

            fp = (*Z)(i, jp);
            fm = (*Z)(i, jm);
            fy01 = (fp - fm) / dx;  // <<<<<<

            fp = (*Z)(i + 1, jp);
            fm = (*Z)(i + 1, jm);
            fy11 = (fp - fm) / dx;  // <<<<<<

            // xy derivatives

            im = std::max(i - 1, 0);
            ip = std::min(i + 1, iN - 1);
            jm = std::max(j - 1, 0);
            jp = std::min(j + 1, jN - 1);

            dx = ((*X)(ip) - (*X)(im)) / xL;

            dy = ((*Y)(jp) - (*Y)(jm)) / yL;

            fp = ((*Z)(ip, jp) - (*Z)(im, jp)) / dx;
            fm = ((*Z)(ip, jm) - (*Z)(im, jm)) / dx;
            fxy00 = (fp - fm) / dy;  // <<<<<<

            jm = std::max(j, 0);
            jp = std::min(j + 2, jN - 1);

            dy = ((*Y)(jp) - (*Y)(jm)) / yL;

            fp = ((*Z)(ip, jp) - (*Z)(im, jp)) / dx;
            fm = ((*Z)(ip, jm) - (*Z)(im, jm)) / dx;
            fxy01 = (fp - fm) / dy;  // <<<<<<

            im = std::max(i, 0);
            ip = std::min(i + 2, iN - 1);
            jm = std::max(j - 1, 0);
            jp = std::min(j + 1, jN - 1);

            dx = ((*X)(ip) - (*X)(im)) / xL;

            dy = ((*Y)(jp) - (*Y)(jm)) / yL;

            fp = ((*Z)(ip, jp) - (*Z)(im, jp)) / dx;
            fm = ((*Z)(ip, jm) - (*Z)(im, jm)) / dx;
            fxy10 = (fp - fm) / dy;  // <<<<<<

            jm = std::max(j, 0);
            jp = std::min(j + 2, jN - 1);

            dy = ((*Y)(jp) - (*Y)(jm)) / yL;

            fp = ((*Z)(ip, jp) - (*Z)(im, jp)) / dx;
            fm = ((*Z)(ip, jm) - (*Z)(im, jm)) / dx;
            fxy11 = (fp - fm) / dy;  // <<<<<<

            F << f00, f01, fy00, fy01, f10, f11, fy10, fy11, fx00, fx01, fxy00,
                fxy01, fx10, fx11, fxy10, fxy11;

            a(i, j) = Left * F * Right;
        }
    }
}

template <class Real>
Real BicubicInterpolator<Real>::operator()(Real x, Real y) const {
    BASE::checkData();

    // no extrapolation...
    if (x < this->xView->minCoeff() || x > this->xView->maxCoeff() ||
        y < this->yView->minCoeff() || y > this->yView->maxCoeff()) {
        return 0;
    }

    int i = this->get_x_index_to_left_of(x);
    int j = this->get_y_index_below(y);

    if (i < 0) i = 0;
    if (j < 0) j = 0;

    Real xL = (*X)(i + 1) - (*X)(i);
    Real yL = (*Y)(j + 1) - (*Y)(j);

    // now, create the coordinate vectors (see Wikipedia
    // https://en.wikipedia.org/wiki/Bicubic_interpolation)
    RowVector4 vx;
    vx[0] = 1;                   // x^0
    vx[1] = (x - (*X)(i)) / xL;  // x^1
    vx[2] = vx[1] * vx[1];       // x^2
    vx[3] = vx[2] * vx[1];       // x^3

    ColVector4 vy;
    vy[0] = 1;                   // y^0
    vy[1] = (y - (*Y)(j)) / yL;  // y^1
    vy[2] = vy[1] * vy[1];       // y^2
    vy[3] = vy[2] * vy[1];       // y^3

    // interpolation is just x*a*y

    return vx * a(i, j) * vy;
}

}  // namespace _2D

#endif  // include protector
