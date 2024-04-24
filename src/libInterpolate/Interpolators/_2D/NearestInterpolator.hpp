#ifndef Interpolators__2D_NearestInterpolator_hpp
#define Interpolators__2D_NearestInterpolator_hpp

/** @file NearestInterpolator.hpp
 * @brief
 * @author Finn Lukas Busch
 * @date 10/10/23
 */

#include <boost/range/adaptor/strided.hpp>
#include <boost/range/algorithm/lower_bound.hpp>

#include "InterpolatorBase.hpp"

/** @class
 * @brief Nearest interpolation for 2D functions.
 * @author Finn Lukas Busch
 */

namespace _2D {
template <class Real>
class NearestInterpolator : public InterpolatorBase<NearestInterpolator<Real>> {
   public:
    using BASE = InterpolatorBase<NearestInterpolator<Real>>;

   protected:
    using BASE::X;
    using BASE::xView;
    using BASE::Y;
    using BASE::yView;
    using BASE::Z;
    using BASE::zView;

   public:
    template <typename I>
    NearestInterpolator(I n, Real* x, Real* y, Real* z) {
        this->setData(n, x, y, z);
    }

    template <typename X, typename Y, typename Z>
    NearestInterpolator(X& x, Y& y, Z& z) {
        this->setData(x, y, z);
    }

    NearestInterpolator() : BASE() {}
    Real operator()(Real x, Real y) const;

    NearestInterpolator(const NearestInterpolator& rhs) : BASE(rhs) {}

    // copy-swap idiom
    friend void swap(NearestInterpolator& lhs, NearestInterpolator& rhs) {
        swap(static_cast<BASE&>(lhs), static_cast<BASE&>(rhs));
    }

    NearestInterpolator& operator=(NearestInterpolator rhs) {
        swap(*this, rhs);
        return *this;
    }
};

template <class Real>
Real NearestInterpolator<Real>::operator()(Real x, Real y) const {
    BASE::checkData();
    if (x < (*X)(0) || x > (*X)(X->size() - 1) || y < (*Y)(0) ||
        y > (*Y)(Y->size() - 1)) {
        return 0;
    }
    int i = this->get_x_index_closest_to(x);
    int j = this->get_y_index_closest_to(y);

    if (i < 0) i = 0;
    if (j < 0) j = 0;
    return Z->operator()(i, j);
}
}  // namespace _2D

#endif  // include protector
