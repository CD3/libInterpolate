#ifndef Interpolators__2D_NearestNeighborInterpolator_hpp
#define Interpolators__2D_NearestNeighborInterpolator_hpp

/** @file NearestNeighborInterpolator.hpp
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
class NearestNeighborInterpolator
    : public InterpolatorBase<NearestNeighborInterpolator<Real>> {
   public:
    using BASE = InterpolatorBase<NearestNeighborInterpolator<Real>>;

   protected:
    using BASE::X;
    using BASE::xView;
    using BASE::Y;
    using BASE::yView;
    using BASE::Z;
    using BASE::zView;

   public:
    template <typename I>
    NearestNeighborInterpolator(I n, Real* x, Real* y, Real* z) {
        this->setData(n, x, y, z);
    }

    template <typename X, typename Y, typename Z>
    NearestNeighborInterpolator(X& x, Y& y, Z& z) {
        this->setData(x, y, z);
    }

    NearestNeighborInterpolator() : BASE() {}
    Real operator()(Real x, Real y) const;

    NearestNeighborInterpolator(const NearestNeighborInterpolator& rhs)
        : BASE(rhs) {}

    // copy-swap idiom
    friend void swap(NearestNeighborInterpolator& lhs,
                     NearestNeighborInterpolator& rhs) {
        swap(static_cast<BASE&>(lhs), static_cast<BASE&>(rhs));
    }

    NearestNeighborInterpolator& operator=(NearestNeighborInterpolator rhs) {
        swap(*this, rhs);
        return *this;
    }
};

template <class Real>
Real NearestNeighborInterpolator<Real>::operator()(Real x, Real y) const {
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
