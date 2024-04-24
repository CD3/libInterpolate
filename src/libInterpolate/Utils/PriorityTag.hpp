#pragma once

namespace libInterpolate {
namespace Utils {
template <int N>
struct priority_tag : priority_tag<N - 1> {};
template <>
struct priority_tag<0> {};
}  // namespace Utils
}  // namespace libInterpolate
