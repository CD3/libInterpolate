// Nonius - C++ benchmarking tool
//
// Written in 2014 by Martinho Fernandes <martinho.fernandes@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related
// and neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy of the CC0 Public Domain Dedication along with this software.
// If not, see <http://creativecommons.org/publicdomain/zero/1.0/>

// Measure

#ifndef NONIUS_DETAIL_MEASURE_HPP
#define NONIUS_DETAIL_MEASURE_HPP

#include <nonius/clock.h++>
#include <nonius/detail/complete_invoke.h++>
#include <nonius/detail/timing.h++>

#include <utility>

namespace nonius {
    namespace detail {
        template <typename Clock = default_clock, typename Fun, typename... Args>
        TimingOf<Clock, Fun(Args...)> measure(Fun&& fun, Args&&... args) {
            auto start = Clock::now();
            auto&& r = detail::complete_invoke(fun, std::forward<Args>(args)...);
            auto end = Clock::now();
            auto delta = end - start;
            return { delta, std::forward<decltype(r)>(r), 1 };
        }
    } // namespace detail
} // namespace nonius

#endif // NONIUS_DETAIL_MEASURE_HPP

