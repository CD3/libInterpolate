#pragma once

/**
MIT License

Copyright (c) 2018 Volodymyr Bilonenko

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

https://github.com/delfrrr/delaunator-cpp
 */

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace delaunator {

//@see
// https://stackoverflow.com/questions/33333363/built-in-mod-vs-custom-mod-function-improve-the-performance-of-modulus-op/33333636#33333636
inline size_t fast_mod(const size_t i, const size_t c) {
    return i >= c ? i % c : i;
}

// Kahan and Babuska summation, Neumaier variant; accumulates less FP error

template <typename T>
T sum(const std::vector<T>& x) {
    T sum = x[0];
    T err = 0.0;

    for (size_t i = 1; i < x.size(); i++) {
        const T k = x[i];
        const T m = sum + k;
        err += std::fabs(sum) >= std::fabs(k) ? sum - m + k : k - m + sum;
        sum = m;
    }
    return sum + err;
}

template <typename T>
T dist(const T ax, const T ay, const T bx, const T by) {
    const T dx = ax - bx;
    const T dy = ay - by;
    return dx * dx + dy * dy;
}

template <typename T>
T circumradius(const T ax, const T ay, const T bx, const T by, const T cx,
               const T cy) {
    const T dx = bx - ax;
    const T dy = by - ay;
    const T ex = cx - ax;
    const T ey = cy - ay;

    const T bl = dx * dx + dy * dy;
    const T cl = ex * ex + ey * ey;
    const T d = dx * ey - dy * ex;

    const T x = (ey * bl - dy * cl) * 0.5 / d;
    const T y = (dx * cl - ex * bl) * 0.5 / d;

    if ((bl > 0.0 || bl < 0.0) && (cl > 0.0 || cl < 0.0) &&
        (d > 0.0 || d < 0.0)) {
        return x * x + y * y;
    } else {
        return std::numeric_limits<T>::max();
    }
}

template <typename T>
bool orient(const T px, const T py, const T qx, const T qy, const T rx,
            const T ry) {
    return (qy - py) * (rx - qx) - (qx - px) * (ry - qy) < 0.0;
}
template <typename T>
std::pair<T, T> circumcenter(const T ax, const T ay, const T bx, const T by,
                             const T cx, const T cy) {
    const T dx = bx - ax;
    const T dy = by - ay;
    const T ex = cx - ax;
    const T ey = cy - ay;

    const T bl = dx * dx + dy * dy;
    const T cl = ex * ex + ey * ey;
    const T d = dx * ey - dy * ex;

    const T x = ax + (ey * bl - dy * cl) * 0.5 / d;
    const T y = ay + (dx * cl - ex * bl) * 0.5 / d;

    return std::make_pair(x, y);
}

template <typename T>
struct compare {
    std::vector<T> const& coords;
    T cx;
    T cy;

    bool operator()(std::size_t i, std::size_t j) {
        const T d1 = dist(coords[2 * i], coords[2 * i + 1], cx, cy);
        const T d2 = dist(coords[2 * j], coords[2 * j + 1], cx, cy);
        const T diff1 = d1 - d2;
        const T diff2 = coords[2 * i] - coords[2 * j];
        const T diff3 = coords[2 * i + 1] - coords[2 * j + 1];

        if (diff1 > 0.0 || diff1 < 0.0) {
            return diff1 < 0;
        } else if (diff2 > 0.0 || diff2 < 0.0) {
            return diff2 < 0;
        } else {
            return diff3 < 0;
        }
    }
};

template <typename T>
bool in_circle(const T ax, const T ay, const T bx, const T by, const T cx,
               const T cy, const T px, const T py) {
    const T dx = ax - px;
    const T dy = ay - py;
    const T ex = bx - px;
    const T ey = by - py;
    const T fx = cx - px;
    const T fy = cy - py;

    const T ap = dx * dx + dy * dy;
    const T bp = ex * ex + ey * ey;
    const T cp = fx * fx + fy * fy;

    return (dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) +
            ap * (ex * fy - ey * fx)) < 0.0;
}

constexpr std::size_t INVALID_INDEX = std::numeric_limits<std::size_t>::max();

template <typename T>
bool check_pts_equal(T x1, T y1, T x2, T y2) {
    return std::fabs(x1 - x2) <= std::numeric_limits<T>::epsilon() &&
           std::fabs(y1 - y2) <= std::numeric_limits<T>::epsilon();
}

// monotonically increases with real angle, but doesn't need expensive
// trigonometry
template <typename T>
T pseudo_angle(const T dx, const T dy) {
    const T p = dx / (std::abs(dx) + std::abs(dy));
    return (dy > 0.0 ? 3.0 - p : 1.0 + p) / 4.0;  // [0..1)
}

template <typename T>
struct DelaunatorPoint {
    std::size_t i;
    T x;
    T y;
    std::size_t t;
    std::size_t prev;
    std::size_t next;
    bool removed;
};

template <typename T>
class Delaunator {
   public:
    std::vector<T> const& coords;
    std::vector<std::size_t> triangles;
    std::vector<std::size_t> halfedges;
    std::vector<std::size_t> hull_prev;
    std::vector<std::size_t> hull_next;
    std::vector<std::size_t> hull_tri;
    std::size_t hull_start;

    Delaunator(std::vector<T> const& in_coords);

    T get_hull_area();

   private:
    std::vector<std::size_t> m_hash;
    T m_center_x;
    T m_center_y;
    std::size_t m_hash_size;
    std::vector<std::size_t> m_edge_stack;

    std::size_t legalize(std::size_t a);
    std::size_t hash_key(T x, T y) const;
    std::size_t add_triangle(std::size_t i0, std::size_t i1, std::size_t i2,
                             std::size_t a, std::size_t b, std::size_t c);
    void link(std::size_t a, std::size_t b);
};

template <typename T>
Delaunator<T>::Delaunator(std::vector<T> const& in_coords)
    : coords(in_coords),
      triangles(),
      halfedges(),
      hull_prev(),
      hull_next(),
      hull_tri(),
      hull_start(),
      m_hash(),
      m_center_x(),
      m_center_y(),
      m_hash_size(),
      m_edge_stack() {
    std::size_t n = coords.size() >> 1;

    T max_x = std::numeric_limits<T>::min();
    T max_y = std::numeric_limits<T>::min();
    T min_x = std::numeric_limits<T>::max();
    T min_y = std::numeric_limits<T>::max();
    std::vector<std::size_t> ids;
    ids.reserve(n);

    for (std::size_t i = 0; i < n; i++) {
        const T x = coords[2 * i];
        const T y = coords[2 * i + 1];

        if (x < min_x) min_x = x;
        if (y < min_y) min_y = y;
        if (x > max_x) max_x = x;
        if (y > max_y) max_y = y;

        ids.push_back(i);
    }
    const T cx = (min_x + max_x) / 2;
    const T cy = (min_y + max_y) / 2;
    T min_dist = std::numeric_limits<T>::max();

    std::size_t i0 = INVALID_INDEX;
    std::size_t i1 = INVALID_INDEX;
    std::size_t i2 = INVALID_INDEX;

    // pick a seed point close to the centroid
    for (std::size_t i = 0; i < n; i++) {
        const T d = dist(cx, cy, coords[2 * i], coords[2 * i + 1]);
        if (d < min_dist) {
            i0 = i;
            min_dist = d;
        }
    }

    const T i0x = coords[2 * i0];
    const T i0y = coords[2 * i0 + 1];

    min_dist = std::numeric_limits<T>::max();

    // find the point closest to the seed
    for (std::size_t i = 0; i < n; i++) {
        if (i == i0) continue;
        const T d = dist(i0x, i0y, coords[2 * i], coords[2 * i + 1]);
        if (d < min_dist && d > 0.0) {
            i1 = i;
            min_dist = d;
        }
    }

    T i1x = coords[2 * i1];
    T i1y = coords[2 * i1 + 1];

    T min_radius = std::numeric_limits<T>::max();

    // find the third point which forms the smallest circumcircle with the first
    // two
    for (std::size_t i = 0; i < n; i++) {
        if (i == i0 || i == i1) continue;

        const T r =
            circumradius(i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1]);

        if (r < min_radius) {
            i2 = i;
            min_radius = r;
        }
    }

    if (!(min_radius < std::numeric_limits<T>::max())) {
        throw std::runtime_error("not triangulation");
    }

    T i2x = coords[2 * i2];
    T i2y = coords[2 * i2 + 1];

    if (orient(i0x, i0y, i1x, i1y, i2x, i2y)) {
        std::swap(i1, i2);
        std::swap(i1x, i2x);
        std::swap(i1y, i2y);
    }

    std::tie(m_center_x, m_center_y) =
        circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);

    // sort the points by distance from the seed triangle circumcenter
    std::sort(ids.begin(), ids.end(),
              compare<T>{coords, m_center_x, m_center_y});

    // initialize a hash table for storing edges of the advancing convex hull
    m_hash_size =
        static_cast<std::size_t>(std::llround(std::ceil(std::sqrt(n))));
    m_hash.resize(m_hash_size);
    std::fill(m_hash.begin(), m_hash.end(), INVALID_INDEX);

    // initialize arrays for tracking the edges of the advancing convex hull
    hull_prev.resize(n);
    hull_next.resize(n);
    hull_tri.resize(n);

    hull_start = i0;

    size_t hull_size = 3;

    hull_next[i0] = hull_prev[i2] = i1;
    hull_next[i1] = hull_prev[i0] = i2;
    hull_next[i2] = hull_prev[i1] = i0;

    hull_tri[i0] = 0;
    hull_tri[i1] = 1;
    hull_tri[i2] = 2;

    m_hash[hash_key(i0x, i0y)] = i0;
    m_hash[hash_key(i1x, i1y)] = i1;
    m_hash[hash_key(i2x, i2y)] = i2;

    std::size_t max_triangles = n < 3 ? 1 : 2 * n - 5;
    triangles.reserve(max_triangles * 3);
    halfedges.reserve(max_triangles * 3);
    add_triangle(i0, i1, i2, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX);
    T xp = std::numeric_limits<T>::quiet_NaN();
    T yp = std::numeric_limits<T>::quiet_NaN();
    for (std::size_t k = 0; k < n; k++) {
        const std::size_t i = ids[k];
        const T x = coords[2 * i];
        const T y = coords[2 * i + 1];

        // skip near-duplicate points
        if (k > 0 && check_pts_equal(x, y, xp, yp)) continue;
        xp = x;
        yp = y;

        // skip seed triangle points
        if (check_pts_equal(x, y, i0x, i0y) ||
            check_pts_equal(x, y, i1x, i1y) || check_pts_equal(x, y, i2x, i2y))
            continue;

        // find a visible edge on the convex hull using edge hash
        std::size_t start = 0;

        size_t key = hash_key(x, y);
        for (size_t j = 0; j < m_hash_size; j++) {
            start = m_hash[fast_mod(key + j, m_hash_size)];
            if (start != INVALID_INDEX && start != hull_next[start]) break;
        }

        start = hull_prev[start];
        size_t e = start;
        size_t q;

        while (q = hull_next[e],
               !orient(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q],
                       coords[2 * q + 1])) {  // TODO: does it works in a same
                                              // way as in JS
            e = q;
            if (e == start) {
                e = INVALID_INDEX;
                break;
            }
        }

        if (e == INVALID_INDEX)
            continue;  // likely a near-duplicate point; skip it

        // add the first triangle from the point
        std::size_t t = add_triangle(e, i, hull_next[e], INVALID_INDEX,
                                     INVALID_INDEX, hull_tri[e]);

        hull_tri[i] = legalize(t + 2);
        hull_tri[e] = t;
        hull_size++;

        // walk forward through the hull, adding more triangles and flipping
        // recursively
        std::size_t next = hull_next[e];
        while (q = hull_next[next],
               orient(x, y, coords[2 * next], coords[2 * next + 1],
                      coords[2 * q], coords[2 * q + 1])) {
            t = add_triangle(next, i, q, hull_tri[i], INVALID_INDEX,
                             hull_tri[next]);
            hull_tri[i] = legalize(t + 2);
            hull_next[next] = next;  // mark as removed
            hull_size--;
            next = q;
        }

        // walk backward from the other side, adding more triangles and flipping
        if (e == start) {
            while (q = hull_prev[e],
                   orient(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e],
                          coords[2 * e + 1])) {
                t = add_triangle(q, i, e, INVALID_INDEX, hull_tri[e],
                                 hull_tri[q]);
                legalize(t + 2);
                hull_tri[q] = t;
                hull_next[e] = e;  // mark as removed
                hull_size--;
                e = q;
            }
        }

        // update the hull indices
        hull_prev[i] = e;
        hull_start = e;
        hull_prev[next] = i;
        hull_next[e] = i;
        hull_next[i] = next;

        m_hash[hash_key(x, y)] = i;
        m_hash[hash_key(coords[2 * e], coords[2 * e + 1])] = e;
    }
}

template <typename T>
T Delaunator<T>::get_hull_area() {
    std::vector<T> hull_area;
    size_t e = hull_start;
    do {
        hull_area.push_back((coords[2 * e] - coords[2 * hull_prev[e]]) *
                            (coords[2 * e + 1] + coords[2 * hull_prev[e] + 1]));
        e = hull_next[e];
    } while (e != hull_start);
    return sum(hull_area);
}

template <typename T>
std::size_t Delaunator<T>::legalize(std::size_t a) {
    std::size_t i = 0;
    std::size_t ar = 0;
    m_edge_stack.clear();

    // recursion eliminated with a fixed-size stack
    while (true) {
        const size_t b = halfedges[a];

        /* if the pair of triangles doesn't satisfy the Delaunay condition
         * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
         * then do the same check/flip recursively for the new pair of triangles
         *
         *           pl                    pl
         *          /||\                  /  \
         *       al/ || \bl            al/    \a
         *        /  ||  \              /      \
         *       /  a||b  \    flip    /___ar___\
         *     p0\   ||   /p1   =>   p0\---bl---/p1
         *        \  ||  /              \      /
         *       ar\ || /br             b\    /br
         *          \||/                  \  /
         *           pr                    pr
         */
        const size_t a0 = 3 * (a / 3);
        ar = a0 + (a + 2) % 3;

        if (b == INVALID_INDEX) {
            if (i > 0) {
                i--;
                a = m_edge_stack[i];
                continue;
            } else {
                // i = INVALID_INDEX;
                break;
            }
        }

        const size_t b0 = 3 * (b / 3);
        const size_t al = a0 + (a + 1) % 3;
        const size_t bl = b0 + (b + 2) % 3;

        const std::size_t p0 = triangles[ar];
        const std::size_t pr = triangles[a];
        const std::size_t pl = triangles[al];
        const std::size_t p1 = triangles[bl];

        const bool illegal =
            in_circle(coords[2 * p0], coords[2 * p0 + 1], coords[2 * pr],
                      coords[2 * pr + 1], coords[2 * pl], coords[2 * pl + 1],
                      coords[2 * p1], coords[2 * p1 + 1]);

        if (illegal) {
            triangles[a] = p1;
            triangles[b] = p0;

            auto hbl = halfedges[bl];

            // edge swapped on the other side of the hull (rare); fix the
            // halfedge reference
            if (hbl == INVALID_INDEX) {
                std::size_t e = hull_start;
                do {
                    if (hull_tri[e] == bl) {
                        hull_tri[e] = a;
                        break;
                    }
                    e = hull_next[e];
                } while (e != hull_start);
            }
            link(a, hbl);
            link(b, halfedges[ar]);
            link(ar, bl);
            std::size_t br = b0 + (b + 1) % 3;

            if (i < m_edge_stack.size()) {
                m_edge_stack[i] = br;
            } else {
                m_edge_stack.push_back(br);
            }
            i++;

        } else {
            if (i > 0) {
                i--;
                a = m_edge_stack[i];
                continue;
            } else {
                break;
            }
        }
    }
    return ar;
}

template <typename T>
std::size_t Delaunator<T>::hash_key(const T x, const T y) const {
    const T dx = x - m_center_x;
    const T dy = y - m_center_y;
    return fast_mod(static_cast<std::size_t>(std::llround(std::floor(
                        pseudo_angle(dx, dy) * static_cast<T>(m_hash_size)))),
                    m_hash_size);
}

template <typename T>
std::size_t Delaunator<T>::add_triangle(std::size_t i0, std::size_t i1,
                                        std::size_t i2, std::size_t a,
                                        std::size_t b, std::size_t c) {
    std::size_t t = triangles.size();
    triangles.push_back(i0);
    triangles.push_back(i1);
    triangles.push_back(i2);
    link(t, a);
    link(t + 1, b);
    link(t + 2, c);
    return t;
}

template <typename T>
void Delaunator<T>::link(const std::size_t a, const std::size_t b) {
    std::size_t s = halfedges.size();
    if (a == s) {
        halfedges.push_back(b);
    } else if (a < s) {
        halfedges[a] = b;
    } else {
        throw std::runtime_error("Cannot link edge");
    }
    if (b != INVALID_INDEX) {
        std::size_t s2 = halfedges.size();
        if (b == s2) {
            halfedges.push_back(a);
        } else if (b < s2) {
            halfedges[b] = a;
        } else {
            throw std::runtime_error("Cannot link edge");
        }
    }
}

}  // namespace delaunator
