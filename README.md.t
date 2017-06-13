# libInterp

A C++ interpolation library.

This library provides classes to perform various types of function interpolation (linear, spline, etc).

Currently implemented methods are:

- 1D Interpolation
    - linear (https://en.wikipedia.org/wiki/Linear_interpolation)
    - cubic spline (https://en.wikipedia.org/wiki/Spline_interpolation)
    - monotonic spline interpolation (http://adsabs.harvard.edu/full/1990A%26A...239..443S)
- 2D Interpolation
    - bilinear (https://en.wikipedia.org/wiki/Bilinear_interpolation)
    - bilcubic (https://en.wikipedia.org/wiki/Bicubic_interpolation)
    - thin plate spline (https://en.wikipedia.org/wiki/Thin_plate_spline)

## Example

```C++
#include <Interp.hpp>
...
vector<double> x,y;
...
fill x and y with data
...

// Use a cubic spline interpolator to interpolate the data
_1D::CubicSplineInterpolator<double> interp();
interp.setData(x,y);

double val = interp(2.0); // val contains the value of the function y(x) interpolated at x = 2.0

```

## Installing

Currently, `libInterp` is a header-only C++ library. To use it, simply include
the headers you want/need in your source code. If you use `git subrepo`, you
can clone the source into your externals directory and use it from there.

`libInterp` depends on Boost and Eigen3, so you will need to include the directories
containing their header files when compiling.

This simplest way to use the library is to build your project using CMake. You can then
put a copy of this project in a directory named `externals/libInterp` and include it in
your `CMakeLists.txt` file with the `add_subdirectory` command

```CMake
# add libInterp
add_subdirectory(externals/libInterp)
```

The `libInterp` `CMakeLists.txt` will check for its dependencies and export the required include
directories to a cache variable named `libInterp_INCLUDE_DIRS`, so you can just include them rather
than needing to search yourself.

```CMake
# add libInterp
add_subdirectory(externals/libInterp)
# add include dirs required by libInterp
include_directories( ${libInterp_INCLUDE_DIRS} )
```

## Design

`libInterp` aims to provides a simple, flexible, easy to use, interpolation library. 
The library defines an interpolator interface (the `InterpolatorInterface` abstract base
classes), and all interpolation algorithms are implemented as classes that
inherit this interface.  The benefit of this approach is that interpolation
methods can be easily interchanged, even at run-time. The potential downside is that
it requires virtual methods and run-time polymorphism, which *could* have an impact
on performance. In many applications, the flexibility to switch interpolation methods
greatly outweighs any possible performance issues.

## Contributing

To add an interpolation algorithm, simply create a class that inherits from the
`InterpolationInterface` class. You should put your in a file named
after the class name in the `src/Interpolators/_1D/` (or
`src/Interpolators/_2D/` if it is a 2 dimensional interpolator). Your class must
implement the following methods:

`setData(size_t n, Real* x, Real* y)` :
    This function takes the data points that are to be interpolated and sets up the
    interpolator.

`Real operator()(Real x)` :
    Interpolates the data to the point x.

`void operator()(size_t n, Real *x, Real *y)` :
    Interpolates the data to each of the points in `x`, writing the result to `y`.
