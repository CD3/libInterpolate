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

`libInterp` depends on `Boost` and `Eigen3`, so you will need to include the directories
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

`libInterp` uses inheritance for code reuse and implements the "Curiously Recurring Template Patter" (CRTP).
CRTP allows the base class to implement functions that depend on parts that need to be implemented by the
derived class. For example, `InterpolatorBase` implements a function named `setData` that reads the
interpolated data into the interpolator. This way, derived classes do not need to implement functions to load the
interpolated data. However, different interpolation methods may require some additional setup. For example, they
may compute and store some coefficients that are used during interpolation. The `setData` function therefore calls a
function named `setupInterpolator` that is implemented by the derived class, if needed. This is where the CRTP comes in.

If you write a new interpolation method and derive from the `InterpolatorBase` class, you need to pass the derived class
to the base class as a template parameter.
```C++
class MyInterpolator : _1D::InterpolatorBase<MyInterpolator>
{
...
};
```
If your interpolator needs additional setup *after* the interpolated data has been read in, then you should also
implement the function `void setupInterpolator()`.
```C++
class MyInterpolator : _1D::InterpolatorBase<MyInterpolator>
{
  double operator()(double x); // implement the interpolation
  void setupInterpolator();    // implement for additional setup
};
```

CRTP does not allow runtime polymorphism, since each derived type derives from a different base class.
If you need a way to bind the interpolation type at runtime, you should use `std::function` instead

```C++
// OK. now interp can be used like a function and it uses "value semantics".
std::function<double(double)> interp = _1D::CubicSplineInterpolator<double>();

...

```

If you need runtime polymorphism with access to more member functions, you can use type erasure.
