# libInterp

A C++ interpolation library.

This library provides classes to perform various types of function interpolation (linear, spline, etc).

## Installing

Curretnly, `libInterp` is a header-only C++ library. To use it, simply include
the headers you want/need in your source code. If you use `git subrepo`, you
can clone the source into your externals directory and use it from there.

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
