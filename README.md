# libInterp

A C++ library interpolation library.

This library provides classes to perform various types of function interpolation (linear, spline, etc).

## Installing

`libInterp` is a header-only C++ library. To use it, simply include the headers you want/need in your source code. If you use
`git subrepo`, you can clone the source into your externals directory and use it from there. Just add the `src/` directory
to your compiler includes. For example, if you are using the RHDO C++ component template:

    > git subrepo clone http://fermi.fhsu.edu/CD3/libInterp.git externals/libInterp 
    > vim CMakeLists.txt
        Add the following line:
        include_directories( ${PROJECT_SOURCE_DIR}/externals/libInterp/src )


