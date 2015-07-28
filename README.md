## Synopsis

This is a spline interpolator based on the method outlined on the [Wikipedia Page](https://en.wikipedia.org/wiki/Spline_interpolation). It currently supports reading code in from a file, or just interpolating from two std::vector<double> type arrays. Dr.Clark's multicolumn file reader is used for reading data files.

## Code Example

The class can read data from a file name:

	#include "splineInterp/interpLib.h"
	
	//somehow define the path to the filename
	std::string path = "nameOfTwoColumnDataFile.txt"

	//construct the class
	SplineInterp testInterp;
	
	testInterp.setData(path);

	double interpolatedValue = testInterp(pointToInterpolate);

or with two std::vector<double> types:

	#include "splineInterp/interpLib.h"
	#include <vector>

	std::vector<double> x;
	std::vector<double> y;

	SplineInterp testInterp;
	testInterp.setData(x,y);
	
	double interpolatedValue = testInterp(pointToInterpolate);
The class supports calling the setData member function multiple times in order to change out the data set used by the interpolator.

There is a short program in the example folder meant to show usage of setting data and getting interpolated values. It generates its own dataset, and then uses that as the data to interpolate against.

## Member functions
| **Function** | **Description** |
|----------|-------------|
| operator() | return interpolated value |
| integral(Real _a, Real _b) | Return the integral of the given range |
| integral() | Return the integral of the entire range |
| derivative(Real x) | Return the derivative at point x |
| setData(size_t n, Real *_x, Real *_y) | Set data using pointer arrays |
| setData(std::vector<Real> x, std::vector<Real> y) | Set data using std::vectors |
## Tests

A unit test using the Boost test suite is included with the library.