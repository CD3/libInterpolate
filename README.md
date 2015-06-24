## Synopsis

This is a spline interpolator based on the method outlined on the [Wikipedia Page](https://en.wikipedia.org/wiki/Spline_interpolation). It currently supports reading code in from a file, or just interpolating from two std::vector<double> type arrays. Dr.Clark's multicolumn file reader is used for reading data files.

## Code Example

The class can be constructed using a filename as follows:

	#include "splineInterp/interpLib.h"
	
	//somehow define the path to the filename
	std::string path = "nameOfTwoColumnDataFile.txt"

	//construct the class
	SplineInterp testInterp(path);

	double interpolatedValue = testInterp.interpAt( targetInterpValue );

or with two std::vector<double> types:

	#include "splineInterp/interpLib.h"
	#include <vector>

	std::vector<double> x;
	std::vector<double> y;

	SplineInterp testInterp(x, y);

	double interpValue = testInterp.interpAt( targetInterpValue );
At the time of construction, the class builds all of the matrices and calculates coefficients relevant to the data set needed by the interpolation method. The expensive part happens once, and then getting interpolated values afterwards is quick.

## Installation

This library depends on boost being installed. Compilation will require linking boost_system and boost_filesystem because boost::filesystem is used to validate the existance of the file to be read in. Also, the main data types of the calculation are boost::matrix and boost::vector.

The main.cpp and test.txt included with this project are just meant to show an example of basic usage for the library. The text file contains a coarsely defined quadratic from -10 to 10. The program reads this file in and interpolates a value passed by the command line.

	$ make
	$ ./int 3.32
	11.0221

## Tests

To be populated...
