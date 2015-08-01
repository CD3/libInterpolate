## Synopsis

This is a spline interpolator based on the method outlined on the [Wikipedia Page](https://en.wikipedia.org/wiki/Spline_interpolation). It currently supports reading code in from a file, or just interpolating from two std::vector<double> type arrays. Dr.Clark's multicolumn file reader is used for reading data files.

## 1D Interpolator Example

Initializing with two STL vectors:

	//populate these, obviously
	std::vector<double> x, y;

	SplineInterp<double> testInterp;
	testInterp.setData(x,y);
	
	double interpolatedValue = testInterp(pointToInterpolate);

The setData member function can be called multiple times and will functionally make a new interpolator object.

There's a short example for the 1D case in the example directory.

## Member functions (1D)
| **Function** | **Description** |
|----------|-------------|
| operator()(Real _x) | return interpolated value |
| integral(Real _a, Real _b) | Return the integral of the given range |
| integral() | Return the integral of the entire range |
| derivative(Real x) | Return the derivative at point x |
| setData(size_t n, Real *_x, Real *_y) | Set data using pointer arrays |
| setData(std::vector<Real> x, std::vector<Real> y) | Set data using std::vectors |

## 2D Interpolation Method

The 2D interpolator accepts 2D gnuplot formatted data similar to:

	0	0	0
	0	1	2
	0	2	4

	1	0	1
	1	1	3
	1	2	5

The data should be input as three STL vectors or double pointers. The code is only tested for data in which every block is the same lengths (each x point has the same number of y points).

The method separates the data into blocks as follows:

	std::vector<Real> x;
	std::vector< std::vector<Real> > y, z;

So x[0] is the x point for the first block, and y[0] and z[0] are the second two vectors of the first block.

The data then reorders the data so that the second column becomes the first column:

	std::vector<Real> y;
	std::vector< std::vector<Real> > x,z;

Each x,z vector pair is interpolated at the target x point. These interpolated  values are put into a vector, which would correspond to the 3rd column in the original data ordering. Because of the assumption mentioned earlier, the y vector in the second ordering is assumed to be the proper second column for the new block. The y vector from the second ordering and the newly interpolated z vector are interpolated for the targeted y value, returning the final value.

## 2D Interpolator Example

	//populate these
	std::vector<double> x, y, z;

	SplineInterp2D<double> testInterp;
	testInterp.setData(x,y,z);
	
	double interpolatedValue = testInterp(xPoint, yPoint);

## Member functions (2D)
| **Function** | **Description** |
|----------|-------------|
| operator()(Real _x, Real _y) | return interpolated value |
| setData(size_t n, Real \*_x, Real \*_y, Real \*_z) | Set data using pointer arrays |
| setData(std::vector<Real> x, std::vector<Real> y, std::vector<Real> z) | Set data using std::vectors |