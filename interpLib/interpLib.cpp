#include "interpLib.h"
#include "matrixBuild.h"


void SplineInterp::initCoefficients()
{

    //init the matrices that get solved by the lu factorization
    ublas::matrix<double> A(this->n, this->n);
    ublas::vector<double> B(this->n);

    //intermediate values for the solver
    ublas::vector<double> rhs(this->n);
    ublas::permutation_matrix<int> P(this->n);

    //build the matrices that get solved
    A = matrixABuild(this->X, this->Y);
    B = matrixBBuild(this->X, this->Y);

    for (int i = 0; i < this->n; ++i)
    {
       rhs(i) = B(i); 
    }
    
    ublas::lu_factorize(A,P);
    B = rhs;
    ublas::lu_substitute(A,P,rhs);

    this->a.resize(this->n - 1);
    this->b.resize(this->n - 1);

    for (int i = 0; i < this->n - 1; ++i)
    {
        this->a(i) = rhs(i) * (X(i+1)-X(i)) - (Y(i+1) - Y(i));
        this->b(i) = -rhs(i+1) * (X(i+1) - X(i)) + (Y(i+1) - Y(i));
    }

}

void SplineInterp::setData(size_t _n, double* _x, double* _y)
{
    this->n = _n;

    this->X.resize(this->n);
    this->Y.resize(this->n);

    for (int i = 0; i < this->n; ++i)
    {
       this->X(i) = _x[i]; 
       this->Y(i) = _y[i]; 
    }

    this->initCoefficients();
}
  
void SplineInterp::setData(std::vector<double> x, std::vector<double> y)
{
  this->setData( x.size(), x.data(), y.data() );
}

double SplineInterp::derivative(double x)
{
    //No extrapolation
    if( x < X(0) )
        return 0;

    if( x >= X(this->n-1) )
        return 0;

    // find the index that is just to the right of x
    int i = 1;
    while( i < this->n-1 && X(i) < x )
        i++;

    //this should be the same t as in the regular interpolation case
    double t = ( x - X(i-1) ) / ( X(i) - X(i-1) );

    double qprime = ( Y(i) - Y(i-1) )/( X(i)-X(i-1) ) + ( 1 - 2*t )*( a(i-1)*(1-t) + b(i-1)*t )/( X(i) - X(i-1))
                    + t*(1-t)*(b(i-1)-a(i-1))/(X(i)-X(i-1)) ;

    return qprime;
}

double SplineInterp::operator() (double x)
{
  // don't extrapolate at all
  if( x < X(0) )
    return 0;

  if( x >= X(this->n-1) )
    return 0;

  // find the index that is just to the right of the x
  int i = 1;
  while( i < this->n-1 && X(i) < x )
    i++;

  // See the wikipedia page on "Spline interpolation" (https://en.wikipedia.org/wiki/Spline_interpolation)
  // for a derivation this interpolation.
  double t = ( x - X(i-1) ) / ( X(i) - X(i-1) );
  double q = ( 1 - t ) * Y(i-1) + t * Y(i) + t*(1-t)*(a(i-1)*(1-t)+b(i-1)*t);

  return q;
}

double SplineInterp::operator[] (double x)
{
    return this->operator()(x);
}

SplineInterp::SplineInterp()
{
}

SplineInterp::~SplineInterp()
{
}
