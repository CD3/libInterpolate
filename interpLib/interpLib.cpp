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

double SplineInterp::derivative(double targetValue)
{
    //No extrapolation
    if( targetValue < X(0) )
        return 0;

    if( targetValue >= X(this->n-1) )
        return 0;

    // find the index that is just to the right of the targetValue
    int i = 1;
    while( i < this->n-1 && X(i) < targetValue )
        i++;

    //this should be the same t as in the regular interpolation case
    double tmp = ( targetValue - X(i-1) ) / ( X(i) - X(i-1) );

    double outVal = ( Y(i) - Y(i-1) )/( X(i)-X(i-1) ) + ( 1 - 2*tmp )*( a(i-1)*(1-tmp) + b(i-1)*tmp )/( X(i) - X(i-1))
       + tmp*(1-tmp)*(b(i-1)-a(i-1))/(X(i)-X(i-1)) ;

    return outVal;
}

double SplineInterp::operator() (double targetValue)
{
  // don't extrapolate at all
  if( targetValue < X(0) )
    return 0;

  if( targetValue >= X(this->n-1) )
    return 0;

  // find the index that is just to the right of the targetValue
  int i = 1;
  while( i < this->n-1 && X(i) < targetValue )
    i++;

  double tmp = ( targetValue - X(i-1) ) / ( X(i) - X(i-1) );
  double outVal = ( 1 - tmp ) * Y(i-1) + tmp * Y(i) + tmp*(1-tmp)*(a(i-1)*(1-tmp)+b(i-1)*tmp);

  return outVal;
}

double SplineInterp::operator[] (double targetValue)
{
    return this->operator()(targetValue);
}

SplineInterp::SplineInterp()
{
}

SplineInterp::~SplineInterp()
{
}
