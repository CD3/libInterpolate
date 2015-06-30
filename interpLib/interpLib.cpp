#include "interpLib.h"
#include "fileReader/ReadFunction.h"
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

void SplineInterp::setData(std::string dataFile)
{
    if( fs::exists( dataFile ) )
    {
        if( fs::is_regular_file( dataFile ) )
        {
            std::ifstream inputFile(dataFile.c_str());

            double *x, *y;

            RUC::ReadFunction( inputFile, x, y, this->n );

            inputFile.close();

            this->X.resize(this->n);
            this->Y.resize(this->n);

            for ( int i = 0; i < this->n; ++i)
            {
                this->X(i) = x[i];
                this->Y(i) = y[i];
            }

            this->initCoefficients();

            delete[] x;
            delete[] y;
        }
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

SplineInterp::SplineInterp()
{
}

SplineInterp::~SplineInterp()
{
}

double SplineInterp::operator() (double targetValue)
{
    double outVal;
    for (int i = 0; i < this->n; ++i)
    {
        //this is a not very good check to see if you're passing in a value that's already in the data set.
        //this equality should be checked with a better function for comparing doubles
        if( X(i) == targetValue )
        {
            outVal = Y(i);
        } 

        else if(targetValue < X(i) && targetValue > X(i-1))
        {
           double tmp = ( targetValue - X(i-1) ) / ( X(i) - X(i-1) );
           double q = ( 1 - tmp ) * Y(i-1) + tmp * Y(i) + tmp*(1-tmp)*(a(i-1)*(1-tmp)+b(i-1)*tmp);

           outVal = q;
       } 
    }
 
    return outVal;
}

double SplineInterp::operator[] (double targetValue)
{
    double outVal;
    for (int i = 0; i < this->n; ++i)
    {
        //this is a not very good check to see if you're passing in a value that's already in the data set.
        //this equality should be checked with a better function for comparing doubles
        if( X(i) == targetValue )
        {
            outVal = Y(i);
        } 

        else if(targetValue < X(i) && targetValue > X(i-1))
        {
           double tmp = ( targetValue - X(i-1) ) / ( X(i) - X(i-1) );
           double q = ( 1 - tmp ) * Y(i-1) + tmp * Y(i) + tmp*(1-tmp)*(a(i-1)*(1-tmp)+b(i-1)*tmp);

           outVal = q;
       } 
    }
 
    return outVal;
}
