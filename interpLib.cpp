#include <iostream>
#include <string>
#include <fstream>

#include "interpLib.h"
#include "fileReader/ReadFunction.h"
#include "boost/filesystem.hpp"
#include "matrixBuild.h"

namespace fs = boost::filesystem;

void SplineInterp::initialSetup()
{
    ublas::vector<double> rhs(this->n);
    ublas::permutation_matrix<int> P(this->n);
    //this function should do the stuff that is independant of which constructor is used
    this->A = matrixABuild(this->X, this->Y);
    this->B = matrixBBuild(this->X, this->Y);

    for (int i = 0; i < this->n; ++i)
    {
       rhs(i) = B(i); 
    }
    
    ublas::lu_factorize(this->A,P);
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

SplineInterp::SplineInterp(std::string dataFile)
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

            this->initialSetup();
            //this->A = matrixABuild(this->X, this->Y);
            //this->B = matrixBBuild(this->X, this->Y);

            //ublas::permutation_matrix<int> P(this->n);
            //ublas::vector<double> rhs(this->n);

            //for (size_t i = 0; i < this->B.size(); ++i)
            //{
               //rhs(i) = B(i); 
            //}

            //ublas::lu_factorize(this->A,P);
            //B = rhs;
            //ublas::lu_substitute(A,P,rhs);

            //this->a.resize( this->n - 1);
            //this->b.resize( this->n - 1);

            //for (int i = 0; i < this->n-1; ++i)
            //{
                //this->a(i) = rhs(i) * (X(i+1)-X(i)) - (Y(i+1) - Y(i));
                //this->b(i) = -rhs(i+1) * (X(i+1) - X(i)) + (Y(i+1) - Y(i));
            //}


            delete[] x;
            delete[] y;
        }
    }
}

SplineInterp::SplineInterp(std::vector<double> x, std::vector<double> y)
{
    //this should pretty much just map these into ublas::vectors
    this->n = x.size();

    this->X.resize(n);
    this->Y.resize(n);

    for (int i = 0; i < n; ++i)
    {
       this->X(i) = x[i]; 
       this->Y(i) = y[i]; 
    }

    this->initialSetup();
}

SplineInterp::~SplineInterp()
{
}

double SplineInterp::interpAt(double targetValue)
{
    double outVal;
    for (int i = 0; i < this->n; ++i)
    {
        //this is a feable check to see if you're passing in a value that's already in the data set.
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
