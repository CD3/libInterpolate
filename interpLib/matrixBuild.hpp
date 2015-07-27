#include <array>
#include <iostream>

//this is also the first thing that needs to be converted to boost matrices. 
//each matrixABuild will have to be done first and tested. Get an old version to compare against

template<typename Real>
std::vector< std::vector<Real> > matrixABuild(std::vector<Real> x )
{
    const int N = x.size();

    std::vector< std::vector<Real> > A;
    A.resize(N);
    for (size_t i = 0; i < A.size(); ++i)
    {
        A[i].resize(N);
    }


    for(int i = 0; i < N; ++i)
    {
        if( i == 0 )
        {
            A[i][i] = 2 / (x[i+1] - x[i] );
            A[i][i+1] = 1 / (x[i+1] - x[i] );
            for(int j = i + 2; j < N; ++j)
            {
                A[i][j] = 0;
            }
        } else if ( i == N-1)
        {
            A[i][i-1] = 1 / (x[i] - x[i-1]);
            A[i][i] = 2 / (x[i] - x[i-1] );
            for(int j = 0; j < i-1; ++j)
            {
                A[i][j] = 0;
            }
        } else {
            for(int j = 0; j < i -1  ; ++j)
            {
                A[i][j] = 0;
            }
            for(int j = i + 1; j < N; ++j)
            {
                A[i][j] = 0;
            }
            A[i][i-1] = 1 / (x[i] - x[i-1]);
            A[i][i] = ( 1 / (x[i]-x[i-1]) + 1/(x[i+1] - x[i])) * 2;
            A[i][i+1] = 1/(x[i+1] - x[i]);
        }
    }

    return A;
}

template< typename Real>
std::vector<Real> matrixBBuild(std::vector<Real> x, std::vector<Real> y)
{
    const int n = x.size();
    std::vector<Real> _B;
    _B.resize(n);

    for(int i = 0; i < n; ++i)
    {   
        if(i == 0)
        {   
            _B[i] = 3 * ( y[i+1] - y[i] )/pow(x[i+1]-x[i],2);
        } else if( i == n-1 )
        {   
            _B[i] = 3 * (y[i] - y[i-1])/pow(x[i] - x[i-1],2);
        } else { 
            _B[i] = 3 * ( (y[i] - y[i-1])/(pow(x[i]-x[i-1],2)) + (y[i+1] - y[i])/pow(x[i+1] - x[i],2));     
        }
    }
    return _B;

}
