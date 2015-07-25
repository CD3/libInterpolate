#ifndef interplib_imp_hpp
#define interplib_imp_hpp

template<typename Real>
Real SplineInterp<Real>::operator()(Real x)
{/*{{{*/
  // don't extrapolate at all
  if( x < X(0) )
    return 0;
     
  if( x > X(this->n-1) )
    return 0;
  
  // find the index that is just to the right of the x
  int i = 1;
  while( i < this->n-1 && X(i) < x )
    i++;

  // See the wikipedia page on "Spline interpolation" (https://en.wikipedia.org/wiki/Spline_interpolation)
  // for a derivation this interpolation.
  Real t = ( x - X(i-1) ) / ( X(i) - X(i-1) );
  Real q = ( 1 - t ) * Y(i-1) + t * Y(i) + t*(1-t)*(a(i-1)*(1-t)+b(i-1)*t);
  
  return q;
} /*}}}*/

template<typename Real>
Real SplineInterp<Real>::operator[](Real x)
{/*{{{*/
  // don't extrapolate at all
  if( x < X(0) )
    return 0;
     
  if( x > X(this->n-1) )
    return 0;
  
  // find the index that is just to the right of the x
  int i = 1;
  while( i < this->n-1 && X(i) < x )
    i++;

  // See the wikipedia page on "Spline interpolation" (https://en.wikipedia.org/wiki/Spline_interpolation)
  // for a derivation this interpolation.
  Real t = ( x - X(i-1) ) / ( X(i) - X(i-1) );
  Real q = ( 1 - t ) * Y(i-1) + t * Y(i) + t*(1-t)*(a(i-1)*(1-t)+b(i-1)*t);
  
  return q;
}   /*}}}*/

template<typename Real>
Real SplineInterp<Real>::derivative(Real x)
{/*{{{*/
    //No extrapolation
    if( x < X(0) )
        return 0;

    if( x > X(this->n-1) )
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
}/*}}}*/

template<typename Real>
Real SplineInterp<Real>::integral(Real _a, Real _b)
{/*{{{*/
    // allow b to be less than a
    double sign = 1;
    if( _a > _b )
    {
      std::swap( a, b );
      sign = -1;
    }
    //No extrapolation
    _a = std::max( _a, X(0) );
    _b = std::min( _b, X(this->n-1) );

    // find the indexes that is just to the right of a and b
    int ai = 1;
    while( ai < this->n-1 && X(ai) < _a )
        ai++;
    int bi = 1;
    while( bi < this->n-1 && X(bi) < _b )
        bi++;

    /**
     *
     * We can integrate the function directly using its cubic spline representation.
     *
     * from wikipedia:
     *
     * q(x) = ( 1 - t )*y_1 + t*y_2 + t*( 1 - t )( a*(1 - t) + b*t )
     * 
     * t = (x - x_1) / (x_2 - x_1)
     *
     * I = \int_a^b q(x) dx = \int_a^b ( 1 - t )*y_1 + t*y_2 + t*( 1 - t )( a*(1 - t) + b*t ) dx
     *
     * variable substitution: x -> t
     *
     * dt = dx / (x_2 - x_1)
     * t_a = (a - x_1) / (x_2 - x_1)
     * t_b = (b - x_1) / (x_2 - x_1)
     *
     * I = (x_2 - x_1) \int_t_a^t_b ( 1 - t )*y_1 + t*y_2 + t*( 1 - t )( a*(1 - t) + b*t ) dt
     *
     *   = (x_2 - x_1) [ ( t - t^2/2 )*y_1 + t^2/2*y_2 + a*(t^2 - 2*t^3/3 + t^4/4) + b*(t^3/3 - t^4/4) ] |_t_a^t_b
     *
     * if we integrate over the entire element, i.e. x -> [x1,x2], then we will have
     * t_a = 0, t_b = 1. This gives
     *
     * I = (x_2 - x_1) [ ( 1 - 1/2 )*y_1 + 1/2*y_2 + a*(1 - 2/3 + 1/4) + b*(1/3 - 1/4) ]
     *
     */

    double x_1, x_2, t;
    double y_1, y_2;
    double sum = 0;
    for( int i = ai; i < bi-1; i++)
    {
      // x_1 -> X(i)
      // x_2 -> X(i+1)
      // y_1 -> Y(i)
      // y_2 -> Y(i+1)
      x_1 = X(i);
      x_2 = X(i+1);
      y_1 = Y(i);
      y_2 = Y(i+1);
      // X(ai) is to the RIGHT of _a
      // X(bi) is to the RIGHT of _b, but i only goes up to bi-2 and
      // X(bi-1) is to the LEFT of _b
      // therefore, we are just handling interior elements in this loop.
      sum += (x_2 - x_1)*( 0.5*(y_1 + y_2) + (1./12)*(a(i) + b(i)) );
    }


    // now we need to handle the area between [_a,X(ai)] and [X(bi-1),_b]


    // [X(0),_b]
    // x_1 -> X(bi-1)
    // x_2 -> X(bi)
    // y_1 -> Y(bi-1)
    // y_2 -> Y(bi)
    x_1 = X(bi-1);
    x_2 = X(bi);
    y_1 = Y(bi-1);
    y_2 = Y(bi);
    t   = (_b - x_1)/(x_2 - x_1);

    // adding area between x_1 and _b
    sum += (x_2 - x_1) * ( ( t - pow(t,2)/2 )*y_1 + pow(t,2)/2.*y_2 + a(bi-1)*(pow(t,2) - 2.*pow(t,3)/3. + pow(t,4)/4.) + b(bi-1)*(pow(t,3)/3. - pow(t,4)/4.) );

    //
    // [_a,X(0)]
    // x_1 -> X(ai-1)
    // x_2 -> X(ai)
    // y_1 -> Y(ai-1)
    // y_2 -> Y(ai)
    x_1 = X(ai-1);
    x_2 = X(ai);
    y_1 = Y(ai-1);
    y_2 = Y(ai);
    t   = (_a - x_1)/(x_2 - x_1);

    // subtracting area from x_1 to _a
    sum -= (x_2 - x_1) * ( ( t - pow(t,2)/2 )*y_1 + pow(t,2)/2.*y_2 + a(ai-1)*(pow(t,2) - 2.*pow(t,3)/3. + pow(t,4)/4.) + b(ai-1)*(pow(t,3)/3. - pow(t,4)/4.) );

    if( ai != bi ) // _a and _b are not in the in the same element, need to add area of element containing _a
      sum += (x_2 - x_1)*( 0.5*(y_1 + y_2) + (1./12)*(a(ai-1) + b(ai-1)) );

    return sign*sum;
}/*}}}*/

template<typename Real>
Real SplineInterp<Real>::integral()
{/*{{{*/
    return this->integral( X(0), X(this->n-1) );
}/*}}}*/

template<typename Real>
void SplineInterp<Real>::setData( size_t _n, Real *_x, Real *_y )
{/*{{{*/
    this->n = _n;

    this->X.resize(this->n);
    this->Y.resize(this->n);

    for (int i = 0; i < this->n; ++i)
    {
       this->X(i) = _x[i];
       this->Y(i) = _y[i];
    }

    this->initCoefficients();
}/*}}}*/

template<typename Real>
void SplineInterp<Real>::setData( std::vector<Real> x, std::vector<Real> y )
{/*{{{*/
    this->setData( x.size(), x.data(), y.data() );
}/*}}}*/

template<typename Real>
void SplineInterp<Real>::initCoefficients()
{

    /*
     * This now solves the problem Ax=B using the Thomas algorithm, because the matrix A will be tridiagonal.
     */

    //init the matrices that get solved by the lu factorization
    ublas::matrix<Real> A(this->n, this->n);
    ublas::vector<Real> B(this->n);

    //intermediate values for the solver
    ublas::vector<Real> rhs(this->n);
    ublas::permutation_matrix<int> P(this->n);

    //build the matrices that get solved
    A = matrixABuild<Real>(this->X);
    B = matrixBBuild(this->X, this->Y);

    /********************************** everything above this point is fine. ******************************************/

    std::vector<Real> f;
    f.resize( this->n );

    for (int i = 0; i < this->n; ++i)
    {
       f[i] = B(i);
    }

    for (size_t i = 0; i < A.size1(); ++i)
    {
        for (size_t j = 0; j < A.size2(); ++j)
        {
            std::cout << A(i,j) << "    ";
        }
        std::cout << std::endl;
    }

    //c is a vector of the upper diagonals of matrix A
    std::vector<Real> c;
    for (size_t i = 0; i < this->n-1; ++i)
    {/*{{{*/
        c.push_back( A(i,i+1) );
     
    }/*}}}*/

    std::cout << std::endl;
    //b is a vector of the diagnoals of matrix A
    std::vector<Real> b;
    for (size_t i = 0; i < this->n; ++i)
    {
        b.push_back(A(i,i));
    }


    //a is a vector of the lower diagonals of matrix A
    std::vector<Real> a;
    for (size_t i = 1; i < this->n; ++i)
    {
        a.push_back(A(i,i-1));
    }

    c[0] /= b[0];
    f[0] /= b[0];

    for (size_t i = 0; i < this->n; ++i)
    {
        c[i] /= b[i] - a[i]*c[i-1];
        f[i] = (f[i] - a[i]*f[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    int N = this->n - 1;
    f[N] = (f[N] - a[N]*f[N-1]) / (b[N] - a[N]*c[N-1]);

    for( int i = n; i-- > 0;)
    {
        f[i] -= c[i]*f[i+1];
    }



    /********************************** everything below this point is fine. ******************************************/
    this->a.resize(this->n - 1);
    this->b.resize(this->n - 1);

    for (int i = 0; i < this->n - 1; ++i)
    {
        this->a(i) = f[i] * (X(i+1)-X(i)) - (Y(i+1) - Y(i));
        this->b(i) = -f[i+1] * (X(i+1) - X(i)) + (Y(i+1) - Y(i));
    }


    for (int i = 0; i < this->a.size(); ++i)
    {
        std::cout << this->a(i) << std::endl;
    }
}

#endif
