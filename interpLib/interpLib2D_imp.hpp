template<typename Real>
void SplineInterp2D<Real>::setData( std::vector<Real> &x, std::vector<Real> &y, std::vector<Real> &z )
{/*{{{*/

    // This puts the input data into blocks similar to gnuplot 2D data
    Real xlast = x[0];
    this->x.push_back(x[0]);
    std::vector<Real> tmpy, tmpz;
    for (size_t i = 0; i < x.size(); ++i)
    {
        if( x[i] != xlast )
        {
            this->x.push_back(x[i]);
            this->y.push_back(tmpy);
            this->z.push_back(tmpz);

            tmpy.clear();
            tmpz.clear();
        }

        tmpy.push_back(y[i]);
        tmpz.push_back(z[i]);

        if( i == x.size() - 1 )
        {
            this->y.push_back(tmpy);
            this->z.push_back(tmpz);
        }

        xlast = x[i];
    }
    
    // This rearranges that data to being blocked out by the original second column.
    //
    // Ideally, this would just be done directly instead of blocking out in X first.
    this->d2y = this->y[0];
    std::vector<Real> tmpd2x, tmpd2z;
    for (size_t i = 0; i < this->x.size(); ++i)
    {
        for (size_t j = 0; j < this->z[i].size(); ++j)
        {
           tmpd2x.push_back(this->x[j]);  
           tmpd2z.push_back(this->z[j][i]);
        }
        
        this->d2x.push_back(tmpd2x);
        this->d2z.push_back(tmpd2z);

        tmpd2x.clear();
        tmpd2z.clear();
    }

    this->x.clear();
    this->y.clear();
    this->z.clear();
}/*}}}*/

template<typename Real>
Real SplineInterp2D<Real>::operator()( Real _x, Real _y )
{/*{{{*/

    SplineInterp<Real> interp1D;

    //interpolate in x for each of the y blocks
    std::vector<Real> newz;
    for (size_t i = 0; i < this->d2x.size(); ++i)
    {
        interp1D.setData(d2x[i], d2z[i]);
        newz.push_back( interp1D(_x) );
    }

    //interpolate the new block at the y point
    interp1D.setData( this->d2y, newz );
    return interp1D(_y);


}/*}}}*/
