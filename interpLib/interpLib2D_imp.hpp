/**
 * Reads data from x, y, and z arrays. Arrays should all be the same length with each element corresponding to a data point.
 * Basically, x[i], y[i], and z[i] should correspond to the first, second, and third columns of a gnuplot file.
 */
template<typename Real>
void SplineInterp2D<Real>::setData( std::vector<Real> &x, std::vector<Real> &y, std::vector<Real> &z )
{/*{{{*/

    //Each x block will be the same size, so we need to find what that size is.
    //
    //This will also be the number of the new y blocks
    int ii = 0;
    Real xlast = x[0];
    while( fabs(x[ii]-xlast) < 1e-40 )
    {
        this->d2y.push_back(y[ii]);
        xlast = x[ii];
        ii++;
    }

    //this will give the length of the new y blocks
    int b = x.size() / ii;

    //This makes a vector of each of the unique x points by grabbing the first line of each x block
    std::vector<Real> tmpx;
    for (int i = 0; i < b; ++i)
    {
        int a = i*ii;
        tmpx.push_back(x[a]);
    }

    //This generates the new z blocks as well as making sure there are the right number of x blocks
    for (int i = 0; i < this->d2y.size(); ++i)
    {
        std::vector<Real> tmpz;
        int a = i*ii;

        for (int j = 0; j < tmpx.size(); ++j)
        {
            int b = i + j*ii;
            tmpz.push_back(z[b]);
        }
        this->d2z.push_back(tmpz);
        this->d2x.push_back(tmpx);
        tmpz.clear();
    }

    //old method
    // This puts the input data into blocks similar to gnuplot 2D data/*{{{*/
    //Real xlast = x[0];
    //this->x.push_back(x[0]);
    //std::vector<Real> tmpy, tmpz;
    //for (size_t i = 0; i < x.size(); ++i)
    //{
        //if( x[i] != xlast )
        //{
            //this->x.push_back(x[i]);
            //this->y.push_back(tmpy);
            //this->z.push_back(tmpz);

            //tmpy.clear();
            //tmpz.clear();
        //}

        //tmpy.push_back(y[i]);
        //tmpz.push_back(z[i]);

        //if( i == x.size() - 1 )
        //{
            //this->y.push_back(tmpy);
            //this->z.push_back(tmpz);
        //}

        //xlast = x[i];
    //}

    
    // This rearranges that data to being blocked out by the original second column.
    //
    // Ideally, this would just be done directly instead of blocking out in X first.
    //this->d2y = this->y[0];
    //std::vector<Real> tmpd2x, tmpd2z;
    //for (size_t i = 0; i < this->x.size(); ++i)
    //{
        //for (size_t j = 0; j < this->z[i].size(); ++j)
        //{
           //tmpd2x.push_back(this->x[j]);  
           //tmpd2z.push_back(this->z[j][i]);
        //}
        
        //this->d2x.push_back(tmpd2x);
        //this->d2z.push_back(tmpd2z);

        //tmpd2x.clear();
        //tmpd2z.clear();
    //}

    //this->x.clear();
    //this->y.clear();
    //this->z.clear();/*}}}*/
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




/**
 * Reads data from x, y, and z arrays. Arrays should all be the same length with each element corresponding to a data point.
 * Basically, x[i], y[i], and z[i] should correspond to the first, second, and third columns of a gnuplot file.
 *
 * WARNING: we a
 */
template<typename Real>
void BilinearInterp2D<Real>::setData( std::vector<Real> &_x, std::vector<Real> &_y, std::vector<Real> &_z )
{

  // we need to unpack the points in _x, _y, and _z
  // if all three vectors are the same length, then it means that together the specify (x,y,z) for a set of n points.
  // we need to split out the x and y coordinates and then map the z values onto a 2D array
  
  // figure out what the dimensions of the grid are first
  if( _x.size() == _y.size() && _x.size() == _z.size() )
  {
    // first we will determine the number of y coordinates by finding the index at
    // which the x value changes
    int Ny = 0;
    while( std::abs( _x[Ny] - _x[0] ) < std::numeric_limits<Real>::min() )
      Ny++;

    // Nz should be the size of _x, _y, and _z
    int Nz = _z.size();
    
    // Now we should have Nx * Ny == Nz
    // so Nx = Nz / Ny
    int Nx = Nz / Ny;

    this->X.resize( Nx );
    this->Y.resize( Ny );
    this->Z.resize( Nx, Ny );

    // grab x-coordinates. we are ASSUMING grid is regular
    this->X = VectorMap( _x.data(), Nx, InnerStrideType(Ny) );

    // grab y-coordinates. we are ASSUMING grid is regular
    this->Y = VectorMap( _y.data(), Ny, InnerStrideType(1) );

    // grab z values. we are ASSUMING grid is regular
    this->Z = MatrixMap(_z.data(), Nx, Ny, StrideType(Ny,1));

  }

  this->C = ArrayArray22( X.size()-1, Y.size()-1 );

  // we are going to precompute the interpolation coefficients so
  // that we can interpolate quickly
  for(int i = 0; i < X.size() - 1; i++)
  {
    for( int j = 0; j < Y.size() - 1; j++)
    {
      Real tmp = ( (X(i+1) - X(i) )*( Y(j+1) - Y(j) ) );
      this->C(i,j) = Z.block(i,j,2,2)/tmp;
    }
  }


}

template<typename Real>
Real BilinearInterp2D<Real>::operator()( Real _x, Real _y )
{

  // we want to interpolate to a point inside of a rectangle.
  // first, determine what element the point (_x,_y) is in.
  // the indices of the of bottom-left corner are the indecies of the element,
  // so we want to find the grid point to the left and below (_x,_y)


  // find the x index that is just to the LEFT of _x
  int i = 0;
  while( i < this->X.size()-2 && this->X(i+1) < _x )
      i++;

  // find the y index that is just BELOW _y
  int j = 0;
  while( j < this->Y.size()-2 && this->Y(j+1) < _y )
      j++;
  

  Array22 Q;  // coordinate matrix (see Wikipedia)
  Q << (this->X(i+1) -       _x  ) * (this->Y(j+1) - _y  ) ,  (this->X(i+1) -       _x  ) * (    _y - this->Y(j))
    ,  (          _x - this->X(i)) * (this->Y(j+1) - _y  ) ,  (          _x - this->X(i)) * (    _y - this->Y(j));

  return (Q*this->C(i,j)).sum();

}
