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



template<typename Real>
void BilinearInterp2D<Real>::setData( std::vector<Real> &_x, std::vector<Real> &_y, std::vector<Real> &_z )
{
}

/**
 * Reads data from x, y, and z arrays. Arrays should all be the same length with each element corresponding to a data point.
 * Basically, x[i], y[i], and z[i] should correspond to the first, second, and third columns of a gnuplot file.
 *
 * WARNING: we are assuming that data is on a regular grid and not a scatter plot
 */
template<typename Real>
void BilinearInterp2D<Real>::setData( size_t _n, Real *_x, Real *_y, Real *_z )
{

  // we need to unpack the points in _x, _y, and _z
  // if all three vectors are the same length, then it means that together the specify (x,y,z) for a set of n points.
  // we need to split out the x and y coordinates and then map the z values onto a 2D array
  
  // figure out what the dimensions of the grid are first

  // first we will determine the number of y coordinates by finding the index at
  // which the x value changes
  int Ny = 0;
  while( Ny < _n && std::abs( _x[Ny] - _x[0] ) < std::numeric_limits<Real>::min() )
    Ny++;

  // Nz should be the size of _x, _y, and _z
  int Nz = _n;
  
  // Now we should have Nx * Ny == Nz
  // so Nx = Nz / Ny
  int Nx = Nz / Ny;

  this->X.resize( Nx );
  this->Y.resize( Ny );
  this->Z.resize( Nx, Ny );

  // grab x-coordinates. we are ASSUMING grid is regular
  this->X = VectorMap( _x, Nx, InnerStrideType(Ny) );

  // grab y-coordinates. we are ASSUMING grid is regular
  this->Y = VectorMap( _y, Ny, InnerStrideType(1) );

  // grab z values. we are ASSUMING grid is regular
  this->Z = MatrixMap( _z, Nx, Ny, StrideType(Ny,1));


  // Interpolation will be done by multiplying the coordinates by coefficients.

  this->C = Matrix22Array( X.size()-1, Y.size()-1 );

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
  

  // first, create the coordinate vectors (see Wikipedia https://en.wikipedia.org/wiki/Bilinear_interpolation)
  RowVector2 vx;
  ColVector2 vy;
  vx << (this->X(i+1) - _x), (_x - this->X(i));
  vy << (this->Y(j+1) - _y), (_y - this->Y(j));

  // interpolation will now just be x*Q*y

  return vx*this->C(i,j)*vy;

}

template<typename Real>
Real BilinearInterp2D<Real>::integral( Real _xa, Real _xb, Real _ya, Real _yb )
{
  // No extrapolation
  int nx = X.size();
  int ny = Y.size();
  _xa = std::max( _xa, X(0) );
  _xb = std::min( _xb, X(nx-1) );
  _ya = std::max( _ya, Y(0) );
  _yb = std::min( _yb, Y(ny-1) );

  // find bottom-left and top-right corners
  // ia and ib will be to LEFT of _xa and _xb
  // ja and jb will be to BELOW _ya and _yb

  // bottom-left corner
  int ia = 0;
  while( ia < nx-2 && this->X(ia+1) < _xa )
      ia++;
  int ja = 0;
  while( ja < ny-2 && this->Y(ja+1) < _ya )
      ja++;

  // top-right corner
  int ib = 0;
  while( ib < nx-2 && this->X(ib+1) < _xb )
      ib++;
  int jb = 0;
  while( jb < ny-2 && this->Y(jb+1) < _yb )
      jb++;

  /*
   * We can integrate the function directly from the interpolation polynomial.
   * In Matrix form the polynomial looks like this
   *
   * /               \ /             \ /       \
   * | x_2-x   x-x_1 | | F_11   F_12 | | y_2-x |
   * \               / |             | |       |
   *                   | F_21   F_22 | | y-y_1 |
   *                   \             / \       /
   *
   * We can integrate this matrix equation to get the
   * matrix equation for the integral.
   *
   * /                                   \ /             \ /                 \
   * | x_2 x - 0.5 x^2   0.5 x^2 - x_1 x | | F_11   F_12 | | y_2 y - 0.5 y^2 |
   * \                                   / |             | |                 |
   *                                       | F_21   F_22 | | 0.5 y^2 - y_1 y |
   *                                       \             / \                 /
   * 
   * This must be evaluated at the limits for x and y. Note that if the limits on x are [a,b],
   * and the limits on y are [c,d], then to evaluate a function f(x,y) at the limits
   *
   *        |xb |xb
   * f(x,y) |   |   =  ( f(xb,yb) - f(xb,ya) ) - ( f(xa,yb) - f(xa,ya) )
   *        |xa |xa
   */

   // integrate the whole elements first...
   
   RowVector2 vxa,vxb;
   ColVector2 vya,vyb;
   Real xa,xb,ya,yb;
   Real sum = 0;
   int i,j;
   for(i = ia+1; i < ib; i++)
   {
     xa = X(i);
     xb = X(i+1);
     vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
     vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;
     for(j = ja+1; j < jb; j++)
     {
       ya = Y(j);
       yb = Y(j+1);
       vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
       vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

       sum += vxb*this->C(i,j)*vyb;
       sum -= vxb*this->C(i,j)*vya;
       sum -= vxa*this->C(i,j)*vyb;
       sum += vxa*this->C(i,j)*vya;
     }
   }

   // The double loop above will take care of all whole elements,
   // now we need to add the partial elements along the edge.
   //
   // However, there are some special cases we need to consider.
   // If we are integrating over the domain between the points a, b, c, and d
   // there is 1 whole element, and 8 partial elements.
   //
   // +------+-------------+---------------+-----+
   // |   c  |             |        d      |     |
   // |   +  |             |        +      |     |
   // |      |             |               |     |
   // +------+-------------+---------------+-----+
   // |      |             |               |     |
   // |      |             |               |     |
   // |      |             |               |     |
   // +------+-------------+---------------+-----+
   // |   +  |             |        +      |     |
   // |   a  |             |        b      |     |
   // +------+-------------+---------------+-----+
   //
   // 
   // We can integrate the function over this domain by:
   //
   // 1. integrate over all whole elements
   //    - limits of integration for each element will be:
   //      - lower x: x coordinate of the left edge of the element
   //      - upper x: x coordinate of the right edge of the element
   //      - lower y: y coordinate of bottom edge of the element
   //      - upper y: y coordinate of top edge of the element
   // 2. integrate over the bottom-left corner
   //    - limits of integration for will be:
   //      - lower x: x coordinate of point a
   //      - upper x: x coordinate of the right edge of the element
   //      - lower y: y coordinate of point a
   //      - upper y: y coordinate of top edge of the element
   // 3. integrate over the bottom-right corner
   //    - limits of integration for will be:
   //      - lower x: x coordinate of left edge of the element
   //      - upper x: x coordinate of point b
   //      - lower y: y coordinate of point b
   //      - upper y: y coordinate of top edge of the element
   // 4. integrate over the top-left corner
   //    - limits of integration for will be:
   //      - lower x: x coordinate of point c
   //      - upper x: x coordinate of the right edge of the element
   //      - lower y: y coordinate of the bottom edge of the element
   //      - upper y: y coordinate of point c
   // 5. integrate over the top-right corner
   //    - limits of integration for will be:
   //      - lower x: x coordinate of left edge of the element
   //      - upper x: x coordinate of point d
   //      - lower y: y coordinate of bottom edge of the elment
   //      - upper y: y coordinate of point d
   // 6. integrate over the left side elements
   //    - limits of integration for each element will be:
   //      - lower x: x coordinate of point a (same as x coordinate of point c)
   //      - upper x: x coordinate of right edge of the element
   //      - lower y: y coordinate of bottom edge of the elment
   //      - upper y: y coordinate of top edge of the element
   // 7. integrate over the right side elements
   //    - limits of integration for each element will be:
   //      - lower x: x coordinate of left edge of the element
   //      - upper x: x coordinate of point b (same as x coordinate of point d)
   //      - lower y: y coordinate of bottom edge of the elment
   //      - upper y: y coordinate of top edge of the element
   // 8. integrate over the bottom side elements
   //    - limits of integration for each element will be:
   //      - lower x: x coordinate of left edge of the element
   //      - upper x: x coordinate of right edge of the element
   //      - lower y: y coordinate of point a (same as y coordinate of point b)
   //      - upper y: y coordinate of top edge of the element
   // 9. integrate over the top side elements
   //    - limits of integration for each element will be:
   //      - lower x: x coordinate of left edge of the element
   //      - upper x: x coordinate of right edge of the element
   //      - lower y: y coordinate of bottom edge of the elment
   //      - upper y: y coordinate of point c (same as y coordinate of point d)
   //
   //
   // As long as the points a, b, c, and d are all in different elements, the above
   // algorithm will work. Steps 2 - 5 will alwasy be necessary, but steps 1, 6-9 will not be necessary
   // if two or more points are in ajacent elements. For example,
   // +------+-------------+---------------+-----+
   // |      |       c     |        d      |     |
   // |      |       +     |        +      |     |
   // |      |             |               |     |
   // +------+-------------+---------------+-----+
   // |      |             |               |     |
   // |      |             |               |     |
   // |      |             |               |     |
   // +------+-------------+---------------+-----+
   // |      |       +     |        +      |     |
   // |      |       a     |        b      |     |
   // +------+-------------+---------------+-----+
   // in this case, steps 1, 8, and 9 would not be required. The logic for determining when these steps
   // are required will automatically be handled by the for loops.
   //
   // However, if two or more points are in the same elements, we need to be careful. Consider:
   // +------+-------------+---------------+-----+
   // |      |             |    c   d      |     |
   // |      |             |    +   +      |     |
   // |      |             |               |     |
   // +------+-------------+---------------+-----+
   // |      |             |               |     |
   // |      |             |               |     |
   // |      |             |               |     |
   // +------+-------------+---------------+-----+
   // |      |             |    +   +      |     |
   // |      |             |    a   b      |     |
   // +------+-------------+---------------+-----+
   //
   // The method described above will not work because it will include the area between the left side
   // of each element and d AND the area between c and the right side of the element. Not only would it
   // be including an area that should not be included (between the left side of the element and c, and between
   // the right side of the element and d) it would double count the area that should be included (the area between c and d)
   //
   // Consider the middle element first (the one that does not contain any of the points a - d) 
   // Let the area between the left side and c be A_lc, the area between c and d be A_cd, and the area between d and the right
   // side be A_dr. The area that we want then is A_cd. The method above will give us A = A_lc + A_cd + A_cd + A_dr. If we note that
   // A_lc + A_cd + A_dr = A_E, the area of the element, then A_cd can be written as A_cd = A - A_E. So, we can subtract the integral
   // over the entire element from the integral that was computed above to get the correct integral.
   //
   // Subtracting the element area will work for the interior elements, but not for the elements containing the corners. In the top element
   // of the example above, we need to subtract the area of the element between the bottom edge and c instead of the entire element.
   // For the bottom element, we need to subtract the area between the top edge and a.
   //
   // Now consider an element that contains all 4 elements.
   // +------+-------------+---------------+-----+
   // |      |             |    c   d      |     |
   // |      |             |    +   +      |     |
   // |      |             |    +   +      |     |
   // |      |             |    a   b      |     |
   // +------+-------------+---------------+-----+
   // |      |             |               |     |
   // |      |             |               |     |
   // +------+-------------+---------------+-----+
   // |      |             |               |     |
   // |      |             |               |     |
   // +------+-------------+---------------+-----+
   //
   // In this case, the method above will only require steps 2 - 5.
   // The four points defined 9 different areas:
   // +-+-+-+
   // |1|2|3|
   // +-+-+-+
   // |4|5|6|
   // +-+-+-+
   // |7|8|9|
   // +-+-+-+
   //
   // Step 2 (bottom-left corner) will A2, A3, A5, and A6
   // Step 3 (bottom-right corner) will A1, A2, A4, and A5
   // Step 4 (top-left corner) will A5, A6, A8, and A9
   // Step 5 (top-left corner) will A4, A5, A7, and A8
   //
   // So, the total area covered will be
   // A = A1 + 2*A2 + A3 + 2*A4 + 4*A5 + 2*A6 + A7 + 2*A8 + A9 = AT + A2 + A4 + 3*A5 + A6 + A8
   //
   // The area we need is A5,
   //
   // A5 = (A - AT - A2 - A4 - A6 - A8)/3
   //
   // The cases with only two corners in an element would be different, so we probably are better
   // off just handling the corner cases directly.
   //
   // So, after all this...
   
   // add contributions from left side (excluding corners)
   i = ia;
   xa = _xa;
   xb = X(i+1);
   vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
   vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;
   for(j = ja+1; j < jb; j++)
   {
     ya = Y(j);
     yb = Y(j+1);
     vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
     vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

     sum += vxb*this->C(i,j)*vyb;
     sum -= vxb*this->C(i,j)*vya;
     sum -= vxa*this->C(i,j)*vyb;
     sum += vxa*this->C(i,j)*vya;
   }

   // add contributions from right side (excluding corners)
   i = ib;
   xa = X(i);
   xb = _xb;
   vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
   vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;
   for(j = ja+1; j < jb; j++)
   {
     ya = Y(j);
     yb = Y(j+1);
     vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
     vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

     sum += vxb*this->C(i,j)*vyb;
     sum -= vxb*this->C(i,j)*vya;
     sum -= vxa*this->C(i,j)*vyb;
     sum += vxa*this->C(i,j)*vya;
   }

   if( ia == ib )
   {
     // the min and max limits on x are in the same element.
     // we need to subtract the integral over the entire element (see above)
     i  = ia;
     xa = X(i);
     xb = X(i+1);
     vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
     vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;
     for(j = ja+1; j < jb; j++)
     {
       ya = Y(j);
       yb = Y(j+1);
       vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
       vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

       sum += -vxb*this->C(i,j)*vyb;
       sum -= -vxb*this->C(i,j)*vya;
       sum -= -vxa*this->C(i,j)*vyb;
       sum += -vxa*this->C(i,j)*vya;
     }
   }


   // add contributions from bottom side (excluding corners)
   j = ja;
   ya = _ya;
   yb = Y(j+1);
   vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
   vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

   for(i = ia+1; i < ib; i++)
   {
     xa = X(i);
     xb = X(i+1);
     vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
     vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;
     sum += vxb*this->C(i,j)*vyb;
     sum -= vxb*this->C(i,j)*vya;
     sum -= vxa*this->C(i,j)*vyb;
     sum += vxa*this->C(i,j)*vya;
   }

   // add contributions from top side (excluding corners)
   j = jb;
   ya = Y(j);
   yb = _yb;
   vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
   vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

   for(i = ia+1; i < ib; i++)
   {
     xa = X(i);
     xb = X(i+1);
     vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
     vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;
     sum += vxb*this->C(i,j)*vyb;
     sum -= vxb*this->C(i,j)*vya;
     sum -= vxa*this->C(i,j)*vyb;
     sum += vxa*this->C(i,j)*vya;
   }

   if( ja == jb )
   {
     // the min and max limits on y are in the same element.
     // we need to subtract the integral over the entire element (see above)
     j = ja;
     ya = Y(j);
     yb = Y(j+1);
     vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
     vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

     for(i = ia+1; i < ib; i++)
     {
       xa = X(i);
       xb = X(i+1);
       vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
       vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;
       sum += -vxb*this->C(i,j)*vyb;
       sum -= -vxb*this->C(i,j)*vya;
       sum -= -vxa*this->C(i,j)*vyb;
       sum += -vxa*this->C(i,j)*vya;
     }
   }





   if( ia != ib && ja != jb )
   {

   // bottom-left corner
   i = ia;
   xa = _xa;
   xb = X(i+1);
   vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
   vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;

   j = ja;
   ya = _ya;
   yb = Y(j+1);
   vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
   vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

   sum += vxb*this->C(i,j)*vyb;
   sum -= vxb*this->C(i,j)*vya;
   sum -= vxa*this->C(i,j)*vyb;
   sum += vxa*this->C(i,j)*vya;
   

   // bottom-right corner
   i = ib;
   xa = X(i);
   xb = _xb;
   vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
   vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;

   j = ja;
   ya = _ya;
   yb = Y(j+1);
   vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
   vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

   sum += vxb*this->C(i,j)*vyb;
   sum -= vxb*this->C(i,j)*vya;
   sum -= vxa*this->C(i,j)*vyb;
   sum += vxa*this->C(i,j)*vya;


   // top-left corner
   i = ia;
   xa = _xa;
   xb = X(i+1);
   vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
   vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;

   j = jb;
   ya = Y(j);
   yb = _yb;
   vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
   vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

   sum += vxb*this->C(i,j)*vyb;
   sum -= vxb*this->C(i,j)*vya;
   sum -= vxa*this->C(i,j)*vyb;
   sum += vxa*this->C(i,j)*vya;


   // top-right corner
   i = ib;
   xa = X(i);
   xb = _xb;
   vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
   vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;

   j = jb;
   ya = Y(j);
   yb = _yb;
   vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
   vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

   sum += vxb*this->C(i,j)*vyb;
   sum -= vxb*this->C(i,j)*vya;
   sum -= vxa*this->C(i,j)*vyb;
   sum += vxa*this->C(i,j)*vya;

   }

   if( ia == ib )
   {
     xa = _xa;
     xb = _xb;
     vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
     vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;

     // bottom block
     j = ja;
     ya = _ya;
     yb = Y(j+1);
     vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
     vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

     sum += vxb*this->C(i,j)*vyb;
     sum -= vxb*this->C(i,j)*vya;
     sum -= vxa*this->C(i,j)*vyb;
     sum += vxa*this->C(i,j)*vya;

     // top block
     j = jb;
     ya = Y(j);
     yb = _yb;
     vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
     vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

     sum += vxb*this->C(i,j)*vyb;
     sum -= vxb*this->C(i,j)*vya;
     sum -= vxa*this->C(i,j)*vyb;
     sum += vxa*this->C(i,j)*vya;
   }

   if( ja == jb )
   {
     ya = _ya;
     yb = _yb;
     vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
     vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

     // left block
     i = ia;
     xa = _xa;
     xb = X(i+1);
     vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
     vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;

     sum += vxb*this->C(i,j)*vyb;
     sum -= vxb*this->C(i,j)*vya;
     sum -= vxa*this->C(i,j)*vyb;
     sum += vxa*this->C(i,j)*vya;

     // right block
     i = ia;
     xa = _xa;
     xb = X(i+1);
     vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
     vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;

     sum += vxb*this->C(i,j)*vyb;
     sum -= vxb*this->C(i,j)*vya;
     sum -= vxa*this->C(i,j)*vyb;
     sum += vxa*this->C(i,j)*vya;

   }

   if( ia == ib && ja == jb )
   {
     xa = _xa;
     xb = _xb;
     vxa << this->X(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - this->X(i)*xa;
     vxb << this->X(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - this->X(i)*xb;

     ya = _ya;
     yb = _yb;
     vya << this->Y(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - this->Y(j)*ya;
     vyb << this->Y(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - this->Y(j)*yb;

     sum += vxb*this->C(i,j)*vyb;
     sum -= vxb*this->C(i,j)*vya;
     sum -= vxa*this->C(i,j)*vyb;
     sum += vxa*this->C(i,j)*vya;
   }






   return sum;
}
