#ifndef Interpolators__2D_BilinearInterpolator_hpp
#define Interpolators__2D_BilinearInterpolator_hpp

/** @file BilinearInterpolator.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 12/27/16
  */

#include "InterpolatorBase.hpp"

/** @class 
  * @brief Linear interpolation for 2D functions.
  * @author C.D. Clark III, Aaron Hoffman
  */

namespace _2D {

template<class Real>
class BilinearInterpolator : public InterpolatorBase<Real>
{
  public:
    // typedefs
    typedef typename InterpolatorBase<Real>::VectorType VectorType;
    typedef typename InterpolatorBase<Real>::MapType MapType;
    typedef typename InterpolatorBase<Real>::GradientType GradientType;

    // types used to view data as 2D coordinates
    typedef typename Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> MatrixType;
    typedef Eigen::Map<VectorType,Eigen::Unaligned,Eigen::InnerStride<Eigen::Dynamic               >> _2DVectorView;
    typedef Eigen::Map<MatrixType,Eigen::Unaligned,     Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>> _2DMatrixView;

    // types used for 2x2 matrix algebra
    typedef Eigen::Matrix<Real,2,2 > Matrix22;
    typedef Eigen::Array< Matrix22, Eigen::Dynamic, Eigen::Dynamic > Matrix22Array;
    typedef Eigen::Matrix<Real,2,1 > ColVector2;
    typedef Eigen::Matrix<Real,1,2 > RowVector2;

    // methods required by the interface
    virtual Real operator()( Real x, Real y ) const;
    virtual Real integral( Real a, Real b, Real c, Real d ) const;
    virtual GradientType gradient( Real x, Real y ) const;
    virtual void setData( size_t _n, Real *x, Real *y, Real *z, bool deep_copy = true );
    using InterpolatorBase<Real>::setData;

  protected:
    using InterpolatorBase<Real>::xv;
    using InterpolatorBase<Real>::yv;
    using InterpolatorBase<Real>::zv;
    // these maps are used to view the x,y,z data as two coordinate vectors and a function matrix, instead of three vectors.
    std::shared_ptr<_2DVectorView> X,Y;
    std::shared_ptr<_2DMatrixView> Z;
    
    Matrix22Array Q; // naming convention used by wikipedia article (see Wikipedia https://en.wikipedia.org/wiki/Bilinear_interpolation)


};

template<class Real>
void
BilinearInterpolator<Real>::setData( size_t n, Real *x, Real *y, Real *z, bool deep_copy )
{
  InterpolatorBase<Real>::setData(n,x,y,z,deep_copy);

  // setup 2D view of the data
  // We need to figure out what the x and y dimensions are.
  int Nx = 0, Ny = 0;
  // Ny will be the number of elements that have the same x coordinate
  Real xlast = (*xv)(0);
  while( Ny < n-1 && fabs((*xv)(Ny)-xlast) < 1e-40 )
    Ny++;
  Nx = n/Ny;

  // consecutive values in the x data are separated by Ny, so this is the inner stride for X
  X.reset( new _2DVectorView( &(*xv)(0), Nx, Eigen::InnerStride<Eigen::Dynamic>(Ny) ) );

  // consecutive values in the y data are next to each other, so the stride is just 1
  Y.reset( new _2DVectorView( &(*yv)(0), Ny, Eigen::InnerStride<Eigen::Dynamic>(1) ) );

  // Eigen defaults to COLUMN MAJOR
  // consecutive elements in a column are separated by Ny (this is the inner stride)
  // consecutive elements in a row are located next to each other (this is the outer stride)
  // Stride object takes outer,inner as arguments.
  Z.reset( new _2DMatrixView( &(*zv)(0), Nx, Ny, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(1,Ny) ) );

  // Interpolation will be done by multiplying the coordinates by coefficients.

  Q = Matrix22Array( X->size()-1, Y->size()-1 );

  // We are going to pre-compute the interpolation coefficients so
  // that we can interpolate quickly
  for(int i = 0; i < X->size() - 1; i++)
  {
    for( int j = 0; j < Y->size() - 1; j++)
    {
      Real tmp = ( ((*X)(i+1) - (*X)(i) )*( (*Y)(j+1) - (*Y)(j) ) );
      Q(i,j) = Z->block(i,j,2,2)/tmp;
    }
  }
}

template<class Real>
Real
BilinearInterpolator<Real>::operator()( Real x, Real y ) const
{
  InterpolatorBase<Real>::checkData();
  
  // no extrapolation...
  if( x < (*X)(0)
   || x > (*X)(X->size()-1)
   || y < (*Y)(0)
   || y > (*Y)(Y->size()-1) )
  {
    return 0;
  }
  // find the x index that is just to the LEFT of x
  int i  = Utils::index_last_lt( x, *X, X->size() );
  if(i < 0)
    i = 0;

  // find the y index that is just BELOW y
  int j  = Utils::index_last_lt( y, *Y, Y->size() );
  if(j < 0)
    j = 0;
  


  // now, create the coordinate vectors (see Wikipedia https://en.wikipedia.org/wiki/Bilinear_interpolation)
  RowVector2 vx;
  ColVector2 vy;
  vx << ((*X)(i+1) - x), (x - (*X)(i));
  vy << ((*Y)(j+1) - y), (y - (*Y)(j));

  // interpolation is just x*Q*y

  return vx*Q(i,j)*vy;
}

template<class Real>
auto
BilinearInterpolator<Real>::gradient( Real x, Real y ) const -> GradientType
{
  // no extrapolation...
  if( x < (*X)(0)
   || x > (*X)(X->size()-1)
   || y < (*Y)(0)
   || y > (*Y)(Y->size()-1) )
  {
    return {0,0};
  }

  // find the x index that is just to the LEFT of x
  int i  = Utils::index_last_lt( x, *X, X->size() );
  if(i < 0)
    i = 0;

  // find the y index that is just BELOW y
  int j  = Utils::index_last_lt( y, *Y, Y->size() );
  if(j < 0)
    j = 0;


  // now, create the coordinate vectors (see Wikipedia https://en.wikipedia.org/wiki/Bilinear_interpolation)
  RowVector2 vx;
  ColVector2 vy;

  Real dzdx, dzdy;

  vx << -1,1;
  vy << ((*Y)(j+1) - y), (y - (*Y)(j));
  dzdx = vx*Q(i,j)*vy;

  vx << ((*X)(i+1) - x), (x - (*X)(i));
  vy << -1,1;
  dzdy = vx*Q(i,j)*vy;

  return {dzdx,dzdy};

}

template<class Real>
Real
BilinearInterpolator<Real>::integral( Real _xa, Real _xb, Real _ya, Real _yb ) const
{
  // No extrapolation
  int nx = X->size();
  int ny = Y->size();
  _xa = std::max( _xa, (*X)(0) );
  _xb = std::min( _xb, (*X)(nx-1) );
  _ya = std::max( _ya, (*Y)(0) );
  _yb = std::min( _yb, (*Y)(ny-1) );

  // find bottom-left and top-right corners
  // ia and ib will be to LEFT of _xa and _xb
  // ja and jb will be to BELOW _ya and _yb

  // bottom-left corner
  int ia = Utils::index_last_lt( _xa, *X, X->size() );
  int ja = Utils::index_last_lt( _ya, *Y, Y->size() );

  // top-right corner
  int ib = Utils::index_last_lt( _xb, *X, X->size() );
  int jb = Utils::index_last_lt( _yb, *Y, Y->size() );

   //We can integrate the function directly from the interpolation polynomial.
   //In Matrix form the polynomial looks like this
  
   ///               \ /             \ /       \
   //| x_2-x   x-x_1 | | F_11   F_12 | | y_2-x |
   //\               / |             | |       |
                     //| F_21   F_22 | | y-y_1 |
                     //\             / \       /
  
   //We can integrate this matrix equation to get the
   //matrix equation for the integral.
  
   ///                                   \ /             \ /                 \
   //| x_2 x - 0.5 x^2   0.5 x^2 - x_1 x | | F_11   F_12 | | y_2 y - 0.5 y^2 |
   //\                                   / |             | |                 |
                                         //| F_21   F_22 | | 0.5 y^2 - y_1 y |
                                         //\             / \                 /
   
   //This must be evaluated at the limits for x and y. Note that if the limits on x are [a,b],
   //and the limits on y are [c,d], then to evaluate a function f(x,y) at the limits
  
          //|xb |xb
   //f(x,y) |   |   =  ( f(xb,yb) - f(xb,ya) ) - ( f(xa,yb) - f(xa,ya) )
          //|xa |xa

   // integrate the whole elements first...
   
   RowVector2 vxa,vxb;
   ColVector2 vya,vyb;
   Real xa,xb,ya,yb;
   Real sum = 0;
   int i,j;
   for(i = ia+1; i < ib; i++)
   {
     xa = (*X)(i);
     xb = (*X)(i+1);
     vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
     vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;
     for(j = ja+1; j < jb; j++)
     {
       ya = (*Y)(j);
       yb = (*Y)(j+1);
       vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
       vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

       sum += vxb*Q(i,j)*vyb;
       sum -= vxb*Q(i,j)*vya;
       sum -= vxa*Q(i,j)*vyb;
       sum += vxa*Q(i,j)*vya;
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
   xb = (*X)(i+1);
   vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
   vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;
   for(j = ja+1; j < jb; j++)
   {
     ya = (*Y)(j);
     yb = (*Y)(j+1);
     vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
     vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

     sum += vxb*Q(i,j)*vyb;
     sum -= vxb*Q(i,j)*vya;
     sum -= vxa*Q(i,j)*vyb;
     sum += vxa*Q(i,j)*vya;
   }

   // add contributions from right side (excluding corners)
   i = ib;
   xa = (*X)(i);
   xb = _xb;
   vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
   vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;
   for(j = ja+1; j < jb; j++)
   {
     ya = (*Y)(j);
     yb = (*Y)(j+1);
     vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
     vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

     sum += vxb*Q(i,j)*vyb;
     sum -= vxb*Q(i,j)*vya;
     sum -= vxa*Q(i,j)*vyb;
     sum += vxa*Q(i,j)*vya;
   }

   if( ia == ib )
   {
     // the min and max limits on x are in the same element.
     // we need to subtract the integral over the entire element (see above)
     i  = ia;
     xa = (*X)(i);
     xb = (*X)(i+1);
     vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
     vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;
     for(j = ja+1; j < jb; j++)
     {
       ya = (*Y)(j);
       yb = (*Y)(j+1);
       vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
       vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

       sum += -vxb*Q(i,j)*vyb;
       sum -= -vxb*Q(i,j)*vya;
       sum -= -vxa*Q(i,j)*vyb;
       sum += -vxa*Q(i,j)*vya;
     }
   }


   // add contributions from bottom side (excluding corners)
   j = ja;
   ya = _ya;
   yb = (*Y)(j+1);
   vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
   vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

   for(i = ia+1; i < ib; i++)
   {
     xa = (*X)(i);
     xb = (*X)(i+1);
     vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
     vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;
     sum += vxb*Q(i,j)*vyb;
     sum -= vxb*Q(i,j)*vya;
     sum -= vxa*Q(i,j)*vyb;
     sum += vxa*Q(i,j)*vya;
   }

   // add contributions from top side (excluding corners)
   j = jb;
   ya = (*Y)(j);
   yb = _yb;
   vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
   vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

   for(i = ia+1; i < ib; i++)
   {
     xa = (*X)(i);
     xb = (*X)(i+1);
     vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
     vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;
     sum += vxb*Q(i,j)*vyb;
     sum -= vxb*Q(i,j)*vya;
     sum -= vxa*Q(i,j)*vyb;
     sum += vxa*Q(i,j)*vya;
   }

   if( ja == jb )
   {
     // the min and max limits on y are in the same element.
     // we need to subtract the integral over the entire element (see above)
     j = ja;
     ya = (*Y)(j);
     yb = (*Y)(j+1);
     vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
     vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

     for(i = ia+1; i < ib; i++)
     {
       xa = (*X)(i);
       xb = (*X)(i+1);
       vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
       vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;
       sum += -vxb*Q(i,j)*vyb;
       sum -= -vxb*Q(i,j)*vya;
       sum -= -vxa*Q(i,j)*vyb;
       sum += -vxa*Q(i,j)*vya;
     }
   }





   if( ia != ib && ja != jb )
   {

   // bottom-left corner
   i = ia;
   xa = _xa;
   xb = (*X)(i+1);
   vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
   vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;

   j = ja;
   ya = _ya;
   yb = (*Y)(j+1);
   vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
   vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

   sum += vxb*Q(i,j)*vyb;
   sum -= vxb*Q(i,j)*vya;
   sum -= vxa*Q(i,j)*vyb;
   sum += vxa*Q(i,j)*vya;
   

   // bottom-right corner
   i = ib;
   xa = (*X)(i);
   xb = _xb;
   vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
   vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;

   j = ja;
   ya = _ya;
   yb = (*Y)(j+1);
   vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
   vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

   sum += vxb*Q(i,j)*vyb;
   sum -= vxb*Q(i,j)*vya;
   sum -= vxa*Q(i,j)*vyb;
   sum += vxa*Q(i,j)*vya;


   // top-left corner
   i = ia;
   xa = _xa;
   xb = (*X)(i+1);
   vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
   vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;

   j = jb;
   ya = (*Y)(j);
   yb = _yb;
   vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
   vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

   sum += vxb*Q(i,j)*vyb;
   sum -= vxb*Q(i,j)*vya;
   sum -= vxa*Q(i,j)*vyb;
   sum += vxa*Q(i,j)*vya;


   // top-right corner
   i = ib;
   xa = (*X)(i);
   xb = _xb;
   vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
   vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;

   j = jb;
   ya = (*Y)(j);
   yb = _yb;
   vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
   vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

   sum += vxb*Q(i,j)*vyb;
   sum -= vxb*Q(i,j)*vya;
   sum -= vxa*Q(i,j)*vyb;
   sum += vxa*Q(i,j)*vya;

   }

   if( ia == ib && ja != jb )
   {
     xa = _xa;
     xb = _xb;
     vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
     vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;

     // bottom block
     j = ja;
     ya = _ya;
     yb = (*Y)(j+1);
     vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
     vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

     sum += vxb*Q(i,j)*vyb;
     sum -= vxb*Q(i,j)*vya;
     sum -= vxa*Q(i,j)*vyb;
     sum += vxa*Q(i,j)*vya;

     // top block
     j = jb;
     ya = (*Y)(j);
     yb = _yb;
     vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
     vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

     sum += vxb*Q(i,j)*vyb;
     sum -= vxb*Q(i,j)*vya;
     sum -= vxa*Q(i,j)*vyb;
     sum += vxa*Q(i,j)*vya;
   }

   if( ia != ib && ja == jb )
   {
     j = ja;
     ya = _ya;
     yb = _yb;
     vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
     vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

     // left block
     i = ia;
     xa = _xa;
     xb = (*X)(i+1);
     vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
     vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;

     sum += vxb*Q(i,j)*vyb;
     sum -= vxb*Q(i,j)*vya;
     sum -= vxa*Q(i,j)*vyb;
     sum += vxa*Q(i,j)*vya;

     // right block
     i = ib;
     xa = (*X)(i);
     xb = _xb;
     vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
     vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;

     sum += vxb*Q(i,j)*vyb;
     sum -= vxb*Q(i,j)*vya;
     sum -= vxa*Q(i,j)*vyb;
     sum += vxa*Q(i,j)*vya;

   }

   if( ia == ib && ja == jb )
   {

     i = ia;
     xa = _xa;
     xb = _xb;
     vxa << (*X)(i+1)*xa - 0.5*xa*xa, 0.5*xa*xa - (*X)(i)*xa;
     vxb << (*X)(i+1)*xb - 0.5*xb*xb, 0.5*xb*xb - (*X)(i)*xb;

     j = ja;
     ya = _ya;
     yb = _yb;
     vya << (*Y)(j+1)*ya - 0.5*ya*ya, 0.5*ya*ya - (*Y)(j)*ya;
     vyb << (*Y)(j+1)*yb - 0.5*yb*yb, 0.5*yb*yb - (*Y)(j)*yb;

     sum += vxb*Q(i,j)*vyb;
     sum -= vxb*Q(i,j)*vya;
     sum -= vxa*Q(i,j)*vyb;
     sum += vxa*Q(i,j)*vya;

     
   }






   return sum;
}



}

#endif // include protector
