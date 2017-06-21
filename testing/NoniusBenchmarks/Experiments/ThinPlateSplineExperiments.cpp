#define NONIUS_RUNNER
#include "nonius/nonius.h++"
#include "nonius/main.h++"

#include <eigen3/Eigen/Dense>

typedef double Real;

Real
G(Real x1, Real y1, Real x2, Real y2)
{
  if( x1 == x2 && y1 == y2 )
    return 0;

  Real r = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );

  return r*r*log(r);
}

NONIUS_BENCHMARK("EXPERIMENT | TPS Construction | Baseline | No tmp vars",
[](nonius::chronometer meter)
{

  typedef Eigen::Matrix<Real,Eigen::Dynamic,1> VectorType;
  typedef Eigen::Map<VectorType> VectorMapType;
  typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> MatrixType;

  int Nx = 5;
  int Ny = 5;
  VectorType x(Nx*Ny),y(Nx*Ny),z(Nx*Ny);

  for( int i = 0; i < Nx; i++ )
  {
    for( int j = 0; j < Ny; j++ )
    {
      x(i*Ny + j) = i;
      y(i*Ny + j) = j;
      z(i*Ny + j) = sin(i)*sin(j);
    }
  }

  MatrixType a = MatrixType(x.rows(),1);
  MatrixType b = MatrixType(x.rows(),3);

  MatrixType M = MatrixType(x.rows(),x.rows());
  MatrixType N = MatrixType(x.rows(),3);

  // rows
  for( int i = 0; i < x.rows(); i++ )
  {
    // N
    N(i,0) = 1;
    N(i,1) = x(i);
    N(i,2) = y(i);

    for( int j = 0; j < x.rows(); j++ )
      M(i,j) = G(x(i), y(i), x(j), y(j) );


  }

  meter.measure([&](){
  b = (N.transpose()*M.inverse()*N).inverse()*N.transpose()*M.inverse()*z;
  a = M.inverse()*(z - N*b);
  });

})


NONIUS_BENCHMARK("EXPERIMENT | TPS Construction | Test 1 | Tmp vars",
[](nonius::chronometer meter)
{

  typedef Eigen::Matrix<Real,Eigen::Dynamic,1> VectorType;
  typedef Eigen::Map<VectorType> VectorMapType;
  typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> MatrixType;

  int Nx = 5;
  int Ny = 5;
  VectorType x(Nx*Ny),y(Nx*Ny),z(Nx*Ny);

  for( int i = 0; i < Nx; i++ )
  {
    for( int j = 0; j < Ny; j++ )
    {
      x(i*Ny + j) = i;
      y(i*Ny + j) = j;
      z(i*Ny + j) = sin(i)*sin(j);
    }
  }

  MatrixType a = MatrixType(x.rows(),1);
  MatrixType b = MatrixType(x.rows(),3);

  MatrixType M = MatrixType(x.rows(),x.rows());
  MatrixType N = MatrixType(x.rows(),3);

  // rows
  for( int i = 0; i < x.rows(); i++ )
  {
    // N
    N(i,0) = 1;
    N(i,1) = x(i);
    N(i,2) = y(i);

    for( int j = 0; j < x.rows(); j++ )
      M(i,j) = G(x(i), y(i), x(j), y(j) );


  }

  MatrixType Minv = M.inverse();
  MatrixType Ntrans = N.transpose();


  meter.measure([&](){
  b = (Ntrans*Minv*N).inverse()*Ntrans*Minv*z;
  a = Minv*(z - N*b);
  });

})









