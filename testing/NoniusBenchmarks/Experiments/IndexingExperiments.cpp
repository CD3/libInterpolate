#define NONIUS_RUNNER
#include "nonius/nonius.h++"
#include "nonius/main.h++"

#include <algorithm>
#include <eigen3/Eigen/Dense>

namespace IndexingExperiments {

struct Fixture
{
  typedef Eigen::Matrix<double,Eigen::Dynamic,1> VectorType;
  VectorType v;

  Fixture(size_t N): v(N)
  {
    for(int i = 0; i < N; i++)
      v(i) = 2.*i;
  }

};


template<class Val, class Indexable>
int index_first_gt( Val val, const Indexable& vals, size_t N, int i = 0, int stride = 1 )
{
  if(i < 0)
    i = 0;
  while( i < N && vals[i] <= val )
  {
    i += stride;
  }

  return i;
}


NONIUS_BENCHMARK("Search for upper bound | Baseline | stl on pointer",
[](nonius::chronometer meter)
{
  Fixture f(500);

  meter.measure( [&](){

  for(int j = 0; j < 500; j++)
  {
    double *p = std::upper_bound( f.v.data(), f.v.data()+500, 2.5*j );
    int i = p - f.v.data();
  }

  });


})

NONIUS_BENCHMARK("Search for upper bound | Trial 1 | manual on eigen matrix",
[](nonius::chronometer meter)
{
  Fixture f(500);

  meter.measure( [&](){

  for(int j = 0; j < 500; j++)
  {
    int i = index_first_gt( 2.5*j, f.v, 500 );
  }

  });

})

NONIUS_BENCHMARK("Search for upper bound | Trial 1 | manual on pointer",
[](nonius::chronometer meter)
{
  Fixture f(500);

  meter.measure( [&](){

  for(int j = 0; j < 500; j++)
  {
    int i = index_first_gt( 2.5*j, f.v.data(), 500 );
  }

  });

})

}
