#define NONIUS_RUNNER
#include "nonius/nonius.h++"
#include "nonius/main.h++"

#include <algorithm>
#include <eigen3/Eigen/Dense>
#include "Utils/Indexing.hpp"


/* Tese experiments/benchmarks/tests are to help
 * decide if we should switch to the std::upper_bound algorithm,
 * or keep using our custom index finder.
 *
 * 2017/6/20 testing indicates that the std::upper_bound is *slower*
 */

namespace Utils {
  // an iterator that wraps a pointer with
  // a stride.
template <typename T>
struct StridedIterator : public std::iterator<std::forward_iterator_tag, T> {
    StridedIterator(T *t, unsigned stride) : ptr(t), stride(stride) {}
    bool operator==(const StridedIterator<T> &other) const
    {
        return ptr == other.ptr;
    }
    bool operator!=(const StridedIterator<T> &other) const
    {
        return ptr != other.ptr;
    }
    T *operator->() const { return ptr; }
    T &operator*() const { return *ptr; }
    StridedIterator &operator++()
    {
        ptr += stride;
        return *this;
    }
    StridedIterator operator++(int)
    {
        auto ret = *this;
        ++*this;
        return ret;
    }
    StridedIterator operator+(int amt) const
    {
        auto ret = StridedIterator(ptr + amt * stride, stride);
        return ret;
    }
    int operator-(const StridedIterator &right) const
    {
        auto ret = (ptr - right.ptr)/stride;
        return ret;
    }

  private:
    unsigned stride;
    T *ptr;
};
}

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



NONIUS_BENCHMARK("Search for upper bound : sequential | Baseline | stl on pointer",
[](nonius::chronometer meter)
{
  meter.measure( [&](){
  int N = 500*500;
  Fixture f(N);


  for(int j = 0; j < N; j++)
  {
    int i = std::upper_bound( f.v.data(), f.v.data()+N, 2*j ) - f.v.data();
  }

  });

})

NONIUS_BENCHMARK("Search for upper bound : sequential | Trial 1 | manual on eigen matrix",
[](nonius::chronometer meter)
{
  meter.measure( [&](){
  int N = 500*500;
  Fixture f(N);


  for(int j = 0; j < N; j++)
  {
    int i = Utils::index_first_gt( 2*j, f.v, N, 0, 1);
  }

  });

})




NONIUS_BENCHMARK("Search for upper bound : strided | Baseline | stl on sequential data",
[](nonius::chronometer meter)
{
  meter.measure( [&](){
  int Nx = 500, Ny = 500;
  Fixture f(Nx*Ny);


  for(int j = 0; j < Nx; j++)
  {
    int i = std::upper_bound( f.v.data(), f.v.data()+Nx, 2*j ) - f.v.data();
  }

  });
})

NONIUS_BENCHMARK("Search for upper bound : strided | Trail 1 | stl on strided iterator",
[](nonius::chronometer meter)
{
  meter.measure( [&](){
  int Nx = 500, Ny = 500;
  Fixture f(Nx*Ny);

  typedef Eigen::Matrix<double,Eigen::Dynamic,1> VectorType;
  typedef Eigen::Map<VectorType,Eigen::Unaligned,Eigen::InnerStride<Eigen::Dynamic>> _2DVectorView;
  _2DVectorView view( f.v.data(), Nx, Eigen::InnerStride<Eigen::Dynamic>(Ny) );

  Utils::StridedIterator<double> begin(f.v.data(), Ny);
  Utils::StridedIterator<double> end(f.v.data()+Nx*Ny, Ny);


  for(int j = 0; j < Nx; j++)
  {
    int i = std::upper_bound( begin, end, 2*Ny*j ) - begin; 
  }

  });
})

NONIUS_BENCHMARK("Search for upper bound : strided | Trail 2 | manual on srided eigen matrix",
[](nonius::chronometer meter)
{
  meter.measure( [&](){
  int Nx = 500, Ny = 500;
  Fixture f(Nx*Ny);

  typedef Eigen::Matrix<double,Eigen::Dynamic,1> VectorType;
  typedef Eigen::Map<VectorType,Eigen::Unaligned,Eigen::InnerStride<Eigen::Dynamic>> _2DVectorView;
  _2DVectorView view( f.v.data(), Nx, Eigen::InnerStride<Eigen::Dynamic>(Ny) );


  for(int j = 0; j < Nx; j++)
  {
    int i = Utils::index_first_gt( 2*Ny*j, f.v, Nx, 0, 1 );
  }

  });
})




}
