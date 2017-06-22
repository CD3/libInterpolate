#include <celero/Celero.h>

#include<algorithm>
#include<cassert>

CELERO_MAIN

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

template<class Val, class Indexable>
int index_first_gt( Val val, const Indexable& vals, size_t N, int i = 0, int stride = 1 )
{
  if(i < 0) // don't let the user do harm...
    i = 0;
  // to find first element that is greater
  // we have to keep looking until element is not less
  while( i < N && vals[i] <= val )
  {
    i += stride;
  }

  return i/stride;
}


class SerialDataFixture : public celero::TestFixture
{
  public:
    SerialDataFixture()
    {
    }

    virtual std::vector<std::pair<int64_t, uint64_t>> getExperimentValues() const override
    {
      std::vector<std::pair<int64_t, uint64_t>> problemSpace;
      const int runs = 10;
      for( int i = 0; i < runs; i++)
      {
        problemSpace.push_back(std::make_pair(int64_t(pow(2,i+1)),uint64_t(0)));
      }

      return problemSpace;
    }
    virtual void setUp(int64_t N)
    {
      this->size = N;
      this->data_ptr = new double[this->size];
      for(int64_t i = 0; i < this->size; i++)
        this->data_ptr[i] = i;
    }
    virtual void tearDown()
    {
      delete this->data_ptr;
    }

    double *data_ptr;
    int64_t size;


};

BASELINE_F(SerialData_Begin, upper_bound_raw_pointer, SerialDataFixture, 100, 10000)
{
  int i = std::upper_bound( this->data_ptr, this->data_ptr+this->size, 0 ) - this->data_ptr;
  assert(i == 1);
}

BENCHMARK_F(SerialData_Begin, manual_raw_pointer, SerialDataFixture, 100, 10000 )
{
  int i = index_first_gt( 0, this->data_ptr, this->size );
  assert(i == 1);
}


BASELINE_F(SerialData_Middle, upper_bound_raw_pointer, SerialDataFixture, 100, 10000)
{
  int i = std::upper_bound( this->data_ptr, this->data_ptr+this->size, this->size/2 ) - this->data_ptr;
  assert(i == this->size/2+1);
}

BENCHMARK_F(SerialData_Middle, manual_raw_pointer, SerialDataFixture, 100, 10000 )
{
  int i = index_first_gt( this->size/2, this->data_ptr, this->size );
  assert(i == this->size/2+1);
}


BASELINE_F(SerialData_End, upper_bound_raw_pointer, SerialDataFixture, 100, 10000)
{
  int i = std::upper_bound( this->data_ptr, this->data_ptr+this->size, this->size-2 ) - this->data_ptr;
  assert(i == this->size-1);
}

BENCHMARK_F(SerialData_End, manual_raw_pointer, SerialDataFixture, 100, 10000 )
{
  int i = index_first_gt( this->size-2, this->data_ptr, this->size );
  assert(i == this->size-1);
}

