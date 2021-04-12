
#include "catch.hpp"
#include <vector>
#include <algorithm>
#include <numeric>

TEST_CASE("std::copy between different types")
{
  std::vector<double> x1(10);
  std::vector<float> x2(10);
  std::iota(x1.begin(),x1.end(),10);

  std::copy(x1.begin(), x1.end(), x2.begin());

  CHECK(x1[9] == Approx(19));

}

TEST_CASE("2D Point Sorting")
{

  std::vector<double> x,y,z;

  x.push_back(0); y.push_back(1); z.push_back(1);
  x.push_back(0); y.push_back(-1); z.push_back(1);
  x.push_back(1); y.push_back(0); z.push_back(1);
  x.push_back(-1); y.push_back(0); z.push_back(1);

  x.push_back(1); y.push_back(1); z.push_back(1);
  x.push_back(1); y.push_back(-1); z.push_back(1);
  x.push_back(1); y.push_back(1); z.push_back(1);
  x.push_back(-1); y.push_back(1); z.push_back(1);








  


}
