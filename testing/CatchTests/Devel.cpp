
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



template<typename T>
struct can_call_Method {
  private:
   typedef std::true_type  yes;
   typedef std::false_type no;

   template<typename U>
   static auto test(int) -> decltype(std::declval<U>().Method(), yes());
   template<typename>
   static no test(...);

  public:
   static constexpr bool value = std::is_same<decltype(test<T>(0)), yes>::value;

};

template<typename T>
struct defines_Method {
  template<typename U, void (U::*)()>
  struct SFINAE {
  };
  template<typename U>
  static char Test(SFINAE<U, &U::Method>*);
  template<typename U>
  static int        Test(...);
  static const bool value = sizeof(Test<T>(0)) == sizeof(char);
};

struct ClassA{};
struct ClassB : public ClassA{void Method();};
struct ClassC : public ClassB{};


TEST_CASE("Method detection")
{
  CHECK(!can_call_Method<ClassA>::value );
  CHECK( can_call_Method<ClassB>::value );
  CHECK( can_call_Method<ClassC>::value );

  CHECK(!defines_Method<ClassA>::value );
  CHECK( defines_Method<ClassB>::value );
  CHECK(!defines_Method<ClassC>::value );
}





