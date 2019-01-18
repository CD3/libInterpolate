#include "catch.hpp"
#include "fakeit.hpp"

// This file contains examples that demonstrate how to write simple tests.

#include<vector>
#include<string>

#include<Interp.hpp>
#include<AnyInterpolator.hpp>

unsigned int Factorial( unsigned int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}

TEST_CASE( "Simple example from the Catch website", "[example]" ) {
    REQUIRE( Factorial(1) == 1 );
    REQUIRE( Factorial(2) == 2 );
    REQUIRE( Factorial(3) == 6 );
    REQUIRE( Factorial(10) == 3628800 );
}

TEST_CASE( "Example tests using Sections from Catch website", "[example]" ) {

    std::vector<int> v( 5 );

    REQUIRE( v.size() == 5 );
    REQUIRE( v.capacity() >= 5 );

    SECTION( "resizing bigger changes size and capacity" ) {
        v.resize( 10 );

        REQUIRE( v.size() == 10 );
        REQUIRE( v.capacity() >= 10 );
    }
    SECTION( "resizing smaller changes size but not capacity" ) {
        v.resize( 0 );

        REQUIRE( v.size() == 0 );
        REQUIRE( v.capacity() >= 5 );
    }
    SECTION( "reserving bigger changes capacity but not size" ) {
        v.reserve( 10 );

        REQUIRE( v.size() == 5 );
        REQUIRE( v.capacity() >= 10 );
    }
    SECTION( "reserving smaller does not change size or capacity" ) {
        v.reserve( 0 );

        REQUIRE( v.size() == 5 );
        REQUIRE( v.capacity() >= 5 );
    }
}

TEST_CASE("Float comparisons","[example]")
{

  REQUIRE( 1.0 == Approx(1.0) );
  REQUIRE( 1.0 != Approx(1.1) );
  REQUIRE( 1.0 == Approx(1.1).epsilon(0.1) ); // epsilon specifies absolute tolerance
  REQUIRE( 1.2345 != Approx(1) );
  // epsilon examples
  REQUIRE( 1.2345 != Approx(1.2344) );
  REQUIRE( 1.2345 != Approx(1.2344).epsilon(0.000001) );
  REQUIRE( 1.2345 != Approx(1.2344).epsilon(0.00001) );
  REQUIRE( 1.2345 == Approx(1.2344).epsilon(0.0001) );
  REQUIRE( 1.2345 == Approx(1.2344).epsilon(0.001) );
  REQUIRE( 1.2345 == Approx(1.2344).epsilon(0.01) );
  REQUIRE( 1.2345 == Approx(1.2344).epsilon(0.1) );

  REQUIRE( 1.2345 != Approx(1.233) );
  REQUIRE( 1.2345 != Approx(1.233).epsilon(0.000001) );
  REQUIRE( 1.2345 != Approx(1.233).epsilon(0.00001) );
  REQUIRE( 1.2345 != Approx(1.233).epsilon(0.0001) );
  REQUIRE( 1.2345 == Approx(1.233).epsilon(0.001) );
  REQUIRE( 1.2345 == Approx(1.233).epsilon(0.01) );
  REQUIRE( 1.2345 == Approx(1.233).epsilon(0.1) );

  REQUIRE( 1.2345 != Approx(1.22) );
  REQUIRE( 1.2345 != Approx(1.22).epsilon(0.000001) );
  REQUIRE( 1.2345 != Approx(1.22).epsilon(0.00001) );
  REQUIRE( 1.2345 != Approx(1.22).epsilon(0.0001) );
  REQUIRE( 1.2345 != Approx(1.22).epsilon(0.001) );
  REQUIRE( 1.2345 == Approx(1.22).epsilon(0.01) );
  REQUIRE( 1.2345 == Approx(1.22).epsilon(0.1) );

  REQUIRE( 1.2345 != Approx(1.1) );
  REQUIRE( 1.2345 != Approx(1.1).epsilon(0.000001) );
  REQUIRE( 1.2345 != Approx(1.1).epsilon(0.00001) );
  REQUIRE( 1.2345 != Approx(1.1).epsilon(0.0001) );
  REQUIRE( 1.2345 != Approx(1.1).epsilon(0.001) );
  REQUIRE( 1.2345 != Approx(1.1).epsilon(0.01) );
  REQUIRE( 1.2345 == Approx(1.1).epsilon(0.1) );

  REQUIRE( 1.2345 != Approx(1.0) );
  REQUIRE( 1.2345 != Approx(1.0).epsilon(0.000001) );
  REQUIRE( 1.2345 != Approx(1.0).epsilon(0.00001) );
  REQUIRE( 1.2345 != Approx(1.0).epsilon(0.0001) );
  REQUIRE( 1.2345 != Approx(1.0).epsilon(0.001) );
  REQUIRE( 1.2345 != Approx(1.0).epsilon(0.01) );
  REQUIRE( 1.2345 != Approx(1.0).epsilon(0.1) );
}


class A
{
  public:
    int func_a() { return 1; }
    int func_b() { return 2; }
    int func_c() { return 3; }
};

class B
{
  public:
    virtual int func_a() { return 1; }
    virtual int func_b() { return 2; }
            int func_c() { return 3; }
};

class C
{
  public:
    virtual int func_a() { return 1; }
    virtual int func_b() { return 2; }
    virtual int func_c() { return func_a(); }
};

// Mock objects are face objects that mimic your project's objects. They are useful
// if you want to test a function that takes an object as an argument and is going to
// call methods on that object that either take a long time to run, or return un-deterministic
// value. By creating a mock object that returns a known value quickly, you can
// write unit tests for these functions.
TEST_CASE("Mock object", "[example]")
{

  SECTION("Simple Mocking")
  {
    // fakeit::Mock<A> mA;  //Won't Compile. Can only mock virtual methods
    fakeit::Mock<B> mB;

    SECTION("Stub single method")
    {
      fakeit::When(Method(mB,func_a)).AlwaysReturn(10);
      B &b = mB.get();

      CHECK( b.func_a() == 10 );
      //CHECK( b.func_b() == 2 ); //Will Compile, but Fail. Only stubbed functions are defined.
    }

    SECTION("Stub multiple methods")
    {
      fakeit::When(Method(mB,func_a)).AlwaysReturn(10);
      fakeit::When(Method(mB,func_b)).AlwaysReturn(20);
      //fakeit::When(Method(mB,func_c)).AlwaysReturn(30); //Will Compile, but Fail. Can only stub virtual methods.
      B &b = mB.get();

      CHECK( b.func_a() == 10 );
      CHECK( b.func_b() == 20 );
    }
  }

  SECTION("Simple Spying")
  {
    B obj;
    fakeit::Mock<B> mB(obj);

    SECTION("Stub single method")
    {
      fakeit::When(Method(mB,func_a)).AlwaysReturn(10);
      B &b = mB.get();

      CHECK( b.func_a() == 10 );
      CHECK( b.func_b() == 2 ); // Works now, spying uses all methods
      CHECK( b.func_c() == 3 ); // Can even all non-virtual methods
    }

  }


  SECTION("Spying with inter-func calls")
  {
    C obj;
    fakeit::Mock<C> mC(obj);

    SECTION("Replace caller")
    {
      fakeit::When(Method(mC,func_c)).AlwaysReturn(10);
      C &c = mC.get();

      CHECK( c.func_a() == 1 );
      CHECK( c.func_b() == 2 );
      CHECK( c.func_c() == 10 );
    }

    SECTION("Replace callee")
    {
      fakeit::When(Method(mC,func_a)).AlwaysReturn(10);
      C &c = mC.get();

      CHECK( c.func_a() == 10 );
      CHECK( c.func_b() == 2 );
      CHECK( c.func_c() == 10 );
    }

  }
  
}

TEST_CASE("Runtime dispatch with AnyInterpolator")
{

  std::string method;
  std::vector<double> x,y;

  x.push_back(0);
  x.push_back(1);
  x.push_back(2);

  y.push_back(10);
  y.push_back(20);
  y.push_back(30);

  method = "linear";

  _1D::AnyInterpolator<double> interp;

  // select the interpolation method on user input
  if( method == "linear")
    interp = _1D::LinearInterpolator<double>();
  if( method == "cubicspline")
    interp = _1D::CubicSplineInterpolator<double>();
  if( method == "monotonic")
    interp = _1D::MonotonicInterpolator<double>();

  interp.setData( x.size(), x.data(), y.data() );

  // interpolation is done with the operator() method.
  double val = interp(1.5);

  CHECK(val == Approx(25));
}

