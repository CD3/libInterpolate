#include "catch.hpp"
#include "fakeit.hpp"

#include "Utils/Indexing.hpp"


TEST_CASE( "Indexing Utilities", "[utils]" ) {

  std::vector<int> x;
  x.push_back(1);
  x.push_back(2);
  x.push_back(3);

  REQUIRE( Utils::index_first_gt( 0  , x ) == 0 );
  REQUIRE( Utils::index_first_gt( 0.5, x ) == 0 );
  REQUIRE( Utils::index_first_gt( 1  , x ) == 1 );
  REQUIRE( Utils::index_first_gt( 1.5, x ) == 1 );
  REQUIRE( Utils::index_first_gt( 2  , x ) == 2 );
  REQUIRE( Utils::index_first_gt( 2.5, x ) == 2 );
  REQUIRE( Utils::index_first_gt( 3  , x ) == 3 );
  REQUIRE( Utils::index_first_gt( 3.5, x ) == 3 );

  REQUIRE( Utils::index_first_gt( 0  , x, 1 ) == 1 );
  REQUIRE( Utils::index_first_gt( 0.5, x, 1 ) == 1 );
  REQUIRE( Utils::index_first_gt( 1  , x, 1 ) == 1 );
  REQUIRE( Utils::index_first_gt( 1.5, x, 1 ) == 1 );
  REQUIRE( Utils::index_first_gt( 2  , x, 1 ) == 2 );
  REQUIRE( Utils::index_first_gt( 2.5, x, 1 ) == 2 );
  REQUIRE( Utils::index_first_gt( 3  , x, 1 ) == 3 );
  REQUIRE( Utils::index_first_gt( 3.5, x, 1 ) == 3 );



  REQUIRE( Utils::index_last_lt( 0  , x ) == -1 );
  REQUIRE( Utils::index_last_lt( 0.5, x ) == -1 );
  REQUIRE( Utils::index_last_lt( 1  , x ) == -1 );
  REQUIRE( Utils::index_last_lt( 1.5, x ) ==  0 );
  REQUIRE( Utils::index_last_lt( 2  , x ) ==  0 );
  REQUIRE( Utils::index_last_lt( 2.5, x ) ==  1 );
  REQUIRE( Utils::index_last_lt( 3  , x ) ==  1 );
  REQUIRE( Utils::index_last_lt( 3.5, x ) ==  2 );

  REQUIRE( Utils::index_last_lt( 0  , x, 1 ) == 1 );
  REQUIRE( Utils::index_last_lt( 0.5, x, 1 ) == 1 );
  REQUIRE( Utils::index_last_lt( 1  , x, 1 ) == 1 );
  REQUIRE( Utils::index_last_lt( 1.5, x, 1 ) == 1 );
  REQUIRE( Utils::index_last_lt( 2  , x, 1 ) == 1 );
  REQUIRE( Utils::index_last_lt( 2.5, x, 1 ) == 1 );
  REQUIRE( Utils::index_last_lt( 3  , x, 1 ) == 1 );
  REQUIRE( Utils::index_last_lt( 3.5, x, 1 ) == 2 );



  REQUIRE( Utils::index_first_ge( 0  , x ) == 0 );
  REQUIRE( Utils::index_first_ge( 0.5, x ) == 0 );
  REQUIRE( Utils::index_first_ge( 1  , x ) == 0 );
  REQUIRE( Utils::index_first_ge( 1.5, x ) == 1 );
  REQUIRE( Utils::index_first_ge( 2  , x ) == 1 );
  REQUIRE( Utils::index_first_ge( 2.5, x ) == 2 );
  REQUIRE( Utils::index_first_ge( 3  , x ) == 2 );
  REQUIRE( Utils::index_first_ge( 3.5, x ) == 3 );

  REQUIRE( Utils::index_first_ge( 0  , x, 1 ) == 1 );
  REQUIRE( Utils::index_first_ge( 0.5, x, 1 ) == 1 );
  REQUIRE( Utils::index_first_ge( 1  , x, 1 ) == 1 );
  REQUIRE( Utils::index_first_ge( 1.5, x, 1 ) == 1 );
  REQUIRE( Utils::index_first_ge( 2  , x, 1 ) == 1 );
  REQUIRE( Utils::index_first_ge( 2.5, x, 1 ) == 2 );
  REQUIRE( Utils::index_first_ge( 3  , x, 1 ) == 2 );
  REQUIRE( Utils::index_first_ge( 3.5, x, 1 ) == 3 );



  REQUIRE( Utils::index_last_le( 0  , x ) == -1 );
  REQUIRE( Utils::index_last_le( 0.5, x ) == -1 );
  REQUIRE( Utils::index_last_le( 1  , x ) ==  0 );
  REQUIRE( Utils::index_last_le( 1.5, x ) ==  0 );
  REQUIRE( Utils::index_last_le( 2  , x ) ==  1 );
  REQUIRE( Utils::index_last_le( 2.5, x ) ==  1 );
  REQUIRE( Utils::index_last_le( 3  , x ) ==  2 );
  REQUIRE( Utils::index_last_le( 3.5, x ) ==  2 );

  REQUIRE( Utils::index_last_le( 0  , x, 1 ) == 1 );
  REQUIRE( Utils::index_last_le( 0.5, x, 1 ) == 1 );
  REQUIRE( Utils::index_last_le( 1  , x, 1 ) == 1 );
  REQUIRE( Utils::index_last_le( 1.5, x, 1 ) == 1 );
  REQUIRE( Utils::index_last_le( 2  , x, 1 ) == 1 );
  REQUIRE( Utils::index_last_le( 2.5, x, 1 ) == 1 );
  REQUIRE( Utils::index_last_le( 3  , x, 1 ) == 2 );
  REQUIRE( Utils::index_last_le( 3.5, x, 1 ) == 2 );

}


