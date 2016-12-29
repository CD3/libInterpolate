#ifndef libInterp_Utils_Indexing_hpp
#define libInterp_Utils_Indexing_hpp

/** @file Utils.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 12/29/16
  */

namespace Utils {

/**
 * @brief find index of element in an array who's value is just greater than
 *
 * This function searches for the position of a value in an array, returing
 * the index of the array element that is just to the "right" (has a greater value)
 * of the value.
 *
 * It assumes that the array is sorted.
 */
template<class Val, class Indexable>
int index_just_greater( Val val, const Indexable& vals, size_t N, int i = 0 )
{
  while( i < N-1 && vals[i] < val )
      i++;

  return i;
}

template<class Val, class Indexable>
int index_just_greater( Val val, const Indexable& vals, int i = 0 )
{
  return index_just_greater( val, vals, vals.size(), i );
}

template<class Val, class Indexable>
int index_just_less( Val val, const Indexable& vals, size_t N, int i = 0 )
{
  i = index_just_greater( val, vals, N, i );
  i--;
  if( i < 0 )
    i = 0;

  return i;
}

template<class Val, class Indexable>
int index_just_less( Val val, const Indexable& vals, int i = 0 )
{
  return index_just_less( val, vals, vals.size(), i );
}

}

#endif // include protector
