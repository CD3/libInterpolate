#ifndef libInterp_Utils_Indexing_hpp
#define libInterp_Utils_Indexing_hpp

/** @file Utils.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 12/29/16
  *
  * @todo index search functions could be improved. perhaps using the mid-point method.
  */

namespace Utils {

/**
 * @brief find index of the first value in a sorted array that is greater than a given value.
 *
 * @param val the value to search for.
 * @param vals the array to be searched.
 * @param N the size of the array.
 * @param i the initial index to start searching from.
 */
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

  return i;
}

template<class Val, class Indexable>
int index_first_gt( Val val, const Indexable& vals, int i = 0, int stride = 1 )
{
  return index_first_gt( val, vals, vals.size(), i, stride );
}

/**
 * @brief Find index of last element in a sorted array with value less than a given value.
 *
 *
 * @param val the value to search for.
 * @param vals the array to be searched.
 * @param N the size of the array.
 * @param i the initial index to start searching from.
 */
template<class Val, class Indexable>
int index_last_lt( Val val, const Indexable& vals, size_t N, int i = -1, int stride = 1 )
{
  // optimization: if val is larger than largest value, just return N-1
  if( vals[N-1] < val )
    return N-1;

  // N is unsigned, so -1 < N will always be false
  // to find the last element that is less than val,
  // we keep looking until the next element is not less than val.
  while( (i < 0 || i < N-1) && vals[i+1] < val )
  {
    i += stride;
  }

  return i;
}

template<class Val, class Indexable>
int index_last_lt( Val val, const Indexable& vals, int i = -1, int stride = 1 )
{
  return index_last_lt( val, vals, vals.size(), i, stride );
}


template<class Val, class Indexable>
int index_first_ge( Val val, const Indexable& vals, size_t N, int i = 0, int stride = 1 )
{
  if(i < 0)
    i = 0;
  while( i < N && vals[i] < val )
  {
    i += stride;
  }

  return i;
}

template<class Val, class Indexable>
int index_first_ge( Val val, const Indexable& vals, int i = 0, int stride = 1 )
{
  return index_first_ge( val, vals, vals.size(), i, stride );
}

template<class Val, class Indexable>
int index_last_le( Val val, const Indexable& vals, size_t N, int i = -1, int stride = 1 )
{
  // optimization: if val is larger than largest value, just return N-1
  if( vals[N-1] < val )
    return N-1;

  // N is unsigned, so -1 < N will always be false
  while( (i < 0 || i < N-1) && vals[i+1] <= val )
  {
    i += stride;
  }

  return i;
}

template<class Val, class Indexable>
int index_last_le( Val val, const Indexable& vals, int i = -1, int stride = 1 )
{
  return index_last_le( val, vals, vals.size(), i, stride );
}

}

#endif // include protector
