#ifndef READ_FUNCTION_H
#define READ_FUNCTION_H
/*
 * =====================================================================================
 *
 *       Filename:  ReadFunction.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/18/10 14:11:24
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  C.D. Clark III (), clifton.clark.ctr@brooks.af.mil
 *        Company:  TASC, Inc.
 *
 * =====================================================================================
 */
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <regex>

#include "Tokenizer.h"
#include "DataConverter.h"


namespace RUC
{
  /** Reads a function from an input stream.
   * Memory for _x and _y is allocated, and _n is set, based
   * on the number of data points found in the stream _in. _x and _y should NOT point
   * to allocated memory before being passed to ReadFunction() because the memory will
   * be lost (a memory leak). the multiplicity and dimensions can also be given, which specify
   * how many columns are expected to be coordinates (dimension) and how many are values (multiplicity). so, for example, a
   * complex function f(x), has _multiplicity = 2 and _dimension = 1.
   * \param _in           input stream to read data from. data is while state of _in is good
   * \param _x            pointer through which coordinates are returned
   * \param _y            pointer through which values are returned
   * \param _n            pointer through which number of points in each dimension is returned
   * \param _multiplicity specifies the number of columns that make up a single value of the function
   * \param _dimensions   specifies the number of columns that are coordinates
   * */
  template < typename ArgType, typename ValType >
  void ReadFunction(std::istream &_in, ArgType  *&_x, ValType *&_y, int *&_n, int _multiplicity = 1, int _dimensions = 1)
  {

    ArgType *xbuffer1, *xbuffer2;
    ValType *ybuffer1, *ybuffer2;

    //////////////////////////////
    //  READ IN DATA FROM FILE  //
    //////////////////////////////

    const int chunck = 1000;
    int buffsize;
    int i,j;
    int n;

    std::string line;

    buffsize = chunck;


    xbuffer1 = new ArgType[buffsize*_dimensions  ];
    ybuffer1 = new ValType[buffsize*_multiplicity];
    
    std::regex re("(#.*)");

    i = 0;
    while( getline( _in, line ) )
    {
        std::smatch match;
        std::string result;
        //if(!std::regex_search(line, match, re))
        if( !std::regex_match(line,  re) )
        {
          std::vector<std::string> columns = Tokenize(line);

          if(i == buffsize - 1)  // this "grows" the buffer
          {
            buffsize = buffsize + chunck;

            xbuffer2 = new ArgType[buffsize*_dimensions  ];
            ybuffer2 = new ValType[buffsize*_multiplicity];
            for(j = 0; j < buffsize*_dimensions; j++)
              xbuffer2[j] = xbuffer1[j];
            for(j = 0; j < buffsize*_multiplicity; j++)
              ybuffer2[j] = ybuffer1[j];

            delete[] xbuffer1;
            delete[] ybuffer1;

            xbuffer1 = xbuffer2;
            ybuffer1 = ybuffer2;
          }

          for(j = 0; j < _dimensions; j++)
            xbuffer1[i*_dimensions + j]   = fromString<ArgType>(columns[j]);
          for(j = 0; j < _multiplicity; j++)
            ybuffer1[i*_multiplicity + j] = fromString<ValType>(columns[j+_dimensions]);
          i++;
        }
    }

    n = i;


    /** \todo add actual multi-coordinate support */

    _n = new int[_dimensions];
    _x = new ArgType[n * _dimensions  ];
    _y = new ValType[n * _multiplicity];

    for(j = 0; j < n * _dimensions  ; j++)
      _x[j] = xbuffer1[j];
    for(j = 0; j < n * _multiplicity; j++)
      _y[j] = ybuffer1[j];
    for(j = 0; j < _dimensions; j++)
      _n[j] = n;

    delete[] xbuffer1;
    delete[] ybuffer1;

    }

  /** specialized version of ReadFunction that takes an integer for _n (rather than an integer array)
   * and reads in a 1D function
   */
  template < typename ArgType, typename ValType >
  void ReadFunction(std::istream &_in, ArgType  *&_x, ValType *&_y, int &_n, int _multiplicity = 1)
  {
    int *n;
    ReadFunction(_in, _x, _y, n, _multiplicity);
    _n = n[0];
    delete[] n;
  }
}


#endif
