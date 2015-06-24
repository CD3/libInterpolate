#ifndef DATA_CONVERTER_H
#define DATA_CONVERTER_H
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <map>
#include <string>
#include <sstream>
#include <stdexcept>

/** \file
 * a set of functions to convert numbers to strings and strings to numbers
*/


namespace RUC
{
  template<typename T>
  std::string toString(T data)
  {
    std::stringstream tmpStream;
    tmpStream << data;

    return tmpStream.str();
  }

  template<typename T>
  std::string toString()
  {
  }



  /** some functions to return the string of types */
  template<>
  inline std::string toString<int>()
  {
    return "int";
  }

  template<>
  inline std::string toString<float>()
  {
    return "float";
  }

  template<>
  inline std::string toString<double>()
  {
    return "double";
  }

  template<>
  inline std::string toString<long int>()
  {
    return "long int";
  }

  template<>
  inline std::string toString<long double>()
  {
    return "long double";
  }



  template<typename T>
  T fromString(std::string data)
  {
    std::stringstream tmpStream(data);
    T tmpData;
    tmpStream >> tmpData;

    return tmpData;
  }

  template<typename T>
  T fromMap( std::map<std::string, std::string>& _vars, std::string _key, std::string _calling_function)
  {
    if( _vars.find(_key) == _vars.end() )
      throw std::runtime_error("Key Not Found: "+_calling_function+" expected map to contain \""+_key+"\" key");

    return fromString<T>( _vars.find(_key)->second );
  }

}


#endif
