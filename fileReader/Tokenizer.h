#ifndef TOKENIZER_H
#define TOKENIZER_H
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <string>

/** \file
*/

namespace RUC
{
  /**
   * tokenizes a string based on a user defined set of delimiters and returns the tokens back through a vector reference
   * \param _str [in]                 the string to tokenize.
   * \param _tokens [out]             the vector of tokens returned. tokens are appended to the vector.
   * \param _delimiters [in]          list of characters that will delimit the tokens
   * \param _allow_empty_tokens [in]  flag that determines if empty tokens will be skipped or returned as empy strings
   */
  inline  void Tokenize(const std::string& _str, std::vector<std::string>& _tokens, const std::string& _delimiters = " ", bool _allow_empty_tokens = false)
  {
      std::string::size_type pos, lastPos;

      /* skips the possible delimiters at beginning */
      lastPos = 0;
      /* find first "delimiter" */
      pos     = _str.find_first_of(_delimiters, lastPos);

      while(pos != std::string::npos)
      {
        if( pos != lastPos )
        { /* found a token, add it to the vector */
          _tokens.push_back(_str.substr(lastPos, pos - lastPos));
        }
        else
        { /* we have found an empty token. only add it if this is allowed */
          if(_allow_empty_tokens)
            _tokens.push_back("");
        }

        /* find next delimiter */
        lastPos = pos+1;
        /* check that we aren't already at the end of the string
         * (this would happen if the last char is a delimter */
        if(lastPos == _str.size())
        { // last character is a delimiter
          if(_allow_empty_tokens)
            _tokens.push_back("");
          pos = lastPos;
          break;
        }

        pos = _str.find_first_of(_delimiters, lastPos);
      }

      if( pos != lastPos )
      {
        _tokens.push_back(_str.substr(lastPos, pos - lastPos));
      }
  }

  /**
   * tokenizes a string and returns the tokens in a vector
   * \param _str [in]                 the string to tokenize.
   * \param _delimiters [in]          list of characters that will delimit the tokens
   * \param _allow_empty_tokens [in]  flag that determines if empty tokens will be skipped or returned as empy strings
   * \return                  a vector of tokens.
   */
  inline std::vector<std::string> Tokenize(const std::string& _str, const std::string& _delimiters = " ", bool _allow_empty_tokens = false)
  {

    std::vector<std::string> tokens;
    Tokenize(_str, tokens, _delimiters, _allow_empty_tokens);
    return tokens;

  }
}


#endif
