//Stringmanip.cpp

#include "Stringmanip.h"
namespace mdr{
///
/// checks that the string passed in is a number
/// @param str string to check as number
/// @return true (1)  or false (0)
///
int Stringmanip::is_number(string str){
  for(unsigned int i=0; i<str.length(); i++){
    if((str[i]<'0'||str[i]>'9') && str[i] != '.'){
      return 0;
    }
  }
  return 1;
}


/// 
/// Converts string to uppercase
/// @param convert string to convert to uppercase
/// @return converted string
///
string Stringmanip::to_upper(string convert){
  // convert keyword to uppercase
  unsigned int strpos, length = convert.length();
  for(strpos=0; strpos<length; strpos++){
    convert[strpos] = toupper(convert[strpos]);
  }
  return convert;
}


/// 
/// Converts string to lowercase
/// @param convert string to convert to lowercase
/// @return converted string
///
string Stringmanip::to_lower(string convert){
  // convert keyword to uppercase
  unsigned int strpos, length = convert.length();
  for(strpos=0; strpos<length; strpos++){
    convert[strpos] = tolower(convert[strpos]);
  }
  return convert;
}


///
/// Split string into tokens as defined by the delimiter
/// <br> delim1 and delim2 are concatenated to
/// create the token to split on
/// @param input string to split
/// @param delim1 first delimiter to use
/// @param delim2 second delimiter to use
/// @param results [out] stores string tokens
/// @return number of tokens
///
int Stringmanip::split_string(const string &input, 
  string delim1, string delim2, vector<string> &results){
  int iPos = 0;
  int newPos = -1;
  string delimiter = delim1 + delim2;;

  int sizeS2 = 1;
  int isize = input.size();
  vector<unsigned int> positions;
  newPos = input.find_first_of (delimiter, 0);
  if( newPos < 0 ) { return 0; }

  int numFound = 0;

  while( newPos >= iPos ){
    numFound++;
    positions.push_back(newPos);
    iPos = newPos;
    newPos = input.find_first_of (delimiter, iPos+sizeS2);
  }

  for(unsigned int i=0; i <= positions.size(); i++ ){
    string s;
    if( i == 0 ) { s = input.substr( i, positions[i] ); }
    int offset = positions[i-1] + sizeS2;
    if( offset < isize )
    {
      if( i == positions.size() )
      {
        s = input.substr(offset);
      }
      else if( i > 0 )
      {
        s = input.substr( positions[i-1] + sizeS2, 
          positions[i] - positions[i-1] - sizeS2 );
      }
    }
    if( s.size() > 0 )
    {
      results.push_back(s);
    }
  }

  return results.size();
}

///
/// Converts a string into the int representation
/// @param number to convert
/// @return number 
///
int Stringmanip::stoi(string number){
    int value;
    stringstream ss(number);
    ss >> value;
    return value;
}

///
/// Converts a string into the unsigned int representation
/// @param number to convert
/// @return number 
///
unsigned Stringmanip::stouint(string number){
  unsigned int value;
  stringstream ss(number);
  ss >> value;
  return value;
}

///
/// Converts a string into the int representation
/// @param number to convert
/// @return number 
///
int Stringmanip::stodata(string number){
  int value;
  stringstream ss(number);
  ss >> value;
  return value;
}

///
/// Converts a string into the double representation
/// @param number to convert
/// @return number 
///
double Stringmanip::stodouble(string number){
  double value;
  stringstream ss(number);
  ss >> value;
  return value;
}

///
/// Converts an unsigned int into the string representation
/// @param number to convert
/// @return number as a string
///
string Stringmanip::itos(unsigned int number){
  stringstream oss;
  oss << number;
  return oss.str();
}

///
/// Converts a long into the string representation
/// @param number to convert
/// @return number as a string
///
string Stringmanip::itos(long number){ stringstream oss;
  oss << number;
  return oss.str();
}

///
/// Converts an int into the string representation
/// @param number to convert
/// @return number as a string
///
string Stringmanip::itos(int number){
  stringstream oss;
  oss << number;
  return oss.str();
}

///
/// Converts a float into the string representation
/// @param number to convert
/// @return number as a string
///
string Stringmanip::itos(float number){
  stringstream oss;
  oss << number;
  return oss.str();
}

///
/// Converts a double into the string representation
/// @param number to convert
/// @return number as a string
///
string Stringmanip::itos(double number){
  stringstream oss;
  oss << number;
  return oss.str();
}

///
/// Trims spaces from string
///
void Stringmanip::trim(string& str)
{
  string::size_type pos = str.find_last_not_of(' ');
  if(pos != string::npos) {
    str.erase(pos + 1);
    pos = str.find_first_not_of(' ');
    if(pos != string::npos) str.erase(0, pos);
  }
  else str.erase(str.begin(), str.end());
}

}
