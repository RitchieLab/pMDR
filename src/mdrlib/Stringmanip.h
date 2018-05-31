//Stringmanip.h

#ifndef datahandling__STRINGMANIP_H__
#define datahandling__STRINGMANIP_H__

#include <sstream>
#include <vector>

using namespace std;

namespace mdr{
const int LINEMAX = 512;

///
/// Performs various string manipulation functions
///

/// String manipulator
class Stringmanip{
  
  public:
    Stringmanip(); 
    /// check if number contained in string
    static int is_number(string str);
    /// convert to uppercase
    static string to_upper(string convert);
    /// convert to lowercase
    static string to_lower(string convert);
    
    /// splits string 
    static int split_string(const string& input, 
      string delim1, string delim2, vector<string>& results);
    static int stoi(string number); 
    static unsigned stouint(string number);
    static int stodata(string number);
    static string itos(long number);
    static string itos(float number);
    static string itos(int number);
    static string itos(unsigned int number);
    static string itos(double number);
    static double stodouble(string number);
    static void trim(string& str);
 
  private:
    
  
};
}

#endif   

