//KnuthKnuthComboGenerator.cpp

//     Copyright 2008 Marylyn Ritchie

//     This file is part of MDR.
// 
//     MDR is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
// 
//     MDR is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License
//     along with MDR.  If not, see <http://www.gnu.org/licenses/>.


#include<iostream>
using namespace std;
#include "KnuthComboGenerator.h"
#include <fstream>

namespace mdr{

///
/// Constructor initializes parameters
///
KnuthComboGenerator::KnuthComboGenerator(){
  initialize();
}

///
/// Initializes variables for generator
/// @return none
///
void KnuthComboGenerator::initialize(){
  c = NULL;
}


///
/// Creates initial state 
///
void KnuthComboGenerator::initialize_state(){
  kdec = ComboEnd;
  x = 0;
  j = kdec;

  c = new int[ComboEnd+3];

  for(int i=1; i <= kdec; i++){
    c[i] = i;
  }
  c[kdec+1] = NumLoci+1;
  c[kdec+2] = 0;
}


///
/// Uses modified code by Donald Knuth and Glenn C. Rhoads to generate
/// all possible combinations from COMBOSTART to COMBOEND in intervals specified in COMBOINTERVAL
/// That is, N = COMBOEND, C = COMBOSTART. Summation of N choose C  -->  N choose N.
/// Code is confusing but fast
/// @return true if all combinations completed
///
bool KnuthComboGenerator::GenerateCombinations(){
  // Clear the ComboList!!! Very important, or Memory will overflow
  ComboList.clear();

  counter = 0;  // counter for the number of combinations created

  if(AlreadyStarted)   // If combination generator has already been started
  {
    counter--;
    goto resume;    // Ugh.. a GOTO... but it works!
  }
  else{
    kdec = ComboEnd;
    x = 0;
    j = kdec;
    
    if(NumLoci == ComboEnd){
      j=0;
    }
    
    c = new int[ComboEnd+3];

    for(int i=1; i <= kdec; i++)
    {
      c[i] = i;
    }
    c[kdec+1] = NumLoci+1;
    c[kdec+2] = 0;
  }

  AlreadyStarted = true;  // If it wasn't already started, it is now

  init:

      if(kdec == (ComboStart-1))
      {
        delete [] c;
        AlreadyStarted = false;
        counter--;
        return true;
      }

      for (int i=1; i <= kdec; i++)
      {  c[i] = i;  }
      c[kdec+1] = NumLoci+1;
      c[kdec+2] = 0;
      j = kdec;
      
      if(NumLoci == ComboEnd){
        j=0;
      }
      
      // Add a new combination to the ComboList
      ComboList.push_back(std::vector <unsigned int>());

  visit:
      for (int i=kdec; i >= 1; i--)
      {
    // Add an element to the new combination
    ComboList[counter].push_back(c[i]);
      }


      if(counter >= ComboInterval)  // If you exceed the interval limit
      {
          return false;
      }

  resume:
            // Add a new combination to the ComboList
      ComboList.push_back(std::vector <unsigned int>());

      counter++;
      if (j > 0) {x = j+1; goto incr;}

      if (c[1] + 1 < c[2])
         {
         c[1] += 1;
         goto visit;
         }

      j = 2;

   do_more:
      c[j-1] = j-1;
      x = c[j] + 1;
      if (x == c[j+1]) {j++; goto do_more;}

      if (j > kdec)
      {
          kdec--;

    // Remove the last empty combination before quitting
          ComboList.pop_back();
          goto init;
      }

   incr:
      c[j] = x;
      j--;
      goto visit;
}

///
/// Generates combinations in amount passed
/// @param new_interval Combinations to create
/// 
bool KnuthComboGenerator::GenerateCombinations(int new_interval){
  int old_interval = getComboInterval();
  SetComboInterval(new_interval);
  bool done = GenerateCombinations();
  SetComboInterval(old_interval);
  return done;
}


///
/// Advances parameters without adding any combinations
/// to the combination list
///
bool KnuthComboGenerator::AdvanceParameters(){
  counter = 0;  // counter for the number of combinations created

  if(AlreadyStarted)   // If combination generator has already been started
  {
    counter--;
    goto resume;    // Ugh.. a GOTO... but it works!
  }
  else{
    kdec = ComboEnd;
    x = 0;
    j = kdec;
    
    if(NumLoci == ComboEnd){
      j=0;
    }
    
    c = new int[ComboEnd+3];

    for(int i=1; i <= kdec; i++)
    {
      c[i] = i;
    }
    c[kdec+1] = NumLoci+1;
    c[kdec+2] = 0;
  }

  AlreadyStarted = true;  // If it wasn't already started, it is now

  init:

      if(kdec == (ComboStart-1))
      {
        delete [] c;
        AlreadyStarted = false;
        --counter;
        return true;
      }

      for (int i=1; i <= kdec; i++)
      {  c[i] = i;  }
      c[kdec+1] = NumLoci+1;
      c[kdec+2] = 0;
      j = kdec;
      
      if(NumLoci == ComboEnd){
        j=0;
      }
      

  visit:


      if(counter >= ComboInterval)  // If you exceed the interval limit
      {
          return false;
      }

  resume:
      counter++;
      if (j > 0) {x = j+1; goto incr;}

      if (c[1] + 1 < c[2])
         {
         c[1] += 1;
         goto visit;
         }

      j = 2;

   do_more:
      c[j-1] = j-1;
      x = c[j] + 1;
      if (x == c[j+1]) {j++; goto do_more;}

      if (j > kdec)
      {
          kdec--;

          goto init;
      }

   incr:
      c[j] = x;
      j--;
      goto visit;  
}

///
/// returns total number of combinations in set
/// @return number of combinations to construct
///
double KnuthComboGenerator::calc_combinations(){
  double total = 0.0;
  
  unsigned int endFact, n = NumLoci;
  double k, numerator, k_fact, result;
  unsigned int i;
  
  for(int combsize=ComboStart; combsize <= ComboEnd;
    combsize++){
    k = combsize;
    numerator = double(n);
    endFact = n-combsize;
    // calculate numerator
    for(i=NumLoci-1; i>endFact; i--){
      numerator *= i;
    }

    k_fact = combsize;
    for(i=combsize-1; i >=1; i--){
      k_fact *= i;
    }

    result = numerator / k_fact;
    total += result;
  }
  return total;
}

///
/// Reads checkpoint information into generator
/// @param instream Input to read from
///
void KnuthComboGenerator::ReadCheckpoint(ifstream& instream){
  int total_c, j, x, kdec, c;
  instream >> j >> x >> kdec;
  param_j(j);
  param_x(x);
  param_kdec(kdec);
  instream >> total_c;
  int i;
  for(i=0; i<total_c; i++){
    instream >> c;
    set_param_c(i,c);
  }
}
 
 ///
 /// Writes checkpoint information from generator
 /// @param outstream Stream for writing
 ///
void KnuthComboGenerator::WriteCheckpoint(ofstream& outstream){
  outstream << param_j() << " " << param_x() << " " << param_kdec() << endl;
  int total_c = get_combo_max()+3;
  outstream << total_c << endl;
  outstream << param_c(0);
  for(int i=1; i<total_c; i++){
    outstream << " " << param_c(i);
  }
  outstream << endl;
}

}

