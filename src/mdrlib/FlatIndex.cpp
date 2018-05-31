//FlatIndex.cpp

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

#include <iostream>
#include "FlatIndex.h"

namespace mdr{

//  int FlatIndex::num_genos_per_locus=0;
// vector<vector< int> > FlatIndex::includedIndexes;

///
///  Convert multi-dimensional array indexes
///   into linear array index
/// @param indexes -- vector containing index values
/// @return index into linear array
///
 int FlatIndex::flatten_indexes(vector<unsigned int> & indexes){

   int linear_index = 0;
  for( unsigned int i = 0; i < indexes.size(); i++){
    linear_index += indexes[i];

    if(i+1 < indexes.size())
      linear_index *= num_genos_per_locus;
  }

  return linear_index;
}


///
///  Convert multi-dimensional array indexes
///   into linear array index
/// @param indexes -- vector containing index values
/// @return index into linear array
///
 int FlatIndex::flatten_indexes(vector<int> & indexes){

   int linear_index = 0;
  for( unsigned int i = 0; i < indexes.size(); i++){
    linear_index += indexes[i];

    if(i+1 < indexes.size())
      linear_index *= num_genos_per_locus;
  }

  return linear_index;
}

///
/// Converts the linear index value into
/// the vector (genotype) corresponding to that index
/// @param index -- linear index for the genotype
/// @param numloci -- number of loci in genotype
/// @return vector where each value is one locus
///
vector< int> FlatIndex::decode_index(int index, int numloci){
  vector< int> genotype(numloci,0);
   long locusIndex;
  for(int i=numloci; i; --i){
    locusIndex = i-1;
    genotype[locusIndex] = index % num_genos_per_locus;
    index /= num_genos_per_locus;
  }
  return genotype;
}


///
/// Returns total size of array needed
/// @return size
///
int FlatIndex::get_size_array( int numLoci){
  int size=1;

  for( int i=1; i <= numLoci; i++)
    size *= num_genos_per_locus;
  return size;
}


///
/// Sets the included indexes.  These are the indexes
///  that are totalled when determining number of
///  individuals classified correctly.  When the option
///  to ignore
/// @param startCombSize
/// @param endCombSize
/// @param ignoreMissing
/// @param missingValue
///
void FlatIndex::set_included_indexes(int startCombSize, int endCombSize,
  bool ignoreMissing, int missingValue){
  // add empty vectors to the 2-D vector
  includedIndexes.clear();
  
  for(int currVec=0; currVec <= endCombSize; currVec++){
    includedIndexes.push_back(vector< int>(0,0));
  }

  // when not ignoring individuals with missing
  // genotypes, all indexes (cells) are included
  if(!ignoreMissing){
    int sizeVector, currCell;
    for(int currCombSize=startCombSize; currCombSize <= endCombSize;
      currCombSize++){
      sizeVector = FlatIndex::get_size_array(currCombSize);
      for(currCell=0; currCell<sizeVector; currCell++)
        includedIndexes[currCombSize].push_back(currCell);
    }
  }
  else{
    for(int currCombSize=startCombSize; currCombSize <= endCombSize;
      currCombSize++)
      set_non_missing_indexes(currCombSize, missingValue);
  }
  
}



void FlatIndex::set_non_missing_indexes(int comboSize,  int missingValue){

   int * lower_index = new  int[comboSize];
   int * upper_index = new  int[comboSize];

  for(int currLoc=0; currLoc<comboSize; currLoc++){
    lower_index[currLoc] = 0;
    upper_index[currLoc] = num_genos_per_locus;
  }

   int * indexes = new  int[comboSize];
  int cur_depth = 0;
  int max_depth = comboSize;

  indexes[cur_depth] = lower_index[cur_depth];
  vector< int> loc_indexes;
  bool containsMissing = false;
  int i;

  while(1){
    if( indexes[cur_depth] < upper_index[cur_depth] ) {
      if( cur_depth == max_depth - 1 ) {
        // need to determine whether to include
        // or not here
        loc_indexes.clear();
        containsMissing = false;
        for(i=0; i<=cur_depth; i++){
          if(missingValue == indexes[i]){
            containsMissing = true;
            break;
          }
          loc_indexes.push_back(indexes[i]);
        }

        if(!containsMissing){
          includedIndexes[comboSize].push_back(FlatIndex::flatten_indexes(loc_indexes));
        }
        indexes[cur_depth]++;
      }
      else{
        cur_depth++;
        indexes[cur_depth]=lower_index[cur_depth];
      }
    }
    else{
      if(cur_depth > 0){
        indexes[--cur_depth]++;
              }
      else{
        break;
      }
    }
  }
  delete [] indexes;

  delete [] lower_index;
  delete [] upper_index;
}



}

