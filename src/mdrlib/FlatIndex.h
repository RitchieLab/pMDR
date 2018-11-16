//FlatIndex.h

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

#ifndef FLATINDEX_H
#define FLATINDEX_H

#include <vector>
#include<iostream>
using namespace std;

namespace mdr{

///
/// Contains functionality for converting a multi-dimensional
/// array into a linear array.
/// Adapted from Nate Barney's C MDR library
///

class FlatIndex{

  public:
    int flatten_indexes(vector<unsigned int> & indexes);
    int flatten_indexes(vector<int> & indexes);

    vector<int> decode_index(int index, int numloci);

    void set_genos_per_locus(int numGenos)
      {num_genos_per_locus = numGenos;}

    int get_genos_per_locus(){return num_genos_per_locus;}

    unsigned int get_size_array(int numLoci);

    vector<int> & get_valid_indexes(int numLoci){return includedIndexes[numLoci];}

    void set_included_indexes(int startLoc, int endLoc,
      bool ignoreMissing, int missingValue);

    void set_non_missing_indexes(int comboSize, int missingValue);

    vector<vector<int> > get_included_indexes(){return includedIndexes;}

  private:
    int num_genos_per_locus;
    vector<vector<int> > includedIndexes;
};

}

#endif // FLATINDEX_H

