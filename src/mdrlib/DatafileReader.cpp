//DatafileReader.cpp

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

#include "DatafileReader.h"

namespace mdr{

/// 
/// Opens and reads datafile
/// @param datafile string containing name of file
/// @param set Dataset to fill with datafile
///
void DatafileReader::read_datafile(string datafile, Dataset& set){
  count_total(datafile, set);
  
  // initialize data within dataset
  set.initialize_data();
  
  // fills dataset
  fill_set(datafile, set);  

}



void DatafileReader::fill_set(string datafile, Dataset& set){
  
  int read_as_int;
  
  int num_loci = set.get_num_loci();
  int num_inds = set.get_num_inds();
  int max_locus_value = set.get_max_locus_value();
  int missing_value = set.get_missing_value();
  int resetmissing = max_locus_value+1;
  
  // Open the data file 
  ifstream infile;
  infile.open(datafile.c_str());
  
  int affectedTotal=0, unaffectedTotal=0;
  
  // Read the data from the file into the data array
  for(int i = 0; i < num_inds; i++){
  
    // if IDs are included skip the first column
    // MDR doesn't need ID information
    if(ids_included){
      infile >> read_as_int;
    }
  
    infile >> read_as_int;
    set.data[i][0] = read_as_int;
    if(set.data[i][0] == 1){
      affectedTotal++;
    }
    else if(set.data[i][0] == 0){
      unaffectedTotal++;
    }
    else{
      throw MDRExcept("All individual status (first column) must be coded as 0 or 1 in " + datafile);
    }
    for(int j = 1; j <= num_loci; j++){
      infile >> read_as_int;

      set.data[i][j] = read_as_int;
      
      //check that all data are between min and maximum values
      if(set.data[i][j] != missing_value && (set.data[i][j] < 0 || set.data[i][j] > max_locus_value)){
        ostringstream maxLoc, missVal;
        maxLoc << max_locus_value;
        missVal << missing_value;
        throw MDRExcept("All genotypes must be between 0 and the MAXLOCUSVALUE (" + 
          maxLoc.str() + ") or equal to the value of MISSING parameter (" +
          missVal.str() + ")");
      }
      
     
      // reset any missing values to be one more than the maxlocivalue
      if(set.data[i][j] == missing_value){
        set.data[i][j] = resetmissing;
        set.set_contains_missing(true);
      }
    }
  }
  infile.close();  
  
  // calculate set threshold and pass to the dataset
  set.set_entire_threshold(float(affectedTotal) / unaffectedTotal);
  
}


/// 
/// Opens sets totals in the Dataset
/// @param datafile string containing name of file
/// @param set Dataset 
///
void DatafileReader::count_total(string datafile, Dataset& set){

  ifstream datastream;
  int linecount = 0;
  string line;      
  
  // default is no IDs in data file
  ids_included = false;
  
    // 8/1/05 smd -- new code for determining rows and columns
  datastream.open(datafile.c_str());
  if(!datastream.is_open()){
    throw MDRExcept("\nFILENAME: " + datafile + " cannot open.\n");
  }

  getline(datastream, line);
  int colcount=0;
  string out;

  string::size_type lastNumberPos = line.find_last_of("0123456789");

  if(lastNumberPos != string::npos){
    line = line.substr(0,lastNumberPos);
  }
  stringstream ss(line);
  while(!ss.eof()){
    ss >> out;
    colcount++;
  }
  set.set_num_loci(colcount-1);
  
  while (!datastream.eof()){
    // don't count blank lines
    if(line.find_first_of("0123456789") != string::npos){
      linecount++;
    }
    stringstream subss(line);
    subss >> out;
    if((out[0] != '0' and out[0] != '1') or out.length() > 1){
      ids_included = true;
    }
    getline(datastream, line);
  }
  
  if(ids_included){
    set.set_num_loci(colcount-2);
  }
  
  set.set_num_inds(linecount);

  datastream.close();
}


}

