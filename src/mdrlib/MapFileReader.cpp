//MapFileReader.cpp

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

#include "MapFileReader.h"
#include "Stringmanip.h"

namespace mdr{

void MapFileReader::parse_map_file(string mapfile, Dataset& set){

  std::ifstream mapStream(mapfile.c_str(), ios::in);
  if(!mapStream.is_open()){
    throw MDRExcept("ERROR: Unable to open " + mapfile + "\n");
  }

  string line;

  string snpID;
  string chrom;
  unsigned int pos;
  int var_count=0;

  while(!mapStream.eof()){
    getline(mapStream, line);

    if(line.find_first_of("0123456789") == string::npos){
      continue;
    }

    stringstream ss(line);

    ss >> chrom >> snpID >> pos;
    set.add_snp_name(snpID);
    set.add_chrom(chrom);
    set.add_pos(pos);
    var_count++;
  }
  mapStream.close();
  if(var_count != set.get_num_loci()){
    throw MDRExcept("ERROR: " + Stringmanip::itos(var_count) + " variables in " + mapfile + 
      " doesn't match number of variables in dataset file (" + Stringmanip::itos(set.get_num_loci()) + ")\n");
  }

}


}

