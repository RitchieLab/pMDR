//BinDataFileReader.cpp

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

#include "BinDataFileReader.h"
#include "dynamic_bitset/dynamic_bitset.hpp"

namespace mdr{

/// 
/// Opens and reads binary datafile
/// @param datafile string containing name of file
/// @param set Dataset to fill with datafile
///
void BinDataFileReader::read_datafile(string datafile, Dataset& set){

  ifstream::pos_type size;
  ifstream binfile(datafile.c_str(), ios::in|ios::binary);
  
  if(!binfile.is_open()){
    throw MDRExcept("Unable to open file: " + datafile);
  }
  
  // get version
  char header_info[5];
  binfile.read(header_info, 4);
  header_info[4] = '\0';
  
  string header = header_info;
  if(header.compare("CC  ") != 0){
    throw MDRExcept("\nERROR:  File '" + datafile + "' does not match expected binary case control format\n");
  }
  
  char endian[3];
  binfile.read(endian, 2);
  endian[2] = '\0';
  
  // get version integer
  unsigned int version;
  binfile.read((char*)(&version), sizeof(int));
  
  // get bits per block
  unsigned int bits_per_block;
  binfile.read((char*)(&bits_per_block), sizeof(unsigned int));
  
  // get flags
  unsigned int flags;
  binfile.read((char*)(&flags), sizeof(unsigned int));
 
  unsigned int num_strings = 2;
  // when first bit is set -- use 
  if(flags && 1){
    num_strings = 3;
    set.set_contains_missing(true);
  }
  else{
    set.set_contains_missing(false);
  }
  
  unsigned int num_inds;
  binfile.read((char*)(&num_inds), sizeof(unsigned int));
  set.set_num_inds(num_inds);
  
  // get number of chromosomes
  unsigned int num_chroms;
  binfile.read((char*)(&num_chroms), sizeof(unsigned int));
  
  unsigned int total_snps=0;
 
  // get chromosome sizes (in SNPs)
  vector<unsigned int> chrom_sizes;
  unsigned int chrom_size;
  for(unsigned int curr_chrom=0; curr_chrom < num_chroms; ++curr_chrom){
    binfile.read((char*)(&chrom_size), sizeof(unsigned int));
    chrom_sizes.push_back(chrom_size);
    total_snps += chrom_size;
  }
  
  set.set_num_loci(total_snps);
  
  // set up data size to be correct in the Dataset object
  set.initialize_data();
  
  unsigned int affectedTotal=0, unaffectedTotal=0;

  set.set_missing_value(3); 
  set.set_max_genotype_value(2);
  
  vector<boost::dynamic_bitset<> > alleles(3,boost::dynamic_bitset<>(0,0));
  
  boost::dynamic_bitset<> missing;
  if(num_strings == 2){
    boost::dynamic_bitset<> temp(total_snps, 0);
    alleles[2] = temp;
  }
  
  unsigned int status, allele_string, curr_ind, current_num_snps, curr_chrom;
  
  for(curr_ind=0; curr_ind < num_inds; ++curr_ind){
    binfile.read((char*)&status, sizeof(unsigned int));
    set.data[curr_ind][0] = status;
  
    if(set.data[curr_ind][0] == 1){
      affectedTotal++;
    }
    else if(set.data[curr_ind][0] == 0){
      unaffectedTotal++;
    }
  
    for(allele_string = 0; allele_string < num_strings; ++allele_string){
      current_num_snps = 0;
      for(curr_chrom=0; curr_chrom < num_chroms; ++curr_chrom){
        unsigned int totalBlocks = (chrom_sizes[curr_chrom] / bits_per_block) + 1;
        boost::dynamic_bitset<>::block_type raw[totalBlocks];
        
        binfile.read((char*)&raw, (sizeof(boost::dynamic_bitset<>::block_type)*totalBlocks));

        alleles[allele_string].append(&raw[0], &raw[totalBlocks]); 
        
        current_num_snps += chrom_sizes[curr_chrom];
        alleles[allele_string].resize(current_num_snps);
      }
    }
  
    // construct genotypes for the indvidual by adding the alleles and missing string together
    for(unsigned int curr_snp=0; curr_snp < total_snps; curr_snp++){
      set.data[curr_ind][curr_snp+1] = alleles[0][curr_snp] + alleles[1][curr_snp] + alleles[2][curr_snp];
    }
    
    for(allele_string=0; allele_string < num_strings; ++allele_string){
      alleles[allele_string].clear();
    }
    
  }
  
  binfile.close();

}


}
