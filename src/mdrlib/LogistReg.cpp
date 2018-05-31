//LogistReg.cpp

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

#include "LogistReg.h"
#include "KnuthComboGenerator.h"
#include <math.h>

namespace mdr{

///
/// Constructor
///
LogistReg::LogistReg(){
  initialize();
}


///
/// initializes values
/// @return none
///
void LogistReg::initialize(){

  modelSize = 0;
  numberGenosPerLocus = 0;

  maxLocusValue = 1;
  missingValue = maxLocusValue+1;
  LociComboLimit = 3;
  LociComboMin = 1;
  maxIterations = 20;
  defaultComboInterval = 10000;
  interaction_included = true;
  
  PiD2 = 3.141592653589793 / 2;
  
  initialize_interactions();
  
}


///
/// Generates and stores interaction lists for producing interaction
/// terms when doing multi-locus models
/// @return
///
void LogistReg::initialize_interactions(){
 
  vector<vector<unsigned int> > blank_vector;

  interaction_lists.clear();
  
  // indexes 0 and 1 contain empty vectors as placeholders
  // so that index of the interaction list will match
  // number of loci in the model
  interaction_lists.push_back(blank_vector);
  interaction_lists.push_back(blank_vector);
  
  unsigned int comboIndex;
  for(unsigned int curr_num=2; curr_num <= LociComboLimit; curr_num++){
    interaction_lists.push_back(blank_vector);
    
    KnuthComboGenerator generator;
    generator.ComboEnds(2, curr_num);
    generator.SetLoci(curr_num);
    generator.SetComboInterval(defaultComboInterval);
    bool done=true;
    
    do{
      done = generator.GenerateCombinations();
      for(comboIndex=0; comboIndex < generator.ComboList.size(); comboIndex++){
for(unsigned int z=0; z<generator.ComboList[comboIndex].size(); z++){
  generator.ComboList[comboIndex][z]--;
}
        
        interaction_lists[curr_num].push_back(generator.ComboList[comboIndex]);
      }      
      
      
    }while(!done);
  }

}



///
/// Sizes and clears arrays for holding data
/// @param currModelSize size of model
/// @return
///
void LogistReg::initialize_summary(unsigned int currModelSize, bool include_full_interaction){

  unsigned int array_size = indexConverter.get_size_array(currModelSize);    
  // check to see if need to resize all parameters
  if(currModelSize != modelSize){
    modelSize = currModelSize;

    // establish vector for holding individual totals
    summary_data.assign(array_size, vector<int>(currModelSize,0));
      zero_summary(array_size, currModelSize, include_full_interaction);
  }
  else{
    // zero out totals for summary data
      zero_summary(array_size, currModelSize, include_full_interaction);
  }

}



///
///  Sets genotype for single array, adds interaction terms if
///  any needed and requested by user, zeroes out totals
///  @param array_size unsigned int number of rows in summary
///  @param currModelSize unsigned int current model size to work with
///  @param include_full_interaction If not true then exclude, 
///  @return 
///
void LogistReg::zero_summary(unsigned int array_size, unsigned int model_size, bool include_full_interaction){

  unsigned int product=1;
  unsigned int number_of_non_interaction;
  
  for(unsigned int sub_array_index=0; sub_array_index < array_size; sub_array_index++){ 
    summary_data[sub_array_index] = indexConverter.decode_index(sub_array_index, model_size);
    number_of_non_interaction = summary_data[sub_array_index].size();
  
    // this version can handle any size model
    if(summary_data[sub_array_index].size() > 1 && interaction_included){
      for(int i=int(interaction_lists[model_size].size())-1; i >= 0; i--){
      product=1;
        for(int j=int(interaction_lists[model_size][i].size())-1; j >= 0; j--){

          product*= summary_data[sub_array_index][interaction_lists[model_size][i][j]];
        }
        summary_data[sub_array_index].push_back(product);
      }
    }
  
    // when not full interaction have to remove the full interaction term
    if(!include_full_interaction){
      summary_data[sub_array_index].pop_back();
    }
      
  
    summary_data[sub_array_index].push_back(0);
    summary_data[sub_array_index].push_back(0);
  }
  
}



///
/// Fills summary vector with totals for use 
/// in logistic regression routine
/// @param dataset vector containing data to analyze
/// @param status vector containing matching status for data
/// @param max_data_value maximum value in the dataset
/// @param any_missing True if data has any missing data
/// @param include_full_interaction True to use full model and false to use the reduced
/// model that excludes the highest order interaction
/// @return
///
void LogistReg::summarize_data(vector<vector<unsigned int> >& dataset, vector<unsigned int>& status,
  int max_data_value, bool any_missing, bool include_full_interaction){
  unsigned int combSize = dataset.size();

  unsigned int nGenosPerLocus = max_data_value+1;

  // set up the converter that can translate multiple dimensional vectors
  // into single positions  
  if(numberGenosPerLocus != nGenosPerLocus){
    if(any_missing)
      indexConverter.set_genos_per_locus(nGenosPerLocus+1); // assume zero based
    else
      indexConverter.set_genos_per_locus(nGenosPerLocus);
    indexConverter.set_included_indexes(LociComboMin, LociComboLimit, 
      !any_missing, nGenosPerLocus);
      
      numberGenosPerLocus = nGenosPerLocus;
  }
  
  includedIndexes = indexConverter.get_included_indexes();
  
  // clear and initialize arrays
  // set model 
  initialize_summary(combSize, include_full_interaction);
  
  unsigned int currInd, currLoc;

  // establish a vector of correct size that can be used to distribute individuals
  vector<int> genotype(combSize,0);

  // determine indexes for unaffected and affected totals
  unsigned int unaffIndex = summary_data[0].size()-2;
  
  // the dataset is orientend with rows being loci and columns the individual
  // genotype value
  unsigned int numInds = dataset[0].size();
  
  // add to summary totals for each genotype
  // use current model to convert genotypes
  for(currInd=0; currInd < numInds; currInd++){
    for(currLoc=0; currLoc < combSize;
      currLoc++){
      genotype[currLoc] = dataset[currLoc][currInd];
    }

    // increment count based on status of individual
    summary_data[indexConverter.flatten_indexes(genotype)][unaffIndex+status[currInd]]++;
  }
}


///
/// Returns -2 Log-likelihood for the vectors passed
/// @param data 2-D vectors containing data to fit
/// @param status vector contains status values for data vectors
/// @return -2 Log-likelihood for the model
/// 
float LogistReg::run_lr(vector<vector<unsigned int> >& data, vector<unsigned int>& status, bool any_missing,
  int max_data_value, bool include_full_interaction){
  
  // clear and initialize arrays
  // set model 
  summarize_data(data, status, max_data_value, any_missing, include_full_interaction);

  // determine number of genotype combination values
  unsigned int numLoci = data.size();
vector<int> includedCells = includedIndexes[numLoci];

  unsigned int nC = includedCells.size();
  // number of coefficients is equal to the main effects + the interactions
  unsigned int nR = summary_data[0].size()-2;
  unsigned int nP = nR + 1;
  unsigned int nP1 = nP + 1;
  unsigned int sY0 = 0;
  unsigned int sY1 = 0;
  unsigned int sC = 0;
  double x;
  double v, xij, s,q;
  
  unsigned int unaffIndex = summary_data[0].size()-2;
  unsigned int affIndex = unaffIndex+1;

  vector<double> X(nC * ( nR + 1 ),0);
  vector<unsigned int> Y0(nC,0);
  vector<unsigned int> Y1(nC,0);
  vector<double> xM(nR+1,0.0);
  vector<double> xSD(nR+1, 0.0);
  vector<double> Par(nP);
  vector<double> SEP(nP);
  vector<double> Arr(nP * nP1); 
  
  unsigned int i,j,k;
  for(i=0; i<nC; i++){
    X[ix(i,0,nR+1)] = 1;
    // store predictor values
    for(j=1; j<=nR; j++){
      X[ix(i,j,nR+1)] = summary_data[includedCells[i]][j-1];
    }
    Y0[i] = summary_data[includedCells[i]][unaffIndex];
    sY0 += Y0[i];
    Y1[i] = summary_data[includedCells[i]][affIndex];
    sY1 += Y1[i];
    
    sC += Y0[i] + Y1[i];
    
    for(j=1; j<= nR; j++){
      x = X[ix(i,j,nR+1)];
      xM[j] += (Y0[i] + Y1[i])*x;
		  xSD[j] += (Y0[i] + Y1[i])*x*x;
    }
  }

  // calculate mean and standard deviation
  for (j = 1; j<=nR; j++) {
    xM[j]  = xM[j]  / sC;
    xSD[j] = xSD[j] / sC;
    xSD[j] = sqrt( fabs( xSD[j] - xM[j] * xM[j] ) );
  }

  xM[0] = 0; 
  xSD[0] = 1;

  for (i = 0; i<nC; i++) {
    for (j = 1; j<=nR; j++) {
      X[ix(i,j,nR+1)] = ( X[ix(i,j,nR+1)] - xM[j] ) / xSD[j];
    }
  }

  Par[0] = log( double(sY1) / sY0 ); // use natural log of the ratio  
  for (j = 1; j<=nR; j++) { // zero out all the others
    Par[j] = 0;
  }
  
  double LnV=0,Ln1mV=0, LLn=0;
  double LLp = 2e+10;
  double LL = 1e+10;
  unsigned int numIterations = 0;  

  while( fabs(LLp-LL)>0.0000001 ) {
    if(++numIterations > maxIterations){
      break;
    }
 
    LLp = LL;
    LL = 0; 
  
    for (j = 0; j<=nR; j++) {
      for (k = j; k<=nR+1; k++) {
        Arr[ix(j,k,nR+2)] = 0;
      }
    }   
	  
    for (i = 0; i<nC; i++) {
      v = Par[0];
      for (j = 1; j<=nR; j++) {
        v = v + Par[j] * X[ix(i,j,nR+1)];
      }
      if( v>15 ) { LnV = -exp(-v); Ln1mV = -v; q = exp(-v); v=exp(LnV); }
      else { 
        if( v<-15 ) { LnV = v; Ln1mV = -exp(v); q = exp(v); v=exp(LnV); }
        else { v = 1 / ( 1 + exp(-v) ); LnV = log(v); Ln1mV = log(1-v); q = v*(1-v); }
      }
      LL = LL - 2*Y1[i]*LnV - 2*Y0[i]*Ln1mV;
      
      // when all samples have same genotype LL will be nan and
      // so return 1.0
      if(isnan(LL))
        return 1.0;

      for (j = 0; j<=nR; j++) {
        xij = X[ix(i,j,nR+1)];
        Arr[ix(j,nR+1,nR+2)] = Arr[ix(j,nR+1,nR+2)] + xij * ( Y1[i] * (1 - v) + Y0[i] * (-v) );
        for (k=j; k<=nR; k++) {
          Arr[ix(j,k,nR+2)] = Arr[ix(j,k,nR+2)] + xij * X[ix(i,k,nR+1)] * q * (Y0[i] + Y1[i]);
        }
      }
    }	    
	  
    if( LLp==1e+10 ) { LLn = LL;}
	  
    for (j = 1; j<=nR; j++) {
      for (k=0; k<j; k++) {
        Arr[ix(j,k,nR+2)] = Arr[ix(k,j,nR+2)];
      }
    }
		
    for (i=0; i<=nR; i++) { 
      s = Arr[ix(i,i,nR+2)]; 
      Arr[ix(i,i,nR+2)] = 1;
      for (k=0; k<=nR+1; k++) { 
        Arr[ix(i,k,nR+2)] = Arr[ix(i,k,nR+2)] / s; 
      }
      for (j=0; j<=nR; j++) {
        if (i!=j) { 
          s = Arr[ix(j,i,nR+2)]; Arr[ix(j,i,nR+2)] = 0;
          for (k=0; k<=nR+1; k++) {
            Arr[ix(j,k,nR+2)] = Arr[ix(j,k,nR+2)] - s * Arr[ix(i,k,nR+2)];
          }
        }
      }
    }
    for( j=0; j<=nR; j++) {
      Par[j] = Par[j] + Arr[ix(j,nR+1,nR+2)];
    }	  
  }

  return LL;
}



///
/// Returns index into array
/// @param j
/// @param k
/// @param nCols
/// @returns index
/// 
unsigned int LogistReg::ix(int j,int k,int nCols){
  return j * nCols + k;
}


} // end namespace mdr

