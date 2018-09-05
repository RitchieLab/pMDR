//SlaveAnalysis.h

//     Copyright 2008 Marylyn Ritchie

//     This file is part of pMDR.
// 
//     pMDR is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
// 
//     pMDR is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License
//     along with pMDR.  If not, see <http://www.gnu.org/licenses/>.


#ifndef __SLAVEANALYSIS_H__
#define __SLAVEANALYSIS_H__

#include "mpi.h"
#include <Analysis.h>
#include "pComboTransfer.h"

using namespace std;
using namespace mdr;

/// Conducts master node portion of analysis
class SlaveAnalysis: public Analysis{

  public:
  
    SlaveAnalysis(){srParameters=5; paramMover=NULL;}
    
    ~SlaveAnalysis();
  
    /// runs analysis
    void run_slave_analysis(Dataset& set, Config& config);
       
    /// receive combination generator parameters from master
    int receive_slave_parameters(int request_number);

    /// send models to master
    void send_models();    
     
    void set_always_included(vector<unsigned int>& loci){always_included = loci;}
    
    void set_combination_limits(int total_loci, int min_size, int max_size);
  
    void receive_continue_signal();
  
    void set_rank(int rank){myrank=rank;}
    
    void set_generator_type(string biofilename, int idSize);
  
  private:
    int myrank, srParameters;
    pComboTransfer* paramMover;

};

#endif

