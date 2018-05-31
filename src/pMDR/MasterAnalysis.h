//MasterAnalysis.h

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


#ifndef __MASTERANALYSIS_H__
#define __MASTERANALYSIS_H__

#include <Analysis.h>
#include "mpi.h"
#include "pComboTransfer.h"

using namespace std;
using namespace mdr;

/// Conducts master node portion of analysis
class MasterAnalysis: public Analysis{

  public:
  
    MasterAnalysis(){paramMover=NULL; srParameters=5;}
    
    ~MasterAnalysis(){}
    
    void run_master_analysis(Dataset& set, Config& config, int nproc,
      LogOutput& log_out, bool log_all);
 
    /// receive models from slaves
    void receive_models(int node);
 
    /// sets generator limits
    void set_combination_limits(int total_loci, int min_size, int max_size);
    
    void receive_checkpoint_models(MPI_Request* RequestParameterspending, int nproc);
 
    void send_continue_signal();
    
    void set_generator_type(string biofilename);
 
  private:

    pComboTransfer* paramMover;
    int srParameters;

};

#endif

