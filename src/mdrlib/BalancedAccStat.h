//BalancedAccStat

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

#ifndef __BALANCEDACCSTAT_H__
#define __BALANCEDACCSTAT_H__

#include "Stat.h"

using namespace std;

namespace mdr{
///
/// Abstract base class that calculates stat to be used in ranking models
///
class BalancedAccStat: public Stat{

  public:
  
    BalancedAccStat(){set_analysis_name("Balanced Acc");}
  
    virtual float analysis_value(unsigned int true_positives, unsigned int false_positives, unsigned int true_negatives, 
      unsigned int false_negatives){        
        return (1 - .5 * (float(true_positives)/(true_positives+false_positives) +
          (float(true_negatives)/(true_negatives+false_negatives)))) * 100;
    }
    
    
    virtual float convert_to_display(float calculated_value){return 100 - calculated_value;}
    
  private:


};



}

#endif
