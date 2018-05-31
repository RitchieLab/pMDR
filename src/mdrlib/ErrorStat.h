//ErrorStat.h

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

#ifndef __ERRORSTAT_H__
#define __ERRORSTAT_H__

#include "Stat.h"
#include <iostream>

using namespace std;

namespace mdr{
///
/// Calculates error rate
///
class ErrorStat: public Stat{

  public:
  
    ErrorStat(){set_analysis_name("Error");}
  
    virtual float analysis_value(unsigned int true_positives, unsigned int false_positives, unsigned int true_negatives, 
      unsigned int false_negatives){        
        return (false_positives + float(false_negatives)) /  (false_positives + false_negatives + true_negatives
          + true_positives) * 100;
    }
    
    virtual float convert_to_display(float calculated_value){return calculated_value;}
    
  private:

};

}

#endif
