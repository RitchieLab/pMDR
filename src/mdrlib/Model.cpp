//Model.cpp

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

#include "Model.h"

namespace mdr{

///
/// Clears stats for new totals
/// 
void Model::reset_stats(){

	training.classhigh=0;
	training.classlow=0;
	training.misclasshigh=0;
	training.misclasslow=0;
	training.correctclass=0;
	training.incorrectclass=0;
	training.totaltiecells=0;
	training.totalnotrecognized=0;
	training.totalmissing=0;
	training.error=0;
	training.balanced_error=0;
	
	testing.classhigh=0;
	testing.classlow=0;
	testing.misclasshigh=0;
	testing.misclasslow=0;
	testing.correctclass=0;
	testing.incorrectclass=0;
	testing.totaltiecells=0;
	testing.totalnotrecognized=0;
	testing.totalmissing=0;
	testing.error=0;
	testing.balanced_error=0;
}


///
/// Resizes risk vector to correct size
/// @param vector_size 
/// @return
///
void Model::resize_cell_status_vector(int vector_size){
  cell_status.assign(vector_size, -1);
}


}

