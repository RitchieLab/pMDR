//pBioComboGenerator.cpp

#include "pBioComboGenerator.h"

namespace mdr{


///
/// For sending combo parameters.  For the biofilter reading, only need
/// the current index to get the next set from.
/// @param node 
/// @param srParameters
///
void pBioComboGenerator::SendParameters(int node, int srParameters){
  unsigned int * params = new unsigned int[2];
  
  params[0] =1;
  params[1] = _current_index;

  MPI_Send(params, 2,  MPI_UNSIGNED, node, srParameters, MPI_COMM_WORLD);
  delete [] params;
}


///
/// For sending stop signal parameters.
/// @param node
/// @param srParameters
///
void pBioComboGenerator::SendStop(int node, int srParameters){
  unsigned int * params = new unsigned int[2];
  
  params[0] = 0;
  params[1] = 0;

  MPI_Send(params, 2,  MPI_UNSIGNED, node, srParameters, MPI_COMM_WORLD);
  delete [] params;
}


///
/// For sending signal to slaves to send models to master.
/// @param node
/// @param srParameters
///
void pBioComboGenerator::SendModelSignal(int node, int srParameters){
  unsigned int * params = new unsigned int[2];
  
  params[0] = 0;
  params[1] = 1;

  MPI_Send(params, 2,  MPI_UNSIGNED, node, srParameters, MPI_COMM_WORLD);
  delete [] params;
}


///
/// For receiving combo parameters
/// @return 1 when stop signal received and -1 for returning models to master
///
int pBioComboGenerator::ReceiveParameters(int srParameters){
  unsigned int * params = new unsigned int[2];
  
  MPI_Status status;
  
  int return_value = 0;
  
  MPI_Recv(params, 2, MPI_UNSIGNED, 0, srParameters, MPI_COMM_WORLD, &status);
  
  
  if(params[0]){
    _current_index = params[1];
  }
  else if(params[1]){
    return_value = -1;
  }
  else{
    return_value = 1;
  }
  
  delete [] params;
  return return_value;
}

}
