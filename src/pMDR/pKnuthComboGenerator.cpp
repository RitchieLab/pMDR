//pKnuthComboGenerator.cpp

#include "pKnuthComboGenerator.h"

namespace mdr{

///
/// For sending combo parameters.  
/// @param node 
/// @param srParameters
///
void pKnuthComboGenerator::SendParameters(int node, int srParameters){

  int numParams = get_combo_max() + 8;
  int * temp = new int[numParams];

  temp[0] = 0;
  temp[1] = param_AlreadyStarted();
  temp[2] = param_j();
  temp[3] = param_x();
  temp[4] = param_kdec();

  for(int i=0; i < get_combo_max()+3; i++){
    temp[i+5] = param_c(i);
  }
  
  MPI_Send(temp, numParams, MPI_INT, node, srParameters, MPI_COMM_WORLD);

  delete [] temp;
  
}

///
/// Sends terminate signal
///
void pKnuthComboGenerator::SendStop(int node, int srParameters){
  int numParams = get_combo_max() + 8;

  int* temp = new int[numParams];
  temp[0] = 1;
  temp[1] = 0;
  temp[2] = 0;
  temp[3] = 0;
  temp[4] = 0;
  
  for(int i=0; i < get_combo_max()+3; i++){
    temp[i+5] = 0;
  }

  MPI_Send(temp, numParams, MPI_INT, node, srParameters, MPI_COMM_WORLD);

  delete [] temp;   
}


///
/// Sends terminate signal
///
void pKnuthComboGenerator::SendModelSignal(int node, int srParameters){
  int numParams = get_combo_max() + 8;

  int* temp = new int[numParams];
  temp[0] = -1;
  temp[1] = 0;
  temp[2] = 0;
  temp[3] = 0;
  temp[4] = 0;
  
  for(int i=0; i < get_combo_max()+3; i++){
    temp[i+5] = 0;
  }

  MPI_Send(temp, numParams, MPI_INT, node, srParameters, MPI_COMM_WORLD);

  delete [] temp;   
}



///
/// For receiving combo parameters
/// @return false when stop signal received
///
int pKnuthComboGenerator::ReceiveParameters(int srParameters){
  
  MPI_Status status;

  int numParams = get_combo_max() + 8;

  int * temp = new int[numParams];

  MPI_Recv(temp, numParams, MPI_INT, 0, srParameters, MPI_COMM_WORLD, &status);
  
  int signal = temp[0];
  
  if(signal != 0){
    delete [] temp;
    return signal;
  }
  
  param_AlreadyStarted(temp[1]);
  param_j(temp[2]);
  param_x(temp[3]);
  param_kdec(temp[4]);
  
  for(int i=0; i < get_combo_max()+3; i++){
    set_param_c(i, temp[i+5]);
  }

  delete [] temp;
  
  return signal;
  
}


}
