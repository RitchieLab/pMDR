//FileManager.cpp

#include "FileManager.h"

namespace mdr{

///
///Returns pointer to appropriate FileReader
///based on format passed<br>
///Caller is responsible for freeing memory
///of dynamically created FileReader
///@param format enumeration with format
///@return FileReader pointer
///@throws MDRExcept when no matching format
///
FileReader * FileManager::get_filereader(FileManager::FileType format){
  FileReader * newFileHandler;
  switch(format){
    case Binary:
      newFileHandler = new BinDataFileReader;
      break;
    case Text:
      newFileHandler = new DatafileReader;
      break;
    default:
      throw MDRExcept(format + " is not a valid format.");
  };
  return newFileHandler;
}


///
/// Returns input type based on string
/// @param input_type
/// @return enumerated value matching the string
/// @throws MDRExcept when no matching format
/// 
FileManager::FileType FileManager::set_reader_type(string input_type){

  FileManager::FileType ftype;

  if(input_type.compare("TEXT")==0){
    ftype = Text;
  }
  else if(input_type.compare("BINARY")==0){
    ftype = Binary;
  }
  else{
    throw MDRExcept("No input handling for " + input_type + "\n");
  }

  return ftype;
}


///
///Parses the data file and fills this object with data from file
///@param identifier name of data file or database ID
///@param dataholder DataHolder to fill with
///data from the datafile
///@param format string with format specifier
///@return none
///@throws DataExcept
///
void FileManager::parse_input_file(string identifier, Dataset & set,
  string format){
  FileReader * fr = get_filereader(set_reader_type(format));
  fr->read_datafile(identifier, set);
  delete fr;
}



}
