//FileManager.h

#ifndef __FILEMANAGER_H__
#define __FILEMANAGER_H__

#include "BinDataFileReader.h"
#include "DatafileReader.h"

namespace mdr{

///
/// Manages reading of input data. <br>
/// Each format is specified in the enumeration.
///

/// Reads data files
class FileManager{

  public:
    /// Indicates format of file to process
    enum FileType{
      /// Enum for binary file
      Binary,
      /// Enum for text file
      Text
    };
  
    /// parses input file and adds data to Dataset 
    void parse_input_file(string identifier, Dataset & set,
      string format);
  
  private:
    
    /// returns pointer to specific filehandler for use in this operation
    FileReader* get_filereader(FileType format);
    
    /// returns reader type based on string specifier
    FileManager::FileType set_reader_type(string input_type);
    
};

}

#endif
