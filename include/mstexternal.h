#ifndef mstexternal_h
#define mstexternal_h

#include "mstsystem.h"

// For classes that need to interact with non-MST programs

/**
 STRIDE is a standard program for computing the secondary structure of protein residues. It uses a combination of dihedral angles and
 hydrogen bonding atoms to assign each residue. The following class is really just a wrapper that runs STRIDE through an external call
 and then parses the output. For convenience, the STRIDE secondary structure class key is provided here:
 
 H      Alpha helix
 G      3-10 helix
 I      PI-helix
 E      Extended conformation
 B or  b   Isolated bridge
 T      Turn
 C      Coil (none of the above)
 
 http://webclu.bio.wzw.tum.de/stride/
 Frishman D, Argos P. Knowledge-Based Protein Secondary Structure Assignment Proteins: Structure, Function, and Genetics 23:566-579 (1995)
 */

struct resInfo {
  //line numbers in file (1-indexed)
  string resName; //6-8
  string chainID; //10
  int resNum;     //12-15
  string ssType;  //25-27
};

class strideInterface {
public:
  /**
   @param _strideBin Path to a STRIDE binary
   @param _S Pointer to a structure object
  */
  strideInterface(string _strideBin, Structure* _S = nullptr) : strideBin(_strideBin), S(_S) {};
  
  /**
   @param _S Path to a structure (will reset the object state)
   */
  void setStructure(Structure* _S) {
    S = _S;
    data.clear();
  }
  
  /**
   Compute the STRIDE classifications of the current structure
   
   @param options Additional options that will be passed to the command line call
   */
  void computeSTRIDEClassifications(string options = "");
  
  /**
   @param strict If true, will throw an error if the stride data doesn't match the residues in the structure. Otherwise will interpolate.
   @return The single letter classifications of the residues in the current structure. The order is consistent with the structure.
   */
  vector<string> getSTRIDEClassifications(bool strict = false);
    
private:
  string strideBin = "";
  Structure* S = nullptr;
  vector<resInfo> data;
};

#endif /* mstexternal_h */
