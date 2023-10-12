#include "mstexternal.h"

void strideInterface::computeSTRIDEClassifications(string options) {
  // External system call to run STRIDE
  string pdbPath = S->getName();
  string pdbName = MstSys::splitPath(pdbPath,1);
  string strideOutFile = pdbName+"_stride.out";
  string structurePDBPath = pdbName+".pdb";
  S->writePDB(structurePDBPath);
  MstSys::csystem(strideBin + " " + structurePDBPath + " " + options + " -f" + strideOutFile);
  
  // Read in the STRIDE output file and parse
  vector<string> lines = MstUtils::fileToArray(strideOutFile);
  for (string line : lines) {
    string lineType = line.substr(0,3);
    if (lineType == "ASG") {
      resInfo vals; data.push_back(vals);
      resInfo& currentLineData = data.back();
      currentLineData.resName = line.substr(5,3);
      currentLineData.chainID = line.substr(9,1);
      currentLineData.resNum = MstUtils::toInt(line.substr(11,4));
      currentLineData.ssType = line.substr(24,1);
      
      if (currentLineData.ssType.empty()) {
        cout << "Warning: STRIDE classification field empty, replacing with X" << endl;
        currentLineData.ssType = "X";
      }
    }
  }
  MstUtils::assertCond(data.size() <= S->residueSize(),"There should be no more STRIDE classifications than residues in the structure","computeSTRIDEClassifications");
  // Delete the files
  MstSys::crm(strideOutFile); MstSys::crm(structurePDBPath);
}

vector<string> strideInterface::getSTRIDEClassifications(bool strict) {
  vector<string> ssClassifications;
  
  // Check for different conditions that indicate STRIDE encountered an error
  if (S->residueSize() < 5) {
    if (strict) MstUtils::error("STRIDE will ignore chains with less than 5 residues","getSTRIDEClassifications");
    else ssClassifications.resize(S->residueSize(),"C");
    return ssClassifications;
  }
  if (data.empty()) {
    if (strict) MstUtils::error("There is no STRIDE data for the requested structure","getSTRIDEClassifications");
    else ssClassifications.resize(S->residueSize(),"X"); // a placeholder for unknown types
    return ssClassifications;
  }
  
  int data_idx = 0; //for indexing data;
  for (Residue* R : S->getResidues()) {
    if (R->getNum() == data[data_idx].resNum) {
      // STRIDE classification found
      ssClassifications.push_back(data[data_idx].ssType);
      data_idx++;
    } else {
      if (strict) MstUtils::error("Residue "+R->getChainID()+MstUtils::toString(R->getNum())+" missing from STRIDE data","getSTRIDEClassifications");
      // interpolate the missing label by looking at the flanking residues with labels
      string NtermResLabel = "C";
      string CtermResLabel = "C";
      if (data_idx > 0) NtermResLabel = data[data_idx - 1].ssType;
      if (data_idx < data.size() - 1) CtermResLabel = data[data_idx + 1].ssType;
      if (NtermResLabel == CtermResLabel) ssClassifications.push_back(NtermResLabel);
      else ssClassifications.push_back("C");
    }
  }
  return ssClassifications;
}
