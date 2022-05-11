#include "msttypes.h"
#include "mstexternal.h"
#include "mstsystem.h"

int main(int argc, char *argv[]) {
  string strideBin = "/scratch/users/swans/STRIDE_repo/stride";
//  Structure P("/scratch/users/swans/MST_workspace/MST/testfiles/1DC7.pdb");
//  Structure P("/scratch/users/swans/databases/singlechain_sim30/PDB/mr/4MR0_A.pdb");
  Structure P("/scratch/users/swans/databases/singlechain_sim30/PDB/zy/4ZYS_A.pdb");

  strideInterface strideObj(strideBin,&P);
  strideObj.computeSTRIDEClassifications();
  vector<string> strideSSType = strideObj.getSTRIDEClassifications();
  
  cout << "Num residues in structure: " << P.residueSize() << ", Num classified: " << strideSSType.size() << endl;
  for (string ssType : strideSSType) {
    cout << ssType << endl;
  }
  
  // test
  vector<string> test;
  test.resize(10,"X");
  test.resize(0,"X");
  return 0;
}
