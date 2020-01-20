#include "msttypes.h"
#include "msttermanal.h"
#include "mstfasst.h"
#include "mstoptions.h"
#include <chrono>

using namespace std;
using namespace MST;

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Tests MST::TERMANAL. Options:");
  op.addOption("t", "path to MST test directory.", true);
  op.addOption("b", "path to a binary FASST database.", true);
  op.setOptions(argc, argv);

  MstUtils::setSignalHandlers();

  Structure S(op.getString("t") + "/2ZTA.pdb");
  vector<Residue*> R = S.getResidues(), subregion;
  for (Residue* res : R) if (RotamerLibrary::hasFullBackbone(res)) subregion.push_back(res);
  Structure subS(subregion);
  vector<Residue*> subR = subS.getResidues();

  FASST F;
  F.setSearchType(FASST::searchType::CA);
  F.readDatabase(op.getString("b"));
  TERMANAL T(&F); T.readRotamerLibrary("testfiles/rotlib.bin");
  T.setCompatMode(true);

  auto begin = chrono::high_resolution_clock::now();
  cout << "scoring..." << endl;
  bool verbose = true;
  vector<mstreal> structScores = T.scoreStructure(subS, NULL, verbose);
  for (int i = 0; i < subR.size(); i++) cout << "structure score for " << *(subR[i]) << " = " << structScores[i] << endl;
  auto end = chrono::high_resolution_clock::now();
  cout << "scoring took " << chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms" << endl;
}
