#include "msttypes.h"
#include "msttransforms.h"

using namespace MST;

vector<vector<real> > covarianceMatrix(AtomPointerVector& A, bool normalize = false) {
  vector<vector<real> > C(3, vector<real>(3, 0));
  CartesianPoint c = A.getGeometricCenter();
  for (int i = 0; i < 3; i++) {
    for (int j = i; j < 3; j++) {
      C[i][j] = 0;
      for (int k = 0; k < A.size(); k++) {
        C[i][j] += ((*A[k])[i] - c[i]) * ((*A[k])[j] - c[j]);
      }
      if (normalize) C[i][j] /= A.size();
      C[j][i] = C[i][j];
    }
  }
  return C;
}

int main(int argc, char** argv) {

  if (argc < 3) {
    MstUtils::error("Usage: ./test [PDB file] [output file base]", "main");
  }
  string pdbFile(argv[1]);
  string outBase(argv[2]);

  // read structure
  Structure S(pdbFile);

  // do some manipulations to it
  Transform T = TransformFactory::rotateAroundX(90);
  T.apply(S);
  S.writePDB(outBase + ".X90.pdb");
  T.apply(S);
  S.writePDB(outBase + ".X180.pdb");
  T.apply(S);
  S.writePDB(outBase + ".X270.pdb");
  T.apply(S);
  S.writePDB(outBase + ".X360.pdb");

  // do a custom transform
  Transform rot2(CartesianPoint(0.492, -0.587, 0.643), CartesianPoint(0.870, 0.310, -0.383), CartesianPoint(0.025, 0.748, 0.663), Transform::byRow);
  Transform tr = TransformFactory::translate(0.0, 0.0, 0.0);
  Transform total = tr * rot2;
  cout << "Rotation matrix:\n" << rot2 << "\ntranslation matrix:\n" << tr << "\nand total:\n" << total;
  Structure So(pdbFile);
  Structure Sr = So;
  total.apply(Sr);
  cout << "resulting RMSD: " << RMSDCalculator::rmsd(So, Sr) << endl;
  Sr.writePDB(outBase + ".cust.pdb");

  AtomPointerVector all = So.getAtoms();
  vector<vector<real> > C = covarianceMatrix(all, true);
  cout << "covariance matrix:\n" << endl;
  for (int i = 0; i < C.size(); i++) {
    for (int j = 0; j < C[i].size(); j++) {
      cout << C[i][j] << " ";
    }
    cout << endl;
  }
  
  printf("TEST DONE\n");
}
