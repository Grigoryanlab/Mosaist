#include "msttypes.h"
#include "msttransforms.h"
#include "mstlinalg.h"

using namespace MST;

int main(int argc, char** argv) {

  if (argc < 3) {
    MstUtils::error("Usage: ./test [PDB file] [output file base]", "main");
  }
  string pdbFile(argv[1]);
  string outBase(argv[2]);

  // test construction
  vector<vector<real> > rot(3, vector<real>(3, 0));
  rot[0][0] = 1;
  rot[1][1] = 1;
  rot[2][2] = 1;
  vector<real> trans(3, 1.5);
  Transform TT(rot, trans);
  cout << "constructed:\n" << TT << endl;

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
  Transform tr = TransformFactory::translate(10.0, 20.0, 30.0);
  Transform total = tr * rot2;
  cout << "Rotation matrix:\n" << rot2 << "\ntranslation matrix:\n" << tr << "\nand total:\n" << total << endl;
  Structure So(pdbFile);
  AtomPointerVector all = So.getAtoms();
  all.center();
  Structure Sr = So;
  total.apply(Sr);
  cout << "resulting RMSD: " << RMSDCalculator::rmsd(So, Sr) << endl;
  TransformRMSD rigidBodyRMSD(So);
  double fastRMSD = rigidBodyRMSD.getRMSD(total);
  cout << "computing the same via the fast method in TransformRMSD: " << fastRMSD << endl;
  Sr.writePDB(outBase + ".cust.pdb");

  // dome some simple matrix algebra
  srand(time(NULL));
  Matrix M(4, 4);
  for (int i = 0; i < M.size(1); i++) {
    for (int j = 0; j < M.size(2); j++) {
      M(i,j) = ((double) rand() / (RAND_MAX));
    }
  }
  cout << "created random matrix:\n" << M << endl;
  Matrix Mi = M.inverse();
  cout << "the inverse of which is:\n" << Mi << endl;
  cout << "product of the matrix by its inverse is:\b" << M*Mi << endl;

  cout << "first row of random matrix:\n" << M.row(0) << endl;
  M.row(0) = Matrix(1, 4, 1.0);;
  cout << "overwrote the first row of random matrix:\n" << M << endl;
  Matrix row = M.row(0);
  cout << "the firs row now is:\n" << row << endl;
  printf("TEST DONE\n");
}
