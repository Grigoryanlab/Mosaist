#include <iostream>
#include <fstream>
#include <string>

#include "msttypes.h"
#include "msttransforms.h"
#include "mstoptions.h"


/*
read pdb file, select chain A (the purpose is to move chain A around a pocket on chain B)
move chain A by 3-rotation parameters and 3-translation parameters and generate n different poses of ChainA,
each different chain A position will generate a new pdb file that contain both chainA and chainB

*/

using namespace std;
using namespace MST;

int main(int argc, char** argv)
{
  MstOptions op;
  op.setTitle("Starting from some arbitrary position close to the PPI pocket, sample positions around this initial position. Options:");
  op.addOption("p", "input PDB file.", true);
  op.addOption("s", "select the moving chain/residues");
  op.addOption("o", "output base name.", true);
  op.addOption("n", "number of output structures");
  op.addOption("t", "set translation distances");
  op.addOption("r", "set rotation degrees");
  op.setOptions(argc, argv);

  srand(time(NULL) + (int) getpid());

  Structure I(op.getString("p")),A;
  vector<Residue*> movingResidues;
  //get the info of the moving chain, modified from termify --fs
  if (op.isGiven("s")) {
    selector sel(I);//selector is in msttypes.cpp, for selecting atoms/residues...
    movingResidues = sel.selectRes(op.getString("s"));
    cout << "moving selection gave " << movingResidues.size() << " moving residues, locating..." << endl;
    }

  int numOfStructures = op.getInt("n",10);
  int distance = op.getInt("t", 10);
  int degree = op.getInt("r", 10);

  //move the selected residues
  for (int j = 0; j < numOfStructures; j++) {
    Structure A = I;
    selector sel(A);
    AtomPointerVector apv = sel.select(op.getString("s")); //getAtoms is a vector of atom pointers
    CartesianPoint cent = apv.getGeometricCenter();//get the x, y, z of the centroid of all the atoms in apv
    //move to origin,otherwise rotation looks like translation
    TransformFactory::translate(-cent).apply(apv);

    //rotate and translate
    Transform translation = TransformFactory::translate((MstUtils::randUnit() *2*distance-distance),(MstUtils::randUnit() *2*distance-distance),(MstUtils::randUnit() *2*distance-distance));
    Transform rotationX= TransformFactory::rotateAroundX(MstUtils::randUnit() *2*degree-degree);
    Transform rotationY = TransformFactory::rotateAroundY(MstUtils::randUnit() *2*degree-degree);
    Transform rotationZ = TransformFactory::rotateAroundZ(MstUtils::randUnit() *2*degree-degree);
    Transform transform = rotationX * rotationY * rotationZ * translation;
	for (int i = 0; i < apv.size(); i++) {
        Atom *atom = apv[i];
        transform.apply(*atom);
    }

    //put it back
    TransformFactory::translate(cent).apply(apv);

    //output the changed pose
    string s = to_string(j);
    A.writePDB(op.getString("o") + s +".pdb");
  }


}
