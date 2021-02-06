#include <iostream>
#include <fstream>
#include <string>

#include "msttypes.h"
#include "msttransforms.h"
#include "mstoptions.h"


/*
read pdb file, select chain A (the purpose is to move chain A around a pocket on chain B)
move chain A by 3-rotation parameters and 3-translation parameters and generate 100 different poses of ChainA, 
each different chain A position will generate a new pdb file that contain both chainA and chainB 
run design for each pdb

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
  op.setOptions(argc, argv);

  srand(time(NULL) + (int) getpid());

  Structure I(op.getString("p")),A;//creat an item of the structure class using input structure
  vector<Residue*> movingResidues;
  //get the info of the moving chain, modified from termify --fs
  if (op.isGiven("s")) {
    selector sel(I);//selector is in msttypes.cpp, probably just for selecting atoms/residues...
    movingResidues = sel.selectRes(op.getString("s"));
    cout << "moving selection gave " << movingResidues.size() << " moving residues, locating..." << endl;
    }
  
  //move the selected residues
  for (int j = 0; j < 300; j++) {
    Structure A = I;
    selector sel(A);
	AtomPointerVector apv = sel.select(op.getString("s")); //getAtoms is a vector of atom pointers
    CartesianPoint cent = apv.getGeometricCenter();//get the x, y, z of the centroid of all the atoms in apv
	//move to origin,otherwise rotation looks like translation
    TransformFactory::translate(-cent).apply(apv);
   
	//rotate and translate
    Transform translation = TransformFactory::translate((MstUtils::randUnit() * 20 -10),(MstUtils::randUnit() * 20 - 10),(MstUtils::randUnit() * 20 -10));
    Transform rotationX= TransformFactory::rotateAroundX(MstUtils::randUnit() * 90 -45);
    Transform rotationY = TransformFactory::rotateAroundY(MstUtils::randUnit() * 90 -45);
    Transform rotationZ = TransformFactory::rotateAroundZ(MstUtils::randUnit() * 90 -45);
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