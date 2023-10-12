// a function that creates an ensemble of randomly docked structures, optionally biased towards certain binding residues for either binding partner, and outputs a distribution of their RMSDs when compared to the correct docked conformation

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <random>
#include <cmath>
#include <iostream>
#include "mstcondeg.h"
#include "mstoptions.h"
#include "mstrotlib.h"
#include "msttypes.h"
#include "msttransforms.h"
#include "mstfasst.h"
#include "mstsequence.h"

using std::cout; using std::endl;
using namespace std;
using namespace MST;

void lilWriteTest(string fileName, string toWrite) {
  fstream myFile;
  myFile.open(fileName, ios::out);
  if( !myFile ) {
     cerr << "Error: file could not be opened" << endl;
     exit(1);
  }
  myFile << toWrite << endl;
  myFile.close();
}

void printMapIntVecInt(map<int,vector<int>> myMap) {
  for(auto it = myMap.begin(); it != myMap.end(); ++it) {
    std::cout << "key: " << it->first << endl;
    vector<int> valuesVect = myMap[it->first];
    for (int i = 0; i < valuesVect.size(); i++) {
      cout << "value: " << valuesVect[i] << endl;
    }
  }
}

AtomPointerVector getBBatoms(AtomPointerVector startingAtoms) {
  AtomPointerVector bbAtoms;
  for (int i = 0; i < startingAtoms.size(); i++) {
    Atom* currAtom = startingAtoms[i];
    string cAtomName = currAtom->getName();
    if ((cAtomName == "CA") or (cAtomName == "N") or (cAtomName == "O") or (cAtomName == "C")) {
      bbAtoms.push_back(currAtom);
    }
  }
  return bbAtoms;
}

AtomPointerVector getHeavyAtoms(AtomPointerVector startingAtoms) {
  AtomPointerVector heavyAtoms;
  for (int i = 0; i < startingAtoms.size(); i++) {
    Atom* currAtom = startingAtoms[i];
    string cAtomName = currAtom->getName();
    char firstChar = cAtomName[0];
    if ((firstChar != 'H')) {
      heavyAtoms.push_back(currAtom);
    }
  }
  return heavyAtoms;
}

// ***
vector<vector<int>> getABindexes(AtomPointerVector aAtoms,AtomPointerVector bAtoms,mstreal conDist) {
  ProximitySearch psA(aAtoms, 15);
  vector<int> aIdx = {};
  vector<int> bIdx = {};

  for (int i = 0; i < bAtoms.size(); i++) {
    vector <int> abContactPts = psA.getPointsWithin(bAtoms[i], 0.0, conDist, false); // checks if the given atom in B is w/in [conDist] Angstroms from any in A
    if (abContactPts.size() > 0) {
      bIdx.push_back(i);
      for (int ii = 0; ii < abContactPts.size(); ii++) {
        int aPoint = abContactPts[ii];
        if (std::find(aIdx.begin(), aIdx.end(), aPoint) == aIdx.end()) {
          aIdx.push_back(aPoint);
        }
      }
    }
  }
  vector<vector<int>> abIdxVector = {aIdx,bIdx};
  return abIdxVector;
}

// gets all the residue contacts between a & b as a dict, with the B indexes as keys and A indexes as a list of values
tuple<AtomPointerVector,AtomPointerVector,vector<int>> getABcontactDict(AtomPointerVector aAtoms,AtomPointerVector bAtoms,mstreal conDist, ProximitySearch psA) {

  AtomPointerVector aCons;
  AtomPointerVector bCons;
  vector<int> bBindxes;

  for (int i = 0; i < bAtoms.size(); i++) {

    Atom* bAtom = bAtoms[i];
    string conNameB = bAtom->getName();
    if (conNameB != "CA") {
      continue;
    }

    vector <int> abContactPts = psA.getPointsWithin(bAtom, 0.0, conDist, false); // checks if the given atom in B is w/in X Angstroms from any in A

    if (abContactPts.size() > 0) {
      for (int ii = 0; ii < abContactPts.size(); ii++) {

        Atom* aAtom = aAtoms[abContactPts[ii]];
        string conNameA = aAtom->getName();
        if (conNameA != "CA") {
          continue;
        }
        //int currResB = bAtom->getResidue()->getResidueIndex();
        //int currResA = aAtom->getResidue()->getResidueIndex();

        if (find(aCons.begin(), aCons.end(), aAtom) == aCons.end()) {
          aCons.push_back(aAtom);
        }
        if (find(bCons.begin(), bCons.end(), bAtom) == bCons.end()) {
          bCons.push_back(bAtom);
        }
        if (find(bBindxes.begin(), bBindxes.end(), i) == bBindxes.end()) {
          bBindxes.push_back(i);
        }
        

      }
    }
  }
  tuple<AtomPointerVector,AtomPointerVector,vector<int>> abConTuple = make_tuple(aCons,bCons,bBindxes);
  return abConTuple;
}

int main(int argc, char** argv) {
  /*time_t start_time = time(NULL);
  cout << "start time: " << ctime(&start_time) << endl;*/
  // Setup and get the input arguments
  MstOptions op;
  op.setTitle("Generates randomly docked pairs of structures, with optional binding residue constraints, and computes their RMSDs and rankings as well as saving one structure and the transformation matrixes for each dock to recapitulate the dock from the structure.");
  op.addOption("m", "use model binding partners for random docking instead of crystal binding partners", false);
  op.addOption("ca", "one half of the correctly docked structure, in pdb form, after having its backbone paired with that of the simulated docked structure via match_backbones.py; if using antibody mode, this should be the antibody ", true);
  op.addOption("da", "one half of the simulated docked structure, in pdb form, after having its backbone paired with that of the correctly docked structure via match_backbones.py (see -ca); if using antibody mode, this should be the antibody ", true);
  op.addOption("cb", "the other half of the correctly docked structure, in pdb form, after having its backbone paired with that of the simulated docked structure via match_backbones.py", true);
  op.addOption("db", "the other half of the simulated docked structure, in pdb form, after having its backbone paired with that of the correctly docked structure via match_backbones.py (see -cb); if using antibody mode, this should be the antibody ", true);
  op.addOption("n", "the number of valid docking positions to generate; defaults to 1 million, which should take around 1 hour per every 500 residues in the complex");
  op.addOption("i", "the number of valid interaction residues required to accept the structure; defaults to 1");
  op.addOption("cla", "the number of clashes allowed in an accepted structure; defaults to 3");
  op.addOption("sd", "the standard deviation of the normal distribution used to pull the docking partners apart from each other; defaults to 1.0; in angstroms");
  op.addOption("al", "an optional list of binding residues for the first docking partner (see -ca); this will bias the random docking towards conformations including those residues in the binding site. Should be a list of tuples, where each tuple has the chain followed by the residue number followed by the residue insertion code (or ' ' if no insertion code). Separate each member of the tuple with a comma, and each tuple with a semi-colon, like 'A,100, ;A,100,A;A,100,B'. If you're only giving binding residues for one of the two partners, it must be this one (i.e. you cannot give bl without giving al)");
  op.addOption("bl", "an optional list of binding residues for the second docking partner (see -cb); this will bias the random docking towards conformations including those residues in the binding site. Should be a list of tuples, where each tuple has the chain followed by the residue number followed by the residue insertion code (or ' ' if no insertion code). Separate each member of the tuple with a comma, and each tuple with a semi-colon, like A,100,;A,100,A;A,100,B.");
  op.addOption("nc", "use native contacts for a & b");
  op.addOption("nnc", "use native contacts for a & b, but disallow any docks with them");
  op.addOption("abm", "an optional antibody / TCR mode, which will treat all the loops of the antibody as the binding residues; the antibody structure must use IMGT numbering and its structure must be entered using the -ca argument not the -cb argument.");
  op.addOption("q", "an optional quick mode for the --al or --abm flags, wherein the docking distribution is skewed towards conformations involving the binding residues given on the A side, to make the calculation faster");
  op.addOption("o", "the output file name base: will be used to save the distribution of LRDPs and rotation matrixes to the randomly rotated reference structure, as a .csv, then the first randomly rotated binding partner, as [base]_A.pdb, and the second randomly rotated binding partner as [base]_B.pdb", true);
  op.addOption("t", "create testing files; provide the output directory to save the testing files. This also sets the number of dockings to just 10 so you aren't inundated with files on accident :)");
  op.setOptions(argc, argv);
  MstUtils::seedRandEngine();

  if (!op.isGiven("al") && !op.isGiven("abm") && op.isGiven("bl")) {
    cout << "you cannot give bl without giving al or abm";
    exit(1);
  }

  int numberDockingsRequired;
  if (op.isGiven("t")) {
    numberDockingsRequired = 10;
  }
  else {
    if (op.isGiven("n")) {
      numberDockingsRequired = stoi(op.getString("n"));
    }
    else {
      numberDockingsRequired = 1000000;
    }
  }

  if (op.isGiven("q") && (!(op.isGiven("al")) && !(op.isGiven("abm")) && !(op.isGiven("nc")))) {
    cout << "you cannot give q without giving al or abm";
    exit(1);
  }

  int contactsRequired;
  if (op.isGiven("i")) {
    contactsRequired = stoi(op.getString("i"));
  }
  else {
    contactsRequired = 1;
  }

  int clashesAllowed;
  if (op.isGiven("cla")) {
    clashesAllowed = stoi(op.getString("cla"));
  }
  else {
    clashesAllowed = 3;
  }

  mstreal normalDistBase = 1.0;

  // load backbone or full-atom structures

  Structure CA0(op.getString("ca"));
  Structure CB0(op.getString("cb"));
  Structure DA(op.getString("da"));
  Structure DB(op.getString("db"));

  // make a full structure for C0, which is connected to CA0 / CB0 rather than a copy of them
  Structure C0;
  for (int i = 0; i < CA0.chainSize(); i++) {
    C0.appendChain(&CA0.getChain(i));
  }
  for (int i = 0; i < CB0.chainSize(); i++) {
    C0.appendChain(&CB0.getChain(i));
  }

  Structure CC(C0); // a copy for doing RMSD calculations with
  AtomPointerVector realAtoms = getHeavyAtoms(CC.getAtoms());
  Structure CAC(CA0);
  AtomPointerVector realAtomsCA = getHeavyAtoms(CAC.getAtoms());
  Structure CBC(CB0);
  AtomPointerVector realAtomsCB = getHeavyAtoms(CBC.getAtoms());

  ProximitySearch psA(realAtomsCA, 15);



  // make a full structure for D, which is connected to DA / DB rather than a copy of them
  Structure D;
  for (int i = 0; i < DA.chainSize(); i++) {
    D.appendChain(&DA.getChain(i));
  }
  for (int i = 0; i < DB.chainSize(); i++) {
    D.appendChain(&DB.getChain(i));
  }

  // if you're using model-structure partners for docking, make CA & CB out of them; otherwise stick with the crystal partners (CA0 & CB0) which is default

  Structure CA;
  Structure CB;
  Structure CenteredRandomA;
  Structure CenteredRandomB;



  if (op.isGiven("m")) {
    CA = Structure(DA); // this duplicates / does not connect CA & DA
    CB = Structure(DB);
  }
  else {
    CA = Structure(CA0);
    CB = Structure(CB0);
  }

  // make a full structure for C, which is connected to CA / CB rather than a copy of them
  Structure C;
  for (int i = 0; i < CA.chainSize(); i++) {
    C.appendChain(&CA.getChain(i));
  }
  for (int i = 0; i < CB.chainSize(); i++) {
    C.appendChain(&CB.getChain(i));
  }

  vector <Residue*> CAreses = CA.getResidues();
  AtomPointerVector CAatoms = getHeavyAtoms(CA.getAtoms());

  vector <Residue*> CBreses = CB.getResidues();
  AtomPointerVector CBatoms = getHeavyAtoms(CB.getAtoms());

  vector <Residue*> DAreses = DA.getResidues();
  AtomPointerVector DAatoms = getHeavyAtoms(DA.getAtoms());

  vector <Residue*> DBreses = DB.getResidues();
  AtomPointerVector DBatoms = getHeavyAtoms(DB.getAtoms());

  tuple<AtomPointerVector,AtomPointerVector,vector<int>> abResCons = getABcontactDict(CAatoms,CBatoms,8.0,psA);


  RMSDCalculator rc;
  vector <tuple<mstreal,mstreal,mstreal,mstreal,mstreal,mstreal,mstreal,mstreal,mstreal,mstreal,mstreal,mstreal,mstreal>> rmsdList {};


  // position 1
  if (op.isGiven("t")) {
    C.writePDB(op.getString("t") + "position1" + ".pdb");
  }

  // translate each to the origin

  Transform TCA0 = TransformFactory::translate(-CAatoms.getGeometricCenter());
  TCA0.apply(CA);

  Transform TCB0 = TransformFactory::translate(-CBatoms.getGeometricCenter());
  TCB0.apply(CB);

  // randomly rotate each

  mstreal xarAngle = MstUtils::randUnit(0,360);
  Transform TXAR = TransformFactory::rotateAroundX(xarAngle);
  TXAR.apply(CA);
  mstreal xbrAngle = MstUtils::randUnit(0,360);
  Transform TXBR = TransformFactory::rotateAroundX(xbrAngle);
  TXBR.apply(CB);

  mstreal yarAngle = MstUtils::randUnit(0,360);
  Transform TYAR = TransformFactory::rotateAroundY(yarAngle);
  TYAR.apply(CA);
  mstreal ybrAngle = MstUtils::randUnit(0,360);
  Transform TYBR = TransformFactory::rotateAroundY(ybrAngle);
  TYBR.apply(CB);

  mstreal zarAngle = MstUtils::randUnit(0,360);
  Transform TZAR = TransformFactory::rotateAroundZ(zarAngle);
  TZAR.apply(CA);
  mstreal zbrAngle = MstUtils::randUnit(0,360);
  Transform TZBR = TransformFactory::rotateAroundZ(zbrAngle);
  TZBR.apply(CB);

  CenteredRandomB = Structure(CB);
  CenteredRandomA = Structure(CA);
  AtomPointerVector CRAatoms = CenteredRandomA.getAtoms();
  for (int a = 0; a < CRAatoms.size(); a++) {
    Atom* currAtom = CRAatoms[a];
    currAtom->clearAlternatives();
    currAtom->addAlternative(currAtom);
  }

  Structure CenteredRandom;
  for (int i = 0; i < CenteredRandomA.chainSize(); i++) {
    CenteredRandom.appendChain(&CenteredRandomA.getChain(i));
  }
  for (int i = 0; i < CenteredRandomB.chainSize(); i++) {
    CenteredRandom.appendChain(&CenteredRandomB.getChain(i));
  }

  CenteredRandomA.writePDB(op.getString("o") + "_A.pdb");
  CenteredRandomB.writePDB(op.getString("o") + "_B.pdb");

  // position 2
  if (op.isGiven("t")) {
    C.writePDB(op.getString("t") + "position2" + ".pdb");
  }

  // if binding residues are given for A (or A has binding residues by virtue of this being run in antibody-antigen binding mode) get them; same for B

  AtomPointerVector CAbinderAtoms;
  AtomPointerVector CBbinderAtoms;

  vector<tuple <string, int, char>> aResesDetails; // stores details for binding residues for A
  vector<tuple <string, int, char>> bResesDetails; // stores details for binding residues for B

  Structure aBinderS;
  AtomPointerVector CAdisallowed;
  AtomPointerVector CBdisallowed;
  vector<int> CBindexDisallowed;

  if (op.isGiven("nc")) {
    CAbinderAtoms = get<0>(abResCons);
  }

  else if (op.isGiven("nnc")) {
    CAdisallowed = get<0>(abResCons);
    CBdisallowed = get<1>(abResCons);
    CBindexDisallowed = get<2>(abResCons);
  }

  else if (op.isGiven("abm")) {
    for (int i = 0; i < CAreses.size(); i++) {
      Residue* currR = CAreses[i];
      int resNum = currR->getNum();
      bool isLoop = false;

      //check if loop by normal numbering, unless mm mode, in which case check if loop by crystal numbering

      if ((27 <= resNum && resNum <= 38) || (56 <= resNum && resNum <= 65) || (105 <= resNum && resNum <= 117))
      {
        isLoop = true;
      }

      if (isLoop == true) {
        aBinderS.addResidue(currR);
        AtomPointerVector cResAtoms = currR->getAtoms();
        for (int ca = 0; ca < cResAtoms.size(); ca++) {
          Atom* currA = cResAtoms[ca];
          string cResAtomName = currA->getName();
          if (!(op.isGiven("dq")) && (cResAtomName != "CA")) { // if not DOCKQ mode, and not a CA atom...
            continue;
          } // only get indexes for CAs, as contacts will be defined as inter-CA-distance of 10 angstroms or less, unless using DOCKQ mode in which case it's all atoms 5 angstroms or less
          CAbinderAtoms.push_back(currA);
        }
      }
    }
  }

  else if (op.isGiven("al")) {
    // get the details for A's binding residues
    string aBinders = op.getString("al");
    vector<string> aResSplit1 = MstUtils::split(aBinders, ";");
    for (int i = 0; i < aResSplit1.size(); i++) {
        vector<string> aResSplit2 = MstUtils::split(aResSplit1[i], ",");
        string resChain = aResSplit2[0];
        int resNum = std::stoi(aResSplit2[1]);
        string resIcodeBase = aResSplit2[2];
        char resIcode = resIcodeBase[0];
        auto aResDetails = std::make_tuple (resChain, resNum, resIcode);
        aResesDetails.push_back(aResDetails);
    }

    // go over each residue, and if its chain / number / ID code match, add it to the binding residues list! Also make a list of the binder indexes

    for (int i = 0; i < CAreses.size(); i++) {
      Residue* currR = CAreses[i];
      string resChain = currR->getChainID();
      int resNum = currR->getNum();
      char resIcode = currR->getIcode();
      for (int ii = 0; ii < aResesDetails.size(); ii++) {
        auto aResTuple = aResesDetails[ii];
        bool isbindingRes = false;

        //check if binding residue by normal numbering, unless mm mode, in which case check if binding residue by crystal numbering

        if ((get<0>(aResTuple) == resChain) && (get<1>(aResTuple) == resNum) && (get<2>(aResTuple) == resIcode)) {
          isbindingRes = true;
        }
        //}

        if (isbindingRes) {
          aBinderS.addResidue(currR);
          AtomPointerVector cResAtoms = currR->getAtoms();
          for (int ca = 0; ca < cResAtoms.size(); ca++) {
            Atom* currA = cResAtoms[ca];
            string cResAtomName = currA->getName();
            if (!(op.isGiven("dq")) && (cResAtomName != "CA")) { // if not DOCKQ mode, and not a CA atom...
              continue;
            }
            CAbinderAtoms.push_back(currA);
          }
        }
      }
    }
  }

  if (op.isGiven("q")) {

      if (op.isGiven("t")) {
        C.writePDB(op.getString("t") + "positionB43" + ".pdb");
        Structure CAbinderS(CAbinderAtoms);
        CAbinderS.writePDB(op.getString("t") + "CAbindersB43" + ".pdb");
      }

      

      CartesianPoint geoCenterA = CAbinderAtoms.getGeometricCenter();
      Transform TZ = TransformFactory::alignVectorWithXAxis(geoCenterA);
      TZ.apply(CA);

      if (op.isGiven("t")) {
        C.writePDB(op.getString("t") + "position3" + ".pdb");
      }
    }

  // if binding residues are given for B, pre-store the details to iterate over easily

  if (op.isGiven("nc")) {
    CBbinderAtoms = get<1>(abResCons);
  }

  if (op.isGiven("bl")) {
    string bBinders = op.getString("bl");
    vector<string> bResSplit1 = MstUtils::split(bBinders, ";");
    for (int i = 0; i < bResSplit1.size(); i++) {
        vector<string> bResSplit2 = MstUtils::split(bResSplit1[i], ",");
        string resChain = bResSplit2[0];
        int resNum = std::stoi(bResSplit2[1]);
        char resIcode = bResSplit2[2][0];
        auto bResDetails = std::make_tuple (resChain, resNum, resIcode);
        bResesDetails.push_back(bResDetails);
    }

    for (int i = 0; i < CBreses.size(); i++) {
      Residue* currR = CBreses[i];
      string resChain = currR->getChainID();
      int resNum = currR->getNum();
      char resIcode = currR->getIcode();
      for (int ii = 0; ii < bResesDetails.size(); ii++) {
        auto bResTuple = bResesDetails[ii];

        bool isbindingRes = false;

        //check if binding residue by normal numbering, unless mm mode, in which case check if binding residue by crystal numbering

        if ((get<0>(bResTuple) == resChain) && (get<1>(bResTuple) == resNum) && (get<2>(bResTuple) == resIcode)) {
          isbindingRes = true;
        }
        //}

        if (isbindingRes) {
          AtomPointerVector cResAtoms = currR->getAtoms();
          for (int ca = 0; ca < cResAtoms.size(); ca++) {
            string cResAtomName = cResAtoms[ca]->getName();
            if (cResAtomName == "CA") {
              int cResIndx = cResAtoms[ca]->getIndex();
              CBbinderAtoms.push_back(cResAtoms[ca]);
            }
          }
        }
      }
    }
  }

  /*vector<int> bBinderIndexes; // store the atom indexes of CA binding atoms
  for (int i = 0; i < CBatoms.size(); i++) { // go over all atoms in CB...
    Atom* currAtom = CBatoms[i];
    string currAtomName = currAtom->getName();
    if (currAtomName != "CA") { // if not DOCKQ mode, and not a CA atom...
      continue;
    }
    Residue* currRes = currAtom->getParent();
    string currChain = currRes->getChainID();
    int currResNum = currRes->getNum();
    char currResIcode = currRes->getIcode();
    for (int ii = 0; ii < bResesDetails.size(); ii++) {
      auto bResTuple = bResesDetails[ii];
      if ((get<0>(bResTuple) == currChain) && (get<1>(bResTuple) == currResNum) && (get<2>(bResTuple) == currResIcode)) { // if the residue is a binding residue...
        bBinderIndexes.push_back(i);
      }
    }
  }*/

  // let's get started on the random docking now! Set up variables (d = number of random dockings required, t = a variable used to count the initial 10 structures to be saved for testing mode)

  int d = 0;
  int t = 0;
  int t2 = 0;
  //int t3 = 0;
  int fails = 0;

  // set up important vectors...

  // keep a list of indexes of previous clashes
  vector <int> recentClashes;
  // and have a list of current clashes
  vector <int> currClashes;
  // as well as a list of indexes possible to be involved in contacts (CA atoms)
  vector <int> bPossibleContacts;

  // set the curr coords of CA to be alt cords, so that they can be reset every new random docking

  for (int a = 0; a < CAatoms.size(); a++) {
    Atom* currAtom = CAatoms[a];
    currAtom->clearAlternatives();
    currAtom->addAlternative(currAtom);
  }

  // set the current coords of CB to be alternate coords, so that they can be reset every new random docking - also fill out which indexes are possible to make contacts in B, and get their atoms

  for (int a = 0; a < CBatoms.size(); a++) {
    Atom* currAtom = CBatoms[a];
    currAtom->clearAlternatives();
    currAtom->addAlternative(currAtom);
    // as well as get indexes for the potential contact residues, if not already gotten from abm / al / bl flags being given
    /*if (op.isGiven("bl")) {
      continue;
    }*/

    if (op.isGiven("nc") || op.isGiven("bl")) {
      if (std::find(CBbinderAtoms.begin(), CBbinderAtoms.end(), currAtom) != CBbinderAtoms.end()) {
        bPossibleContacts.push_back(a);
      }
    }

    else {
      string currName = currAtom->getName();
      if (currName == "CA") {
        bPossibleContacts.push_back(a);
        CBbinderAtoms.push_back(currAtom);
      }
    }
  }

  // and fill out which ones are possible to make contacts in A / get their atoms

  for (int a = 0; a < CAatoms.size(); a++) {
    Atom* currAtom = CAatoms[a];
    if (op.isGiven("al") || op.isGiven("abm") || op.isGiven("nc")) {
      continue;
    }
    else {
      string currName = currAtom->getName();
      if (currName != "CA") { // if not DOCKQ mode, and not a CA atom...
        continue;
      }
      CAbinderAtoms.push_back(currAtom);
    }
  }

  /*if (op.isGiven("bl")) {
    bPossibleContacts = bBinderIndexes;
  }*/

  // set up xyz high & low values to be used later in the loop, for checking xLow & xHigh when pulling out proteins

  mstreal xLow;
  mstreal xHigh;
  mstreal yLow;
  mstreal yHigh;
  mstreal zLow;
  mstreal zHigh;

  // and variables for calculating the Y & Z extents for the Y & Z translations

  mstreal ayExtent;
  mstreal azExtent;
  mstreal axLow;
  mstreal axHigh;
  mstreal ayLow;
  mstreal ayHigh;
  mstreal azLow;
  mstreal azHigh;
  mstreal ayLen;
  mstreal azLen;

  mstreal byExtent;
  mstreal bzExtent;
  mstreal bxLow;
  mstreal bxHigh;
  mstreal byLow;
  mstreal byHigh;
  mstreal bzLow;
  mstreal bzHigh;
  mstreal byLen;
  mstreal bzLen;

  //& docking requirements

  mstreal clashDistance = 3.0;
  mstreal contactDistance = 8.0;

  if (op.isGiven("nc")) {
    cout << "# of partner B's binding CAs:" << endl;
    cout << CBbinderAtoms.size() << endl;
    cout << "# of partner A's binding CAs:" << endl;
    cout << CAbinderAtoms.size() << endl;
  }
  
  

  while (d < numberDockingsRequired) { // while the number of docked structures is still insufficient...

    // reset to alt coors

    for (int i = 0; i < CAreses.size(); i++) {
      CAreses[i]->makeAlternativeMain(0); // reset to origin position
    }

    for (int i = 0; i < CBreses.size(); i++) {
      CBreses[i]->makeAlternativeMain(0); // reset to origin position
    }

    //position 4

    if ((op.isGiven("t")) && (t == 0)) {
      C.writePDB(op.getString("t") + "position4" + ".pdb");
    }

    // rotate CA randomly along the x, y, and z axes, if no binding residues given for A, or limited to prevent A's binding residues from facing too far away from B if binding residues specified

    if (!op.isGiven("q")) {
      mstreal xaAngle = MstUtils::randUnit(0,360);
      Transform TXA = TransformFactory::rotateAroundX(xaAngle);
      TXA.apply(CA);
      mstreal yaAngle = MstUtils::randUnit(0,360);
      Transform TYA = TransformFactory::rotateAroundY(yaAngle);
      TYA.apply(CA);
      mstreal zaAngle = MstUtils::randUnit(0,360);
      Transform TZA = TransformFactory::rotateAroundZ(zaAngle);
      TZA.apply(CA);

      //position 5
      if ((op.isGiven("t")) && (t == 0)) {
        C.writePDB(op.getString("t") + "position5" + ".pdb");
      }
    }

    // rotate CB randomly along the x, y, and z axes

    mstreal xbAngle = MstUtils::randUnit(0,360);
    Transform TXB = TransformFactory::rotateAroundX(xbAngle);
    TXB.apply(CB);
    mstreal ybAngle = MstUtils::randUnit(0,360);
    Transform TYB = TransformFactory::rotateAroundY(ybAngle);
    TYB.apply(CB);
    mstreal zbAngle = MstUtils::randUnit(0,360);
    Transform TZB = TransformFactory::rotateAroundZ(zbAngle);
    TZB.apply(CB);

    //position 6
    if ((op.isGiven("t")) && (t == 0)) {
      C.writePDB(op.getString("t") + "position6" + ".pdb");
    }

    // make standard deviations for random translations in the Y & Z directions, based on the Y & Z lengths of CB - if binding residues are given for a

    if ((op.isGiven("al") || op.isGiven("abm") || op.isGiven("nc")) && op.isGiven("q")) {

      ProximitySearch::calculateExtent(CAbinderAtoms,axLow,ayLow,azLow,axHigh,ayHigh,azHigh);
      ayExtent = ayHigh - ayLow;
      azExtent = azHigh - azLow;

      ProximitySearch::calculateExtent(CBatoms,bxLow,byLow,bzLow,bxHigh,byHigh,bzHigh);
      byExtent = byHigh - byLow;
      bzExtent = bzHigh - bzLow;
      
      mstreal yRand = MstUtils::randUnit(0,(1.0 / 2.0) * (ayExtent + byExtent));
      mstreal zRand = MstUtils::randUnit(0,(1.0 / 2.0) * (azExtent + bzExtent));


      // use those to randomly translate CB in the Y & Z axes so that the angle of binding can vary, with a normal distribution with standard deviation of 1/6 the distance between the furthest Y points / Z points (so 1/3 the distance from the middle and furthest Y / Z point, which means 0.13% will be beyond the range of being able to connect, which is fine)

      for (int a = 0; a < CBatoms.size(); a++) {
        mstreal yStart = CBatoms[a]->getY();
        CBatoms[a]->setY(yStart + yRand);
        mstreal zStart = CBatoms[a]->getZ();
        CBatoms[a]->setZ(zStart + zRand);
      }

      //position 7
      if ((op.isGiven("t")) && (t == 0)) {
        C.writePDB(op.getString("t") + "position7" + ".pdb");
      }

      // if there are B binder atoms, check if they're outside the Y or Z ranges of the A binder atoms, and so could never meet; go back to start without attempting to dock if so

      if (op.isGiven("bl") || op.isGiven("nc")) {

        ProximitySearch::calculateExtent(CBbinderAtoms,bxLow,byLow,bzLow,bxHigh,byHigh,bzHigh);
        
        if ((byLow > ayHigh) || (byHigh < ayLow) || (bzLow > azHigh) || (bzHigh < azLow)) {
          continue;
        }

        // the binding residues for A point positive; if those for B also point positive then pulling it in the + direction will never cause them to meet

        if (CBbinderAtoms.getGeometricCenter().getX() > 0) {
          continue;
        }

        // speed up by making a single big pull of the partners apart, equal to 3/4 the smaller x extent

        mstreal xJump = 0.75 * min(axHigh - axLow, bxHigh - bxLow);
        for (int a = 0; a < CBatoms.size(); a++) {
          mstreal xStart = CBatoms[a]->getX();
          CBatoms[a]->setX(xStart + xJump);
        }

        if ((op.isGiven("t")) && (t == 0)) {
          C.writePDB(op.getString("t") + "position7.5" + ".pdb");
        }
      }
    }

    // set up proximity search with A, since B will be the partner moving (if you move an object you need to re-make a prox-search or it won't work / will segfault)

    ProximitySearch ps;
    ps = ProximitySearch(CAatoms, 15);

    // and a proximity search with just A's CA binding residues - this is a ProximitySearch*

    ProximitySearch ps2(CAbinderAtoms, 15); // for valid contacts
    ProximitySearch ps3;
    if (op.isGiven("nnc")) {
      ps3 = ProximitySearch(CAdisallowed, 15); // for invalid contacts
    }
    

    // while a proper docked stage has not been reached, but there are still opporunitites for success

    bool absoluteSuccess = false;
    bool absoluteFailure = false;

    /// reset recent clashes, if left-over from the last round
    recentClashes.resize(0);

    mstreal contactsCount;
    mstreal xRand;

    while ((absoluteSuccess == false) && (absoluteFailure == false)) {

      // reset current clashes to be empty
      currClashes.resize(0);

      // pull B atoms out along the X axis, in small steps with a normal distrubtion of ~1.0 angstrom and mew of 0, except always positive so it's pulling away from A - can try variations on the 1.0 and adjust

      xRand = fabs(MstUtils::randNormal(0,normalDistBase));

      for (int a = 0; a < CBatoms.size(); a++) {
        mstreal xStart = CBatoms[a]->getX();
        CBatoms[a]->setX(xStart + xRand);
      }

      //position 8
      if (op.isGiven("t") && (t == 0)) {
        C.writePDB(op.getString("t") + "position8" + ".pdb");
      }

      // reject / continue depending on if there are too many clashes or not

      int clashCount = 0;
      currClashes.resize(0);

      // first check positions in B that previously had clashes...

      for (int a = 0; a < recentClashes.size(); a++) {
        int currBcheck = recentClashes[a];
        Atom* currBatom = CBatoms[currBcheck];
        CartesianPoint currBcoor = currBatom->getCoor();
        vector <int> currClashPts = ps.getPointsWithin(currBcoor, 0, clashDistance, false);
        int currBsize = currClashPts.size();
        if (currBsize > 0) {
          currClashes.push_back(currBcheck);
          clashCount+=currBsize;
        }
        if (clashCount > clashesAllowed) {
          break;
        }
      }

      if (clashCount > clashesAllowed) {
        recentClashes = currClashes;
        continue;
      }

      // then if you're not over the clash limit, check the rest of the valid residues

      for (int a = 0; a < CBatoms.size(); a++) {
        int currBcheck = a;
        if (std::find(recentClashes.begin(), recentClashes.end(), currBcheck) != recentClashes.end()) {
          continue; // don't double-count if already checked as part of recent clashesf
        }

        Atom* currBatom = CBatoms[currBcheck];
        CartesianPoint currBcoor = currBatom->getCoor();
        vector <int> currClashPts = ps.getPointsWithin(currBcoor, 0, clashDistance, false);
        int currBsize = currClashPts.size();
        if (currBsize > 0) {
          currClashes.push_back(currBcheck);
          clashCount+=currBsize;
        }
        if (clashCount > clashesAllowed) {
          break;
        }
      }

      if (clashCount > clashesAllowed) {
        recentClashes = currClashes;
        continue;
      }

      //position 9
      if (op.isGiven("t") && (t == 0)) {
        C.writePDB(op.getString("t") + "position9" + ".pdb");
        //exit(0);
      }

      // we've only reached this point if there are an aceeptable number of clashes, so accept if there are at least the required number of contacts - possibly limited by antibody binding mode, a binding residues, and/or b binding residues

      if (op.isGiven("nnc")) {
        for (int a = 0; a < CBindexDisallowed.size(); a++) { // go over all CA atoms in B possible to be part of contacts
          int currBcheck = CBindexDisallowed[a];
          //Residue* currBres = CBatoms[currBcheck]->getParent();
          CartesianPoint currBcoor = CBatoms[currBcheck]->getCoor();
          vector <int> currContactPts = ps3.getPointsWithin(currBcoor, 0.0, contactDistance, false); // ps3 is made from the disallowed A atoms
          if (currContactPts.size() > 0) {
            absoluteFailure = true;
            break;
          }
        }
      }

      if (absoluteFailure == true) {
        break;
      }

      contactsCount = 0;
      vector<pair<Residue*,Residue*>> pastContacts = {};

      for (int a = 0; a < bPossibleContacts.size(); a++) { // go over all CA atoms in B possible to be part of contacts
        int currBcheck = bPossibleContacts[a];
        //Residue* currBres = CBatoms[currBcheck]->getParent();
        CartesianPoint currBcoor = CBatoms[currBcheck]->getCoor();
        vector <int> currContactPts = ps2.getPointsWithin(currBcoor, 0.0, contactDistance, false); // ps2 is only checking for CA atoms in A possible to be part of contacts
        contactsCount += currContactPts.size();
        if (contactsCount >= contactsRequired) {
          absoluteSuccess = true;
          break;
        }
      }

      //check if all atoms in CB are further positive than all atoms in CA by over 8 angstroms... b/c if so the docking is irrecoverable

      ProximitySearch::calculateExtent(CB,xLow,yLow,zLow,xHigh,yHigh,zHigh);
      mstreal bX = xLow;
      ProximitySearch::calculateExtent(CA,xLow,yLow,zLow,xHigh,yHigh,zHigh);
      mstreal aX = xHigh;

      if (bX - aX > 8.0) {
        absoluteFailure = true;
      }
    }

    if (absoluteFailure) {
      fails++;
      if (op.isGiven("t") && (t2 == 0)) {
        C.writePDB(op.getString("t") + "position11_failedWithContacts" + to_string(contactsCount) + ".pdb");
        t2++;
      }
      continue;
    }
    else {
      d++; // accepted! :D
      if (op.isGiven("t") && (t2 == 0)) {
        C.writePDB(op.getString("t") + "position11_succeededWithContacts" + to_string(contactsCount) + ".pdb");
        t2++;
      }
    }

    if (d%10000 == 0) {
      cout << d << " total docks accepted..." << endl;
    }

    // for those accepted, compare to the correct structure to calculate the RMSD

    mstreal simulatedRealRMSDorDOCKQ;
    vector <mstreal> translationA;
    //vector <mstreal> translationB;
    vector <vector <mstreal>> rotationsA;
    //vector <vector <mstreal>> rotationsB;

    /*rc.align(CA.getAtoms(),realAtomsCA,C);
    if (op.isGiven("t") && (t == 0)) {
      C.writePDB(op.getString("t") + "Cposition12.pdb");
      CAC.writePDB(op.getString("t") + "CACposition12.pdb");
    }
    mstreal bRMSD = rc.rmsd(CB.getAtoms(),realAtomsCB); */// gets RMSD without alignment

    rc.align(CB.getAtoms(),realAtomsCB,C); // aligns CB with the original location of CB, but transforms the whole structure C

    if (op.isGiven("t") && (d == 10)) {
      C.writePDB(op.getString("t") + "RandDockForComp.pdb");
      CBC.writePDB(op.getString("t") + "CBCposition13.pdb");
    }
    mstreal aRMSD = rc.rmsd(CA.getAtoms(),realAtomsCA); // gets RMSD without alignment

    simulatedRealRMSDorDOCKQ = aRMSD;

    if (op.isGiven("t") && (d == 10)) {
      CenteredRandomB.writePDB(op.getString("t") + "CRB_before.pdb");
    }

    rc.align(CB.getAtoms(),CenteredRandomB.getAtoms(),C); // align the randomly docked structure to the centered random structure by parnter B

    if (op.isGiven("t")) {
      CenteredRandomB.writePDB(op.getString("t") + "CRB_after.pdb");
    }

    if (op.isGiven("t")) {
      C.writePDB(op.getString("t") + "C_B_aligned");
    }

    if (op.isGiven("t")) {
      CenteredRandomA.writePDB(op.getString("t") + "CRA_before.pdb");
    }

    rc.align(CRAatoms,CA.getAtoms(),CenteredRandomA); // align the centered random structure A to the random dock on A, so you can get the translation partner A needs to go to make that random dock

    if (op.isGiven("t")) {
      CenteredRandomA.writePDB(op.getString("t") + "CRA_after.pdb");
    }

    translationA =rc.lastTranslation();
    rotationsA = rc.lastRotation();


    if (op.isGiven("t")) {
      cout << translationA[0] << " " << translationA[1] << " " << translationA[2] << " " << rotationsA[0][0] << endl;
    }

    for (int i = 0; i < CRAatoms.size(); i++) {
      CRAatoms[i]->makeAlternativeMain(0); // reset CRA to origin position
    }

    if (op.isGiven("t") && (d == 10)) {
      Transform RH = Transform(rotationsA,translationA);
      RH.apply(CenteredRandomA);
      CenteredRandomA.writePDB(op.getString("t") + "CRA_rehydrated.pdb");
      cout << "rmsd was:" << endl;
      cout << simulatedRealRMSDorDOCKQ << endl;
      exit(0);
    }

    rmsdList.push_back(make_tuple(simulatedRealRMSDorDOCKQ,translationA[0],translationA[1],translationA[2],rotationsA[0][0],rotationsA[0][1],rotationsA[0][2],rotationsA[1][0],rotationsA[1][1],rotationsA[1][2],rotationsA[2][0],rotationsA[2][1],rotationsA[2][2]));
    
    // if in testing mode, save the pdb files of the first 10 random things accepted, with their scores in their titles

    if (op.isGiven("t")) {
      if (t < 10) {
        C.writePDB(op.getString("t") + "position11_succeededWithContacts" + to_string(contactsCount) + to_string(t) + ".pdb");
        //cout << "angleOff was: " << endl;
        //cout << angleOff << endl;
        t++;
      }
    }
  }

  // save the distribution of RMSDs as just a list of all the RMSDs, in a text file, comma separated - also compare to comparison RMSD (best correct vs model RMSD)

  cout << "# failed: " << fails << " # worked: " << d << endl;


  mstreal comparisonRMSDorDOCKQ;

  if (op.isGiven("t")) {
    mstreal OLDbRMSD = rc.rmsd(DB.getAtoms(),realAtomsCB); // gets RMSD without alignment
    cout << "testing OLDbRMSD: " << OLDbRMSD << endl;
  }
  /*rc.align(DA.getAtoms(),realAtomsCA,D); // aligns DA with the original location of CA, but transforms the whole structure D - if I'm doing this right lol
  if (op.isGiven("t")) {
    DA.writePDB(op.getString("t") + "DAposition14.pdb");
    CAC.writePDB(op.getString("t") + "CACposition14.pdb");
    D.writePDB(op.getString("t") + "Dposition14.pdb");
    CC.writePDB(op.getString("t") + "CCposition14.pdb");
  }
  mstreal bRMSD = rc.rmsd(DB.getAtoms(),realAtomsCB);*/ // gets RMSD without alignment
  rc.align(DB.getAtoms(),realAtomsCB,D); // aligns CB with the original location of CB, but transforms the whole structure C - if I'm doing this right lol
  if (op.isGiven("t")) {
    D.writePDB(op.getString("t") + "Dposition15.pdb");
    CC.writePDB(op.getString("t") + "CCposition15.pdb");
  }
  mstreal aRMSD = rc.rmsd(DA.getAtoms(),realAtomsCA); // gets RMSD without alignment
  comparisonRMSDorDOCKQ = aRMSD;

  if (op.isGiven("t")) {
    cout << "testing aRMSD... " << aRMSD << endl;
    //cout << "testing bRMSD... " << bRMSD << endl;
    cout << "testing comparisonRMSDorDOCKQ... " << comparisonRMSDorDOCKQ << endl;
  }


  mstreal pValueCount = 1;
  mstreal totCount = rmsdList.size();
  cout << "RMSD between correct structure & structure to evaluate was: " << comparisonRMSDorDOCKQ << endl;

  sort(rmsdList.begin(), rmsdList.end());

  string outFileName = op.getString("o");
  std::ofstream outFile(outFileName + ".csv");

  for (int i = 0; i < rmsdList.size(); i++) {
    if (op.isGiven("t")) {
      cout << i << endl;
      cout << to_string(get<0>(rmsdList[i])) << endl;
      cout << to_string(get<1>(rmsdList[i])) << endl;
    }

    outFile << to_string(i); //put the ranking first...
    outFile << "," << to_string(get<0>(rmsdList[i])); //then the RMSD...
    outFile << "," << to_string(get<1>(rmsdList[i])); //then the translation & rotation matrix info
    outFile << "," << to_string(get<2>(rmsdList[i]));
    outFile << "," << to_string(get<3>(rmsdList[i]));
    outFile << "," << to_string(get<4>(rmsdList[i]));
    outFile << "," << to_string(get<5>(rmsdList[i]));
    outFile << "," << to_string(get<6>(rmsdList[i]));
    outFile << "," << to_string(get<7>(rmsdList[i]));
    outFile << "," << to_string(get<8>(rmsdList[i]));
    outFile << "," << to_string(get<9>(rmsdList[i]));
    outFile << "," << to_string(get<10>(rmsdList[i]));
    outFile << "," << to_string(get<11>(rmsdList[i]));
    outFile << "," << to_string(get<12>(rmsdList[i]));
    outFile << "\n";
  }
  outFile.close();

  cout << "testing" << endl;
  return(0);
}

//example command: ./dockingDistribution --ca "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgTERMs/alphafold8080ColabResultsAb_matchedbb/crystal/6yu8_bef6732d-1f2b-41e2-960b-d5c68c4250b2_0.9186929628704293.pdb" --cb "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgTERMs/alphafold8080ColabResultsAg_matchedbb/crystal/6yu8_4acaabcc-e98e-4ddc-a8c1-cece355b44c3_2.555964177325445.pdb" --da "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgTERMs/alphafold8080ColabResultsAb_matchedbb/af/6yu8_bef6732d-1f2b-41e2-960b-d5c68c4250b2_0.9186929628704293.pdb" --db "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgTERMs/alphafold8080ColabResultsAg_matchedbb/af/6yu8_4acaabcc-e98e-4ddc-a8c1-cece355b44c3_2.555964177325445.pdb"  --o "/dartfs/rc/lab/G/Grigoryanlab/home/coy/Dartmouth_PhD_Repo/ColabAf6yu8DD.txt" --n 10000 --abm