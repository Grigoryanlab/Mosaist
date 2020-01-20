#include "msttypes.h"
#include "mstoptions.h"
#include "msttermanal.h"
#include "mstrotlib.h"

int main(int argc, char** argv) {
  // Setup and get the input arguments
  MstOptions op;
  op.setTitle("Scores the input structure using TERMANAL and outputs three structures whose B-factors are annotated with the design, abundance, and structure scores");
  op.addOption("pdb", "input PDB file to score.", true);
  op.addOption("sel", "selection within the input PDB to score.");
  op.addOption("out", "output prefix for the B-factor annotated output structures (prefix.{dsc,abd,ssc}.pdb).", true);
  op.addOption("tab", "output file containing the three sets of scores in tabular format.");
  op.addOption("fdb", "binary FASST database to search for matches in.", true);
  op.addOption("rlib", "rotamer library to compute contact degree with.", true);
  op.addOption("verb", "if true, prints scoring details as they are computed (default=false).");
  op.addOption("cd", "contact degree cutoff when forming TERMs (default=0.1).");
  op.addOption("pad", "number of residues to pad contacting residues with when forming TERMs (default=2).");
  op.addOption("rmsd", "RMSD cutoff when searching for matches (default=2.0).");
  op.addOption("count", "maximum number of matches to search for within the specified RMSD cutoff (default=50).");
  op.addOption("ps", "pseudocount to add to the smoothed scores before computing their logs (default=0.01).");
  op.setOptions(argc, argv);

  // Read in the structure, FASST database, and rotamer library
  Structure S(op.getString("pdb"));
  selector sele(S);
  Structure selS = op.isGiven("sel") ? sele.selectRes(op.getString("sel")) : S;
  vector<Residue*> selR = selS.getResidues();
  vector<Residue*> bbR;
  for (Residue* res : selR) if (RotamerLibrary::hasFullBackbone(res)) bbR.push_back(res);
  Structure bbS(bbR);
  FASST F; F.readDatabase(op.getString("fdb"), 2);
  TERMANAL T(&F); T.readRotamerLibrary(op.getString("rlib"));
  if (op.isGiven("cd")) T.setCDCut(op.getReal("cd"));
  if (op.isGiven("pad")) T.setPad(op.getInt("pad"));
  if (op.isGiven("rmsd")) T.setRMSDCut(op.getReal("rmsd"));
  if (op.isGiven("count")) T.setMatchCount(op.getInt("count"));
  if (op.isGiven("ps")) T.setPseudocount(op.getReal("ps"));

  // Score each residue in the structure
  vector<pair<mstreal, mstreal>> scoreParts;
  vector<mstreal> structScores = T.scoreStructure(bbS, &scoreParts, op.getBool("verb"));

  // Create the B-factor annotated output PDBs using the computed scores
  Structure dscS(bbS), abdS(bbS), sscS(bbS);
  vector<Residue*> dscR = dscS.getResidues(), abdR = abdS.getResidues(), sscR = sscS.getResidues();
  int numRes = sscR.size();
  for (int i = 0; i < numRes; i++) {
    int numAtoms = sscR[i]->atomSize();
    for (int j = 0; j < numAtoms; j++) {
      dscR[i]->getAtom(j).setB(scoreParts[i].first);
      abdR[i]->getAtom(j).setB(scoreParts[i].second);
      sscR[i]->getAtom(j).setB(structScores[i]);
    }
  }
  string prefix = op.getString("out");
  dscS.writePDB(prefix + ".dsc.pdb");
  abdS.writePDB(prefix + ".abd.pdb");
  sscS.writePDB(prefix + ".ssc.pdb");

  // Store the three sets of scores in tabular format, if specified
  if (op.isGiven("tab")) {
    ofstream outputFS(op.getString("tab"));
    outputFS << "Position | design score | abundance score | structure score" << fixed << setprecision(6) << endl;
    for (int i = 0; i < numRes; i++) {
      outputFS << sscR[i]->getChainID() << "," << sscR[i]->getNum() << "\t" << scoreParts[i].first << "\t" << scoreParts[i].second << "\t" << structScores[i] << endl;
    }
    outputFS.close();
  }

  return 0;
}
