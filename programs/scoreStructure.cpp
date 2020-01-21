#include "msttypes.h"
#include "mstoptions.h"
#include "msttermanal.h"
#include "mstrotlib.h"

int main(int argc, char** argv) {
  // Setup and get the input arguments
  MstOptions op;
  op.setTitle("Scores the input structure using TERMANAL and outputs three structures whose B-factors are annotated with the design, abundance, and structure scores");
  op.addOption("p", "input PDB file to score.", true);
  op.addOption("sel", "selection within the input PDB to score.");
  op.addOption("o", "output prefix for the B-factor annotated output structures (prefix.{dsc,abd,ssc}.pdb).", true);
  op.addOption("tab", "output file containing the three sets of scores in tabular format.");
  op.addOption("db", "binary FASST database to search for matches in.", true);
  op.addOption("rLib", "rotamer library to compute contact degree with.", true);
  op.addOption("v", "if true, prints scoring details as they are computed (default=false).");
  op.addOption("cd", "contact degree cutoff when forming TERMs (default=0.1).");
  op.addOption("pad", "number of residues to pad contacting residues with when forming TERMs (default=2).");
  op.addOption("rmsd", "RMSD cutoff when searching for matches (default=2.0).");
  op.addOption("count", "maximum number of matches to search for within the specified RMSD cutoff (default=50).");
  op.addOption("ps", "pseudocount to add to the smoothed scores before computing their logs (default=0.01).");
  op.setOptions(argc, argv);

  // Read in the structure, FASST database, and rotamer library
  Structure So(op.getString("p")), S;
  RotamerLibrary::extractProtein(S, So);
  selector sele(S);
  vector<Residue*> selR = op.isGiven("sel") ? sele.selectRes(op.getString("sel")) : S.getResidues();
  FASST F; F.readDatabase(op.getString("db"));
  TERMANAL T(&F); T.readRotamerLibrary(op.getString("rLib"));
  if (op.isGiven("cd")) T.setCDCut(op.getReal("cd"));
  if (op.isGiven("pad")) T.setPad(op.getInt("pad"));
  if (op.isGiven("rmsd")) T.setRMSDCut(op.getReal("rmsd"));
  if (op.isGiven("count")) T.setMatchCount(op.getInt("count"));
  if (op.isGiven("ps")) T.setPseudocount(op.getReal("ps"));

  // Score each residue in the structure
  vector<pair<mstreal, mstreal>> scoreParts;
  vector<mstreal> structScores = T.scoreStructure(S, selR, &scoreParts, op.isGiven("v"));

  // Zero out B-factors
  vector<Atom*> allAtoms = S.getAtoms();
  for (Atom* a : allAtoms) a->setB(0.0);

  // Create the B-factor annotated output PDBs using the computed scores
  Structure dscS(S), abdS(S), sscS(S);
  for (int i = 0; i < selR.size(); i++) {
    Residue* r = selR[i];
    int ri = r->getResidueIndex();
    Residue& dscR = dscS.getResidue(ri);
    Residue& abdR = abdS.getResidue(ri);
    Residue& sscR = sscS.getResidue(ri);
    for (int j = 0; j < dscR.atomSize(); j++) dscR.getAtom(j).setB(-scoreParts[i].first);
    for (int j = 0; j < abdR.atomSize(); j++) abdR.getAtom(j).setB(-scoreParts[i].second);
    for (int j = 0; j < sscR.atomSize(); j++) sscR.getAtom(j).setB(-structScores[i]);
  }
  string prefix = op.getString("o");
  dscS.writePDB(prefix + ".dsc.pdb");
  abdS.writePDB(prefix + ".abd.pdb");
  sscS.writePDB(prefix + ".ssc.pdb");

  // Store the three sets of scores in tabular format, if specified
  if (op.isGiven("tab")) {
    ofstream outputFS(op.getString("tab"));
    outputFS << "Position | design score | abundance score | structure score" << fixed << setprecision(6) << endl;
    for (int i = 0; i < selR.size(); i++) {
      Residue* r = selR[i];
      outputFS << r->getChainID() << "," << r->getNum() << "\t" << -scoreParts[i].first << "\t" << -scoreParts[i].second << "\t" << -structScores[i] << endl;
    }
    outputFS.close();
  }

  return 0;
}
