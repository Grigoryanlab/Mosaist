#include "msttypes.h"
#include "mstoptions.h"
#include "mstfasst.h"
#include "mstsystem.h"
#include <chrono>

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Implements the FASST (FAst Structure Search Algorithm). Options:");
  op.addOption("q", "query PDB file.", true);
  op.addOption("d", "a database file with a list of PDB files.");
  op.addOption("b", "a binary database file. If both --d and --b are given, will overwrite this file with a corresponding binary database.");
  op.addOption("r", "RMSD cutoff (takes the size-dependent cutoff by default).");
  op.addOption("red", "set redundancy cutoff level in percent (default is 100, so no redundancy filtering).");
  op.addOption("redProp", "set redundancy property name. If defined, will assume the FASST database encodes this relational property and will define redundancy via it.");
  op.addOption("min", "min number of matches.");
  op.addOption("max", "max number of matches.");
  op.addOption("const", "specify sequence constraints as a semicolon-separated list of specifications. E.g., '0 3 LEU,ALA' means residue index 3 from the first query segment must be either LEU or ALA. Or, '0 3 LEU,ALA; 1 2 LYS' additionally specifies that residue index 2 from the second segment must be LYS.");
  op.addOption("outType", "what portion of matching sequences to output. Default is 'region', which refers to just the matching region. Also possible are: 'full' (for full structure) and 'withGaps' (for the matching regions plus any gaps between segments).");
  op.addOption("strOut", "dump structures into this directory.");
  op.addOption("seqOut", "sequence output file.");
  op.addOption("sc", "dump sidechains (not only the backbone).");
  op.addOption("pp", "store phi/psi properties in the database, if creating a new one from PDB files.");
  op.setOptions(argc, argv);
  int memInit = MstSys::memUsage();
  if (op.isGiven("redProp")) MstUtils::assertCond(!op.getString("redProp").empty(), "--redProp must specify a property name");
  FASST::matchType type = FASST::matchType::REGION;
  if (op.isGiven("outType")) {
    if (op.getString("outType").compare("region") == 0) {
    } else if (op.getString("outType").compare("full") == 0) {
      type = FASST::matchType::FULL;
    } else if (op.getString("outType").compare("withGaps") == 0) {
      type = FASST::matchType::WITHGAPS;
    } else {
      MstUtils::error("unknown output type '" + op.getString("outType") + "'");
    }
  }
  if (op.isGiven("strOut") && !MstSys::isDir(op.getString("strOut"))) MstSys::cmkdir(op.getString("strOut"));

  FASST S;
  cout << "Reading the database..." << endl;
  auto begin = chrono::high_resolution_clock::now();
  Structure query(op.getString("q"));
  S.setQuery(query);
  if (op.isGiven("d")) {
    vector<string> pdbFiles = MstUtils::fileToArray(op.getString("d"));
    for (int i = 0; i < pdbFiles.size(); i++) {
      Structure P(pdbFiles[i]);
      S.addTarget(P);
      // compute and add some properties
      if (op.isGiven("pp")) {
        vector<Residue*> residues = P.getResidues();
        vector<mstreal> phi(residues.size()), psi(residues.size());
        for (int ri = 0; ri < residues.size(); ri++) {
          phi[ri] = residues[ri]->getPhi(false);
          psi[ri] = residues[ri]->getPsi(false);
        }
        S.addResidueProperties(S.numTargets() - 1, "phi", phi);
        S.addResidueProperties(S.numTargets() - 1, "psi", psi);
      }
    }
    if (op.isGiven("b")) {
      S.writeDatabase(op.getString("b"));
    }
  } else if (op.isGiven("b")) {
    S.readDatabase(op.getString("b"), 2);
  } else {
    MstUtils::error("either --b or --d must be given!");
  }
  if (op.isGiven("r")) { S.setRMSDCutoff(op.getReal("r")); }
  else {
    cout << "setting RMSD cutoff to " << RMSDCalculator::rmsdCutoff(query) << endl;
    S.setRMSDCutoff(RMSDCalculator::rmsdCutoff(query));
  }
  S.setMaxNumMatches(op.getInt("max", -1));
  S.setMinNumMatches(op.getInt("min", -1));
  // S.setMaxGap(1, 0, 6); S.setMinGap(1, 0, 0);
  S.setRedundancyCut(op.getReal("red", 100.0)/100.0);
  if (op.isGiven("redProp")) S.setRedundancyProperty(op.getString("redProp"));
  fasstSeqConstSimple seqConst(S.getNumQuerySegments());
  if (op.isGiven("const")) {
    vector<string> cons = MstUtils::split(op.getString("const"), ";");
    for (int i = 0; i < cons.size(); i++) {
      vector<string> con = MstUtils::split(MstUtils::trim(cons[i]), " ");
      MstUtils::assertCond(con.size() == 3, "could not parse constraint '" + cons[i] + "' from constraints line " + op.getString("const"));
      int segIdx = MstUtils::toInt(con[0]);
      int resIdx = MstUtils::toInt(con[1]);
      vector<string> aas = MstUtils::split(MstUtils::trim(con[2]), "|");
      seqConst.addConstraint(segIdx, resIdx, aas);
    }
    S.options().setSequenceConstraints(seqConst);
  }
  auto end = chrono::high_resolution_clock::now();
  cout << "DB reading took " << chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms" << endl;
  cout << "memory usage: " << MstSys::memUsage() - memInit << " KB" << endl;
  cout << "Searching..." << endl;
  begin = chrono::high_resolution_clock::now();
  S.search();
  end = chrono::high_resolution_clock::now();
  cout << "Search took " << chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms" << endl;
  cout << "found " << S.numMatches() << " matches:" << endl;
  cout << "memory usage: " << MstSys::memUsage() << " KB" << endl;
  fasstSolutionSet matches = S.getMatches(); int i = 0;
  vector<vector<mstreal> > phi, psi;
  for (auto it = matches.begin(); it != matches.end(); ++it, ++i) {
    cout << S.toString(*it) << endl;
    cout << *it << endl;
    if (op.isGiven("pp")) {
      if (S.isResiduePropertyDefined("phi")) cout << "\tphi: " << MstUtils::vecToString(S.getResidueProperties(*it, "phi")) << endl;
      if (S.isResiduePropertyDefined("psi")) cout << "\tpsi: " << MstUtils::vecToString(S.getResidueProperties(*it, "psi")) << endl;
    }
    if (op.isGiven("strOut")) {
      // Structure match = S.getMatchStructure(*it, true, FASST::matchType::FULL);
      Structure match = S.getMatchStructure(*it, op.isGiven("sc"), type);
      match.writePDB(op.getString("strOut") + "/match" + MstUtils::toString(i) + ".pdb");
    }
  }
  fstream of;
  if (op.isGiven("seqOut")) MstUtils::openFile(of, op.getString("seqOut"), ios::out);
  for (auto it = matches.begin(); it != matches.end(); ++it, ++i) {
    Sequence seq = S.getMatchSequence(*it);
    cout << seq.toString() << endl;
    if (op.isGiven("seqOut")) of << seq.toString() << endl;
  }
  if (op.isGiven("seqOut")) of.close();
}
