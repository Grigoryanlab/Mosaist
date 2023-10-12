#include "msttypes.h"
#include "mstoptions.h"
#include "mstfasst.h"
#include "mstsystem.h"
#include <chrono>

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Command-line acceess to the FASST (FAst Structure Search Algorithm) method. Options:");
  op.addOption("q", "query PDB file.", true);
  op.addOption("d", "a database file with a list of PDB files.");
  op.addOption("b", "a binary database file. If both --d and --b are given, the combined database will be searched.");
  op.addOption("r", "RMSD cutoff (takes the size-dependent cutoff by default).");
  op.addOption("red", "set redundancy cutoff level in percent (default is 100, so no redundancy filtering).");
  op.addOption("redProp", "set redundancy property name. If defined, will assume the FASST database encodes this relational property and will define redundancy via it.");
  op.addOption("min", "min number of matches.");
  op.addOption("max", "max number of matches.");
  op.addOption("seqConst", "specify sequence constraints as a semicolon-separated list of specifications. E.g., '0 3 LEU,ALA' means residue index 3 from the first query segment must be either LEU or ALA. Or, '0 3 LEU,ALA; 1 2 LYS' additionally specifies that residue index 2 from the second segment must be LYS.");
  op.addOption("gapConst", "specify gap constraints as a semicolon-separated list of specifications. E.g., '0 1 0 10; 1 2 4 6' that the gap between segments 1 and 2 must be between 0 and 10 residues and the gap between segments 2 and 3 must be between 4 and 6 residues.");
  op.addOption("outType", "what portion of matching sequences to output. Default is 'region', which refers to just the matching region. Also possible are: 'full' (for full structure) and 'withGaps' (for the matching regions plus any gaps between segments).");
  op.addOption("strOut", "dump structures into this directory.");
  op.addOption("seqOut", "sequence output file.");
  op.addOption("matchOut", "match output file.");
  op.addOption("m", "memory saving mode: 0 means does not do any memory savings; 1 means strip the side-chains; 2 (default) means destroy the original target structure upon reading, and only keep backbone coordinates.");
  op.addOption("sc", "dump sidechains (not only the backbone).");
  op.setOptions(argc, argv);
  int memInit = MstSys::memUsage();
  if (op.isGiven("redProp")) MstUtils::assertCond(!op.getString("redProp").empty(), "--redProp must specify a property name");
  if (!op.isGiven("b") && !op.isGiven("d")) MstUtils::error("either --b or --d must be given!");
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
  if (op.isGiven("b")) {
    S.readDatabase(op.getString("b"), op.getInt("m", 2));
  }
  if (op.isGiven("d")) {
    vector<string> pdbFiles = MstUtils::fileToArray(op.getString("d"));
    for (int i = 0; i < pdbFiles.size(); i++) {
      Structure P(pdbFiles[i]);
      S.addTarget(P);
    }
  }
  if (op.isGiven("r")) { S.setRMSDCutoff(op.getReal("r")); }
  else {
    cout << "setting RMSD cutoff to " << RMSDCalculator::rmsdCutoff(query) << endl;
    S.setRMSDCutoff(RMSDCalculator::rmsdCutoff(query));
  }
  S.setMaxNumMatches(op.getInt("max", -1));
  S.setMinNumMatches(op.getInt("min", -1));
  S.setRedundancyCut(op.getReal("red", 100.0)/100.0);
  if (op.isGiven("redProp")) S.setRedundancyProperty(op.getString("redProp"));
  fasstSeqConstSimple seqConst(S.getNumQuerySegments());
  if (op.isGiven("seqConst")) {
    vector<string> cons = MstUtils::split(op.getString("seqConst"), ";");
    for (int i = 0; i < cons.size(); i++) {
      vector<string> con = MstUtils::split(MstUtils::trim(cons[i]), " ");
      MstUtils::assertCond(con.size() == 3, "could not parse constraint '" + cons[i] + "' from constraints line " + op.getString("seqConst"));
      int segIdx = MstUtils::toInt(con[0]);
      int resIdx = MstUtils::toInt(con[1]);
      vector<string> aas = MstUtils::split(MstUtils::trim(con[2]), "|");
      seqConst.addConstraint(segIdx, resIdx, aas);
    }
    S.options().setSequenceConstraints(seqConst);
  }
  if (op.isGiven("gapConst")) {
    vector<string> cons = MstUtils::split(op.getString("gapConst"), ";");
    for (int i = 0; i < cons.size(); i++) {
      vector<string> con = MstUtils::split(MstUtils::trim(cons[i]), " ");
      MstUtils::assertCond(con.size() == 4, "could not parse constraint '" + cons[i] + "' from constraints line " + op.getString("gapConst"));
      int segI = MstUtils::toInt(con[0]);
      int segJ = MstUtils::toInt(con[1]);
      int minGapLen = MstUtils::toInt(con[2]);
      int maxGapLen = MstUtils::toInt(con[3]);
      S.setMaxGap(segI, segJ, maxGapLen);
      S.setMinGap(segI, segJ, minGapLen);
    }
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
    cout << *it << endl;
    if (op.isGiven("strOut")) {
      // Structure match = S.getMatchStructure(*it, true, FASST::matchType::FULL);
      Structure match = S.getMatchStructure(*it, op.isGiven("sc"), type);
      match.writePDB(op.getString("strOut") + "/match" + MstUtils::toString(i) + ".pdb");
    }
  }
  fstream of, mof;
  if (op.isGiven("seqOut")) MstUtils::openFile(of, op.getString("seqOut"), ios::out);
  if (op.isGiven("matchOut")) MstUtils::openFile(mof, op.getString("matchOut"), ios::out);
  for (auto it = matches.begin(); it != matches.end(); ++it, ++i) {
    Sequence seq = S.getMatchSequence(*it);
    cout << seq.toString() << endl;
    if (op.isGiven("seqOut")) of << seq.toString() << endl;
    if (op.isGiven("matchOut")) mof << S.toString(*it) << endl;
  }
  if (op.isGiven("seqOut")) of.close();
  if (op.isGiven("matchOut")) mof.close();
}
