#include "msttypes.h"
#include "mstoptions.h"

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Creates a complete database for use in TERM-based protein design. Options:");
  op.addOption("p", "root path to local PDB database (e.g., /home/grigoryanlab/library/databases/pdb.org).", true);
  op.addOption("o", "path to output directory, where everything will be placed.", true);
  op.addOption("se", "extract statistical energies, in addition to building a search database.", true);
  op.addOption("r", "X-ray resolution cutoff (default 2.6).", true);
  op.setOptions(argc, argv);
  string pdbBase = op.getString("p");
  double resolCut = op.getReal("r", 2.6);

  /* Filter all biological units */
  vector<string> entriesAll, IDs;
  MstUtils::fileToArray(pdbBase + "/derived_data/pdb_entry_type.txt", entriesAll);
  printf("started with %d entries in the database...\n", (int) entriesAll.size());

  // make sure they contain proteins and are x-ray structures
  for (int i = 0; i < entriesAll.size(); i++) {
    vector<string> line = MstUtils::split(entriesAll[i], " \t");
    MstUtils::assert(line.size() == 3, "could not parse line '" + entriesAll[i] + "' from PDB entry type file.");
    if (line[1].compare("prot") && line[1].compare("prot-nuc")) continue;
    if (line[2].compare("diffraction")) continue;
    IDs.push_back(MstUtils::uc(line[0]));
  }
  printf("after filtering for protein-containing X-ray structures, %d entries left...\n", (int) IDs.size());

  // get resolutions
  map<string, double> resol;
  fstream ifh; string line; bool f = false;
  MstUtils::openFile (ifh, pdbBase + "/derived_data/index/resolu.idx");
  while (getline(ifh, line)) {
    if (line.find("----") == 0) { f = true; continue; }
    if (!f) { continue; }
    vector<string> ents = MstUtils::split(line, " ;");
    MstUtils::assert(ents.size() == 3, "could not parse line '" + line + "' from PDB resolution file.");
    resol[MstUtils::uc(ents[0])] = MstUtils::toReal(ents[1]);
  }
  ifh.close();

  // filter by resolution
  vector<string> oldIDs;
  oldIDs.swap(IDs);
  for (int i = 0; i < oldIDs.size(); i++) {
    if (resol.find(oldIDs[i]) == resol.end()) MstUtils::error("could not find resolution for entry " + oldIDs[i]);
    if (resol[oldIDs[i]] > resolCut) continue;
    IDs.push_back(oldIDs[i]);
  }
  printf("after filtering for resolution better than %f, %d entries left...\n", resolCut, (int) IDs.size());
}
