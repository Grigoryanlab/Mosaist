#include "msttypes.h"
#include "mstoptions.h"
#include "mstsequence.h"

// extracts just the protein part of the input structure and also reports what
// fraction of residues in the entry are "standard" things (i.e., not hetero
// atoms), but not recognized as standard amino acids
double proteinOnly(System& So, System& S, vector<string>& legalNames);

int main(int argc, char *argv[]) {
  // ------- Input parameters
  MstOptions op;
  op.setTitle("Creates a complete database for use in TERM-based protein design. Options:");
  op.addOption("p", "root path to local PDB database (e.g., /home/grigoryanlab/library/databases/pdb.org).", true);
  op.addOption("o", "path to output directory, where everything will be placed.", true);
  op.addOption("se", "extract statistical energies, in addition to building a search database.");
  op.addOption("r", "X-ray resolution cutoff (default 2.6).");
  op.addOption("maxS", "largest biounit size to consider (default 5000).");
  op.setOptions(argc, argv);
  string pdbBase = op.getString("p");
  string outPath = op.getString("o");
  double resolCut = op.getReal("r", 2.6);
  int maxSize = op.getInt("maxS", 5000);
  string tmpDir = "/tmp";

  // ------- Legal protein residue names
  vector<string> legalNames;
  legalNames.push_back("ALA"); legalNames.push_back("CYS"); legalNames.push_back("ASP"); legalNames.push_back("GLU"); legalNames.push_back("PHE"); legalNames.push_back("GLY");
  legalNames.push_back("HIS"); legalNames.push_back("ILE"); legalNames.push_back("LYS"); legalNames.push_back("LEU"); legalNames.push_back("MET"); legalNames.push_back("ASN");
  legalNames.push_back("PRO"); legalNames.push_back("GLN"); legalNames.push_back("ARG"); legalNames.push_back("SER"); legalNames.push_back("THR"); legalNames.push_back("VAL");
  legalNames.push_back("TRP"); legalNames.push_back("TYR"); legalNames.push_back("HSD"); legalNames.push_back("HSE"); legalNames.push_back("HSC"); legalNames.push_back("HSP");
  legalNames.push_back("MSE");
  legalNames.push_back("CSO"); legalNames.push_back("HIP"); legalNames.push_back("SEC"); legalNames.push_back("SEP"); legalNames.push_back("TPO"); legalNames.push_back("PTR");

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
    vector<string> ents = MstUtils::split(line, " ;\t");
    if (ents.size() != 2) {
      cout << "\tskipping line '" << line << "' of PDB resolution file..." << endl;
      continue;
    }
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

  // filter by biounit size and by fraction of protein residues
  string outPdbDir = outPath + "/PDB"; MstUtils::csystem("mkdir " + outPdbDir);
  string tmpFile = tmpDir + "/db-entry-mek.pdb";
  fstream ofh; MstUtils::openFile(ofh, outPath + "/PDB.fasta");
  oldIDs = vector<string>(); oldIDs.swap(IDs);
  for (int i = 0; i < oldIDs.size(); i++) {
    // look at the first bio unit, since don't know any better
    string m2 = oldIDs[i].substr(1, 2);
    string base = m2 + MstUtils::lc(oldIDs[i]);
    string packedFile = pdbBase + "/data/biounit/coordinates/divided/" + base + ".pdb1.gz";
    MstUtils::csystem("gunzip < " + packedFile + " > " + tmpFile);
    Structure S(tmpFile);

    // filter by size and number of chains
    if (S.residueSize() > maxSize) continue;
    if (S.chainSize() > 26) continue; // don't want to get into alphabet soup issues

    // at least 60% protein (not counting water, metals, or unnatural amino acids)
    Structure Sp;
    if (proteinOnly(S, Sp, legalNames) < 0.6) continue;

    // dump the "cleaned" PDB file and the corresponding sequence
    if (!MstUtils::isDir(outPdbDir + "/" + m2)) MstUtils::csystem("mkdir " + outPdbDir + "/" + m2);
    Sp.writePDB(outPdbDir + "/" + m2 + MstUtils::uc(oldIDs[i]) + ".pdb");
    for (int ci = 0; ci < Sp.chainSize(); ci++) {
      Chain& C = Sp[ci];
      ofh << ">" << MstUtils::uc(oldIDs[i]) << "_" << C.getID() << endl;
      for (int ri = 0; ri < C.residueSize(); ri++) {
        ofh << Sequence::toSingle(C[ri].getName());
      }
      ofh << endl;
    }
    IDs.push_back(oldIDs[i]);
  }
  ofh.close();
  printf("after filtering for size and protein content, and cleaning, %d entries left...\n", (int) IDs.size());
  MstUtils::csystem("rm " + tmpFile);
}

double proteinOnly(System& So, System& S, vector<string>& legalNames) {
  AtomPointerVector A;
  int numNonProtRes = 0;
  for (int i = 0; i < So.chainSize(); i++) {
    Chain& chain = So[i];
    for (int j = 0; j < chain.residueSize(); j++) {
      Residue& res = chain[j];
      for (int k = 0; k < legalNames.size(); k++) {
        if (res.isNamed(legalNames[k])) {
          for (int ai = 0; ai < res.atomSize(); ai++) A.push_back(&(res[ai]));
          break;
        } else {
          // these are things that look "normal" (e.g., DNA, RNA, protein), but
          // not recognized as standard protein amino acids
          if ((res.atomSize() > 0) && (!res[0].isHetero())) numNonProtRes++;
        }
      }
    }
  }
  S.addAtoms(&A);
  return numNonProtRes * 1.0 / S.residueSize();
}
