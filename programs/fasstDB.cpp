#include "msttypes.h"
#include "mstoptions.h"
#include "mstfasst.h"
#include "mstrotlib.h"
#include "mstsystem.h"
#include "mstcondeg.h"
#include <chrono>

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Creates a FASST database from input PDB files. Options:");
  op.addOption("pL", "a file with a list of PDB files.");
  op.addOption("dL", "a file with a list of FASST databases (will consolidate into one).");
  op.addOption("o", "output database file name.", true);
  op.addOption("m", "memory save flag (will store backbone only).");
  op.addOption("c", "clean up PDB files, so that only protein residues with enough of a backbone to support rotamer building survive.");
  op.addOption("pp", "store phi/psi properties in the database.");
  op.addOption("env", "store residue freedom property in the database (option value must be a path to an MST rotamer library file).");
  op.addOption("batch", "an integer. If specified, instead of building the database will spread all the work across this many "
                        "jobs, writing corresponding batch files for submission to the cluster. Will also produce a file called "
                        "<out>.fin.sh (where <out> is the base of the name specified in --o), which is to be run after all jobs "
                        "finish to complete the database building process.");
  op.setOptions(argc, argv);
  if (!op.isGiven("pL") && !op.isGiven("dL")) MstUtils::error("either --pL or --dL must be given!");
  if (op.isGiven("env")) MstUtils::assert(!op.getString("env").empty() && MstSys::fileExists(op.getString("env")), "--env must point to a rotamer library file");
  RotamerLibrary RL;
  if (op.isGiven("env")) RL.readRotamerLibrary(op.getString("env"));

  if (!op.isGiven("batch")) {
    FASST S;
    cout << "Reading structures..." << endl;
    S.setMemorySaveMode(op.isGiven("m"));
    if (op.isGiven("pL")) {
      vector<string> pdbFiles = MstUtils::fileToArray(op.getString("pL"));
      for (int i = 0; i < pdbFiles.size(); i++) {
        Structure P(pdbFiles[i]);
        if (op.isGiven("c")) {
          Structure C; RotamerLibrary::extractProtein(C, P);
          if (P.residueSize() != C.residueSize()) {
            cout << pdbFiles[i] << ", had " << P.residueSize() << " residues, and " << C.residueSize() << " residues after cleaning..." << endl;
          }
          C.setName(P.getName()); P = C;
        }
        S.addTarget(P);
      }
    }
    if (op.isGiven("dL")) {
      vector<string> dbFiles = MstUtils::fileToArray(op.getString("dL"));
      for (int i = 0; i < dbFiles.size(); i++) {
        S.readDatabase(dbFiles[i]);
      }
    }
    if (op.isGiven("pp") || op.isGiven("env")) {
      // compute and add some properties
      for (int ti = 0; ti < S.numTargets(); ti++) {
        Structure P = S.getTarget(ti);
        if (op.isGiven("pp")) {
          vector<Residue*> residues = P.getResidues();
          vector<mstreal> phi(residues.size()), psi(residues.size());
          for (int ri = 0; ri < residues.size(); ri++) {
            phi[ri] = residues[ri]->getPhi(false);
            psi[ri] = residues[ri]->getPsi(false);
          }
          S.addResidueProperties(ti, "phi", phi);
          S.addResidueProperties(ti, "psi", psi);
        }
        if (op.isGiven("env")) {
          ConFind C(&RL, P);
          vector<Residue*> residues = P.getResidues();
          vector<mstreal> freedoms = C.getFreedom(residues);
          S.addResidueProperties(ti, "env", freedoms);
        }
      }
    }
    S.writeDatabase(op.getString("o"));
  } else {
    if (!op.isGiven("pL")) MstUtils::error("--pL must be given with --batch");
    if (!op.isInt("batch") || (op.getInt("batch") <= 0)) MstUtils::error("--batch must be a positive integer!");
    int nJ = op.getInt("batch");
    vector<string> pdbFiles = MstUtils::fileToArray(op.getString("pL"));
    srand(time(NULL) + (int) getpid());
    MstUtils::shuffle(pdbFiles);
    vector<pair<int, int> > tasks = MstUtils::splitTasks(pdbFiles.size(), nJ);
    fstream outf;
    string dbListFile = MstSys::pathBase(op.getString("o")) + ".dL";
    fstream dblf; MstUtils::openFile(dblf, dbListFile, ios::out);
    vector<string> toClean;
    for (int i = 0; i < nJ; i++) {
      string base = MstSys::pathBase(op.getString("o")) + "." + MstUtils::toString(i);

      // dump subset of PDB files this job will work on
      string listFile = base + ".list";
      MstUtils::openFile(outf, listFile, ios::out);
      for (int k = tasks[i].first; k <= tasks[i].second; k++) {
        outf << pdbFiles[k] << endl;
      }
      outf.close();

      // create a job script file that will work on this subset
      string batchFile = base + ".sh";
      string dbFile = base + ".db";
      MstUtils::openFile(outf, batchFile, ios::out);
      outf << "#!/bin/bash\n" << "#$ -j y\n" << "#$ -cwd\n" << "#$ -V\n";
      outf << "#$ -l vf=2G\n";
      int hrs = (int) ceil(5*(tasks[i].second - tasks[i].first + 1)/60.0); // five minutes per structure should be plenty
      outf << "#$ -l h_rt=" << hrs << ":00:00\n";
      outf << op.getExecName() << " --pL " << listFile << " --o " << dbFile;
      // keep all other options from the call to self
      vector<string> allOpts = op.getAllGivenOptions();
      for (int j = 0; j < allOpts.size(); j++) {
        if ((allOpts[j].compare("batch") == 0) || (allOpts[j].compare("pL") == 0) || (allOpts[j].compare("o") == 0)) continue;
        outf << " --" << allOpts[j] << " " << op.getString(allOpts[j]);
      }
      outf << endl;
      outf.close();
      dblf << dbFile << endl;
      toClean.push_back(base + ".sh*"); toClean.push_back(base + ".db");
    }
    dblf.close();
    fstream fin; MstUtils::openFile(fin, MstSys::pathBase(op.getString("o")) + ".fin.sh", ios::out);
    fin << op.getExecName() << " --dL " << dbListFile << " --o " << op.getString("o") << endl;
    fin << "if [ $? -eq 0 ]; then # only clean up if database creation was successful" << endl;
    for (int k = 0; k < toClean.size(); k++) {
      fin << "  rm " << toClean[k] << endl;
    }
    fin << "fi" << endl;
    fin.close();
  }
  return 0; // success
}
