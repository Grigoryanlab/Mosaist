#include "msttypes.h"
#include "mstoptions.h"
#include "mstfasst.h"
#include "mstrotlib.h"
#include "mstsystem.h"
#include "mstcondeg.h"
#include "mstsequence.h"
#include <chrono>

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Creates a FASST database from input PDB files. Options:");
  op.addOption("pL", "a file with a list of PDB files.");
  op.addOption("db", "a previously-written FASST database.");
  op.addOption("dL", "a file with a list of FASST databases (will consolidate into one).");
  op.addOption("o", "output database file name.", true);
  op.addOption("m", "memory save flag (will store backbone only).");
  op.addOption("c", "clean up PDB files, so that only protein residues with enough of a backbone to support rotamer building survive.");
  op.addOption("pp", "store phi/psi properties in the database.");
  op.addOption("env", "store residue freedom property in the database. If this is give, --rLib must also be given.");
  op.addOption("cont", "store inter-residue contact information (for all residue pairs with contact degrees below the specified limit). If this is given, --rLib must also be given.");
  op.addOption("sim", "percent sequence identity cutoff. If specified, will store local-window sequence similarity between all pairs of positions in the database, using this cutoff.");
  op.addOption("win", "window size to use with the similarity searching with --sim; must be an odd integer. Default is 31 (i.e., +/- 15 from the residue in question).");
  op.addOption("rLib", "path to an MST rotamer library file.");
  op.addOption("batch", "an integer. If specified, instead of building the database will spread all the work across this many "
                        "jobs, writing corresponding batch files for submission to the cluster. Will also produce a file called "
                        "<out>.fin.sh (where <out> is the base of the name specified in --o), which is to be run after all jobs "
                        "finish to complete the database building process.");
  op.setOptions(argc, argv);
  RotamerLibrary RL;
  if (!op.isGiven("pL") && !op.isGiven("dL") && !op.isGiven("db")) MstUtils::error("either --pL, --dL, or --db must be given!");
  if (op.isGiven("rLib")) MstUtils::assert(MstSys::fileExists(op.getString("rLib")), "--rLib is not a valid file path");
  if (op.isGiven("cont") && (!op.isReal("cont") || (op.getReal("cont") < 0) || (op.getReal("cont") > 1))) MstUtils::error("--cont must be a real in range [0; 1]");
  if ((op.isGiven("env") || op.isGiven("cont"))) {
    if (!op.isGiven("rLib")) MstUtils::error("--rLib is needed, but not given");
    RL.readRotamerLibrary(op.getString("rLib"));
  }
  if (op.isGiven("sim") && (!op.isReal("sim") || (op.getReal("sim") < 0) || (op.getReal("sim") > 100))) MstUtils::error("--sim must be a non-negative value below 100.");
  if (op.isGiven("win") && (!op.isInt("win") || (op.getInt("win") <= 0) || (op.getInt("win") % 2 == 0))) MstUtils::error("--win must be a positive odd integer.");

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
    if (op.isGiven("db")) {
      S.readDatabase(op.getString("db"));
    }
    if (op.isGiven("dL")) {
      vector<string> dbFiles = MstUtils::fileToArray(op.getString("dL"));
      for (int i = 0; i < dbFiles.size(); i++) {
        S.readDatabase(dbFiles[i]);
      }
    }
    if (op.isGiven("pp") || op.isGiven("env") || op.isGiven("cont")) {
      cout << "Computing per-target residue properties..." << endl;
      // compute and add some properties
      for (int ti = 0; ti < S.numTargets(); ti++) {
        cout << "\ttarget " << ti+1 << "/" << S.numTargets() << "..." << endl;
        Structure P = S.getTargetCopy(ti);
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
        if (op.isGiven("env") || op.isGiven("cont")) {
          ConFind C(&RL, P); // both need the confind object
          // environment
          vector<Residue*> residues = P.getResidues();
          vector<mstreal> freedoms = C.getFreedom(residues);
          S.addResidueProperties(ti, "env", freedoms);

          // contacts
          mstreal cdcut = op.getReal("cont");
          contactList list = C.getContacts(P, cdcut);
          map<int, map<int, mstreal> > conts;
          for (int i = 0; i < list.size(); i++) {
            int rA = list.residueA(i)->getResidueIndex();
            int rB = list.residueB(i)->getResidueIndex();
            conts[rA][rB] = list.degree(i);
            conts[rB][rA] = list.degree(i);
          }
          S.addResiduePairProperties(ti, "conts", conts);
        }
      }
    }
    if (op.isGiven("sim")) {
      cout << "Computing local-window sequence similarity..." << endl;
      // first cluster all local windows
      cout << "\tgathering local sequence windows..." << endl;
      int L = op.getInt("win", 31);
      int L2 = (L - 1)/2;
      vector<Sequence> wins;
      vector<int> winTarg, winStart;
      vector<bool> nTerm, cTerm;
      for (int ti = 0; ti < S.numTargets(); ti++) {
        Structure P = S.getTargetCopy(ti);
        int ri = 0;
        for (int i = 0; i < P.chainSize(); i++) {
          Chain& C = P[i];
          for (int j = 0; j < C.residueSize() - L + 1; j++, ri++) {
            Sequence win(L);
            for (int k = 0; k < L; k++) win[k] = SeqTools::aaToIdx(C[j + k].getName());
            wins.push_back(win);
            winTarg.push_back(ti); winStart.push_back(ri);
            nTerm.push_back(j == 0); cTerm.push_back(j == C.residueSize() - L);
          }
        }
      }
      cout << "\tclustering " << wins.size() << " windows at " << op.getInt("sim") << "\% sequence identity..." << endl;
      vector<vector<int> > clusts = SeqTools::rSearch(wins, op.getInt("sim")/100.0, 0.999, true);
      cout << "\tfound " << clusts.size() << " clusters..." << endl;

      // then visit all similar pairs and mark what residues that makes similar
      int symN = 0;
      for (int i = 0; i < clusts.size(); i++) {
        vector<int>& C = clusts[i];
        int ti = winTarg[i];
        int ri = winStart[i] + L2;
        for (int k = 0; k < C.size(); k++) {
          int j = C[k];
          // windows i and j are similar
          int tj = winTarg[j];
          int rj = winStart[j] + L2;
          S.addResidueRelationship(ti, "sim", ri, tj, rj);
          symN++;
          // if either window starts at the N-terminus, then the whole first half is redundant
          if (nTerm[i] || nTerm[j]) {
            for (rj = winStart[j]; rj < winStart[j] + L2 - 1; rj++) S.addResidueRelationship(ti, "sim", ri, tj, rj);
            symN += L2;
          }
          // if either window ends at the C-terminus, then the whole second half is redundant
          if (cTerm[i] || cTerm[j]) {
            for (rj = winStart[j] + L2 + 1; rj < winStart[j] + L; rj++) S.addResidueRelationship(ti, "sim", ri, tj, rj);
            symN += L2;
          }
        }
      }
      cout << "\trecorded " << symN << " similar windows" << endl;
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
      outf << "#$ -l vf=2G\n" << "#$ -l ironfs\n";
      int hrs = (int) ceil(5*(tasks[i].second - tasks[i].first + 1)/60.0); // five minutes per structure should be plenty
      outf << "#$ -l h_rt=" << hrs << ":00:00\n";
      outf << op.getExecName() << " --pL " << listFile << " --o " << dbFile;
      // keep all other options from the call to self
      vector<string> allOpts = op.getAllGivenOptions();
      for (int j = 0; j < allOpts.size(); j++) {
        if ((allOpts[j].compare("batch") == 0) || (allOpts[j].compare("pL") == 0) || (allOpts[j].compare("o") == 0) || (allOpts[j].compare("sim") == 0)) continue;
        outf << " --" << allOpts[j] << " " << op.getString(allOpts[j]);
      }
      outf << endl;
      outf.close();
      dblf << dbFile << endl;
      toClean.push_back(base + ".sh*"); toClean.push_back(base + ".db");
    }
    dblf.close();
    fstream fin; MstUtils::openFile(fin, "fin." + MstSys::pathBase(op.getString("o")) + ".sh", ios::out);
    fin << op.getExecName() << " --dL " << dbListFile << " --o " << op.getString("o");
    if (op.isGiven("sim")) fin << " --sim " << op.getInt("sim");
    fin << endl;
    fin << "if [ $? -eq 0 ]; then # only clean up if database creation was successful" << endl;
    for (int k = 0; k < toClean.size(); k++) {
      fin << "  rm " << toClean[k] << endl;
    }
    fin << "fi" << endl;
    fin.close();
  }
  return 0; // success
}
