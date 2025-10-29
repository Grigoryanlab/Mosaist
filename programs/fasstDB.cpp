#include "msttypes.h"
#include "mstoptions.h"
#include "mstfasst.h"
#include "mstrotlib.h"
#include "mstsystem.h"
#include "mstcondeg.h"
#include "mstsequence.h"
#include "mstexternal.h"
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
  op.addOption("s", "split final PDB files into chains by connectivity. Among other things, this avoids \"gaps\" within chains (where missing residues would go), which may simplify redundancy identification.");
  op.addOption("pp", "store phi/psi/omega properties in the database.");
  op.addOption("env", "store residue freedom property in the database. If this is give, --rLib must also be given.");
  op.addOption("cont", "store inter-residue contact information (for all residue pairs with contact degrees above the specified threshold). If this is given, --rLib must also be given.");
  op.addOption("contSeq", "store a specific version of inter-residue contact information which is calculated with amino acid constraints (for all residue pairs with contact degrees above the specified threshold). If this is given, --rLib must also be given.");
  op.addOption("int", "store residue to backbone contact information (for all residue pairs with contact degrees above the specified threshold). If this is given, --rLib must also be given.");
  op.addOption("bb", "store the minimum distance between backbone atoms of two residues (for all residue pairs below the specified cutoff).");
  op.addOption("stride", "store residue secondary structure classifications computed by STRIDE (external program). Argument must be the path to a STRIDE binary file.");
  op.addOption("sim", "percent sequence identity cutoff. If specified, will store local-window sequence similarity between all pairs of positions in the database, using this cutoff.");
  op.addOption("win", "window size to use with the similarity searching with --sim; must be an odd integer. Default is 31 (i.e., +/- 15 from the residue in question).");
  op.addOption("rLib", "path to an MST rotamer library file.");
  op.addOption("batch", "an integer. If specified, instead of building the database will spread all the work across this many "
                        "jobs, writing corresponding batch files for submission to the cluster. Will also produce a file called "
                        "<out>.fin.sh (where <out> is the base of the name specified in --o), which is to be run after all jobs "
                        "finish to complete the database building process.");
  op.addOption("slurm", "provide this option along with the batch argument to generate batch job files for a SLURM system");

  op.setOptions(argc, argv);
  RotamerLibrary RL;
  if (!op.isGiven("pL") && !op.isGiven("dL") && !op.isGiven("db")) MstUtils::error("either --pL, --dL, or --db must be given!");
  if (op.isGiven("rLib")) MstUtils::assertCond(MstSys::fileExists(op.getString("rLib")), "--rLib is not a valid file path");
  if (op.isGiven("cont") && (!op.isReal("cont") || (op.getReal("cont") < 0) || (op.getReal("cont") > 1))) MstUtils::error("--cont must be a real in range [0; 1]");
  if (op.isGiven("contSeq") && (!op.isReal("contSeq") || (op.getReal("contSeq") < 0) || (op.getReal("contSeq") > 1))) MstUtils::error("--cont must be a real in range [0; 1]");
  if (op.isGiven("int") && (!op.isReal("int") || (op.getReal("int") < 0) || (op.getReal("int") > 1))) MstUtils::error("--int must be a real in range [0; 1]");
  if (op.isGiven("bb") && (!op.isReal("bb") || (op.getReal("bb") < 0))) MstUtils::error("--bb must be a real greater than 0.");
  if ((op.isGiven("env") || op.isGiven("cont") || op.isGiven("int") || op.isGiven("bb")) || (op.isGiven("contSeq"))) {
    if (!op.isGiven("rLib")) MstUtils::error("--rLib is needed, but not given");
    RL.readRotamerLibrary(op.getString("rLib"));
  }
  if (op.isGiven("sim") && (!op.isReal("sim") || (op.getReal("sim") < 0) || (op.getReal("sim") > 100))) MstUtils::error("--sim must be a non-negative value below 100.");
  if (op.isGiven("win") && (!op.isInt("win") || (op.getInt("win") <= 0) || (op.getInt("win") % 2 == 0))) MstUtils::error("--win must be a positive odd integer.");
  short memSave = op.isGiven("m");
  
  // create map of prop names for contact degree with amino acid constraints
  vector<string> aaNames = {"ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU",
    "LYS", "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL", "ALA"}; //aa allowed in rotamer library
  map<string, string> aaToProp;
  for (string aa : aaNames) aaToProp[aa] = "cont"+aa;

  if (!op.isGiven("batch")) {
    FASST S;
    cout << "Reading structures..." << endl;
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
        if (op.isGiven("s")) {
          P = P.reassignChainsByConnectivity();
          P.deleteShortChains();
        }
        if (P.residueSize() != 0) S.addTarget(P, memSave);
        else cout << "skipping " << pdbFiles[i] << " as it ends up having no residues..." << endl;
      }
    }
    if (op.isGiven("db")) {
      S.readDatabase(op.getString("db"), memSave);
      // If adding to an existing DB, print what properties already exist within it
      cout << "Print the current properties of the DB" << endl;
      cout << "Structures in DB: " << S.numTargets() << endl;
      
      vector<string> res_prop = {"phi","psi","omega","env"};
      // Residue properties
      for (string prop : res_prop) cout << "Residue property: " << prop << " = " << S.isResiduePropertyDefined(prop) << endl;
      
      cout << "Residue string property: " << "stride" << " = " << S.isResidueStringPropertyDefined("stride") << endl;
      
      // Residue pair properties
      vector<string> res_pair_prop = {"cont","interfering","interfered","bb"};
      for (string prop: res_pair_prop) cout << "Residue pair property: " << prop << " = " << S.isResiduePairPropertyPopulated(prop) << endl;
      for (auto it: aaToProp) cout << "Residue pair property: " << it.second << " = " << S.isResiduePairPropertyPopulated(it.second) << endl;
      
      // Residue relational properties
      cout << "Residue relationship property: sim = " << S.isResidueRelationshipPopulated("sim") << endl;
    }
    if (op.isGiven("dL")) {
      vector<string> dbFiles = MstUtils::fileToArray(op.getString("dL"));
      for (int i = 0; i < dbFiles.size(); i++) {
        S.readDatabase(dbFiles[i], memSave);
      }
    }
    if (op.isGiven("pp") || op.isGiven("env") || op.isGiven("cont") || op.isGiven("contSeq") || op.isGiven("int") || op.isGiven("bb") || op.isGiven("stride")) {
      cout << "Computing per-target residue properties..." << endl;
      // compute and add some properties
      for (int ti = 0; ti < S.numTargets(); ti++) {
        cout << "\ttarget " << ti+1 << "/" << S.numTargets() << "..." << endl;
        Structure P = S.getTargetCopy(ti);
        if (op.isGiven("pp")) {
          vector<Residue*> residues = P.getResidues();
          vector<mstreal> phi(residues.size()), psi(residues.size()), omega(residues.size());
          for (int ri = 0; ri < residues.size(); ri++) {
            phi[ri] = residues[ri]->getPhi(false);
            psi[ri] = residues[ri]->getPsi(false);
            omega[ri] = residues[ri]->getOmega(false);
          }
          S.addResidueProperties(ti, "phi", phi);
          S.addResidueProperties(ti, "psi", psi);
          S.addResidueProperties(ti, "omega", omega);
        }
        if (op.isGiven("stride")) {
          string strideBin = op.getString("stride","");
          strideInterface stride(strideBin,&P);
          stride.computeSTRIDEClassifications();
          vector<string> strideSSType = stride.getSTRIDEClassifications();
          S.addResidueStringProperties(ti, "stride", strideSSType);
        }
        if (op.isGiven("env") || op.isGiven("cont") || op.isGiven("contSeq") || op.isGiven("int") || op.isGiven("bb")) {
          ConFind C(&RL, P); // both need the confind object
          // environment
          if (op.isGiven("env")) {
            vector<Residue*> residues = P.getResidues();
            vector<mstreal> freedoms = C.getFreedom(residues);
            S.addResidueProperties(ti, "env", freedoms);
          }
          // contact degree
          if (op.isGiven("cont")) {
            mstreal cdcut = op.getReal("cont");
            contactList list = C.getContacts(P, cdcut);
            map<int, map<int, mstreal> > conts;
            for (int i = 0; i < list.size(); i++) {
              int rA = list.residueA(i)->getResidueIndex();
              int rB = list.residueB(i)->getResidueIndex();
              conts[rA][rB] = list.degree(i);
              conts[rB][rA] = list.degree(i);
            }
            S.addResiduePairProperties(ti, "cont", conts);
          }
          // contact degree, with amino acid constraints
          /* Contact degree is calculated between residues i and j, with the rotamers at position i
           restricted to those of a specific amino acid. Note that this form of contact degree is no
           longer directional. To simplify access, the values for each potential amino acid at residue
           i are stored as distinct pair properties, e.g.: contARG, contASP, ..., contVAL.
           */
          if (op.isGiven("contSeq")) {
            mstreal cdcut = op.getReal("contSeq");
            contactList list = C.getConstrainedContacts(P.getResidues(), cdcut);
            
            set<string> aaNames = C.getAANames();
            for (string aa : aaNames) {
              map<int, map<int, mstreal> > conts;
              for (int i = 0; i < list.size(); i++) {
                const set<string>& alphaA = list.alphabetA(i);
                // skip if the A alphabet does not match current aa
                if ((alphaA.size() != 1) || *alphaA.begin() != aa) continue;
                int rA = list.residueA(i)->getResidueIndex();
                int rB = list.residueB(i)->getResidueIndex();
                conts[rA][rB] = list.degree(i);
              }
              S.addResiduePairProperties(ti, aaToProp[aa], conts);
            }
          }
          // interference
          /* To account for the directionality of interference, the property is stored in two maps:
           interfer_ing_ and interfer_ed_. These are named based on the map that is returned when the
           outer map is called with some residue, so calling interfering[resIdx] returns a map with
           residues whose backbones interfere with the sidechain of res and interfered[resIdx] returns a
           map with residues whose sidechains are interfered by the backbone of res. Note that while
           these store the exact same info, they simplify access.
           */
            if (op.isGiven("int")) {
                mstreal incut = op.getReal("int");
                contactList list = C.getInterference(P, incut);
                map<int, map<int, mstreal> > interfering;
                map<int, map<int, mstreal> > interfered;
                for (int i = 0; i < list.size(); i++) {
                    // rB backbone interferes with rA sidechain
                    int rA = list.residueA(i)->getResidueIndex();
                    int rB = list.residueB(i)->getResidueIndex();
                    interfering[rA][rB] = list.degree(i);
                    interfered[rB][rA] = list.degree(i);
                }
                S.addResiduePairProperties(ti, "interfering", interfering);
                S.addResiduePairProperties(ti, "interfered", interfered);
            }
          // backbone-backbone interaction
          if (op.isGiven("bb")) {
              mstreal dcut = op.getReal("bb");
              contactList list = C.getBBInteraction(P,dcut);
              map<int, map<int, mstreal>> bbInteraction;
              for (int i = 0; i < list.size(); i++) {
                  int rA = list.residueA(i)->getResidueIndex();
                  int rB = list.residueB(i)->getResidueIndex();
                  bbInteraction[rA][rB] = list.degree(i);
                  bbInteraction[rB][rA] = list.degree(i);
              }
              S.addResiduePairProperties(ti, "bb", bbInteraction);
          }
        }
      }
    }
    if (op.isGiven("sim")) {
      cout << "Computing local-window sequence similarity..." << endl;
      S.dropResidueRelationship("sim"); // overwrite any previously populated similarity property (could have been calculated differently)
      // first cluster all local windows
      cout << "\tgathering local sequence windows..." << endl;
      int Nr = 0;
      int L = op.getInt("win", 31);
      int L2 = (L - 1)/2;
      vector<Sequence> wins;
      vector<int> winTarg, winStart;
      vector<bool> nTerm, cTerm;
      for (int ti = 0; ti < S.numTargets(); ti++) {
        Structure P = S.getTargetCopy(ti);
        int off = 0;
        for (int i = 0; i < P.chainSize(); i++) {
          Chain& C = P[i];
          Nr += C.residueSize();
          for (int j = 0; j < C.residueSize() - L + 1; j++) {
            Sequence win(L);
            for (int k = 0; k < L; k++) win[k] = SeqTools::aaToIdx(C[j + k].getName());
            wins.push_back(win);
            winTarg.push_back(ti); winStart.push_back(off + j);
            nTerm.push_back(j == 0); cTerm.push_back(j == C.residueSize() - L);
          }
          off += C.residueSize();
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
// cout << "target " << ti << ", residue " << ri << " is similar to target " << tj << " residue " << rj << endl;
// cout << "\t" << wins[i].toString() << endl << "\t" << wins[j].toString() << endl << endl;
          S.addResidueRelationship(ti, "sim", ri, tj, rj); // adding in only one direction, because the list of clusters is bi-directional
          symN++;
          // if either window starts at the N-terminus, then the whole first half is redundant
          if (nTerm[i] || nTerm[j]) {
            for (int d = -L2; d < 0; d++) S.addResidueRelationship(ti, "sim", ri + d, tj, rj + d);
            symN += L2;
          }
          // if either window ends at the C-terminus, then the whole second half is redundant
          if (cTerm[i] || cTerm[j]) {
            for (int d = 1; d <= L/2; d++) S.addResidueRelationship(ti, "sim", ri + d, tj, rj + d);
            symN += L2;
          }
        }
      }
      cout << "\trecorded " << symN << " similar windows, from a total of " << Nr << " residues" << endl;
    }
    S.writeDatabase(op.getString("o"));
  } else {
    if (!op.isGiven("pL")) MstUtils::error("--pL must be given with --batch");
    if (!op.isInt("batch") || (op.getInt("batch") <= 0)) MstUtils::error("--batch must be a positive integer!");
    int nJ = op.getInt("batch");
    vector<string> pdbFiles = MstUtils::fileToArray(op.getString("pL"));
    MstUtils::seedRandEngine();
    MstUtils::shuffle(pdbFiles);
    vector<pair<int, int> > tasks = MstUtils::splitTasks(pdbFiles.size(), nJ);
    fstream outf;
    string dbListFile = MstSys::pathBase(op.getString("o")) + ".dL";
    fstream dblf; MstUtils::openFile(dblf, dbListFile, ios::out);
    vector<string> toClean; toClean.push_back(dbListFile);
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
      // time per structure depends on what properties need to be calculated
      int time_per_structure = 5;
      vector<string> allOpts = op.getAllGivenOptions();
      for (int j = 0; j < allOpts.size(); j++) {
        if ((allOpts[j].compare("cont") == 0) || (allOpts[j].compare("contSeq") == 0) || (allOpts[j].compare("int") == 0) || (allOpts[j].compare("stride") == 0)) time_per_structure += 5 ;
      }
      cout << "Given options, calculating that " << time_per_structure << " min are needed per structure (on average)" << endl;
      int hrs = (int) ceil(time_per_structure*(tasks[i].second - tasks[i].first + 1)/60.0); // fifteen minutes per structure should be plenty
      outf << "#!/bin/bash\n";
      if (op.isGiven("slurm")) {
        outf << "#SBATCH -J fasstDB.%A\n" << "#SBATCH -o fasstDB.%A.log\n";
        outf << "#SBATCH -p defq\n" << "#SBATCH -n 1\n" << "#SBATCH --mem=2G\n";
        outf << "#SBATCH -t 0-" << hrs << ":00:00\n";
      } else {
        outf << "#$ -j y\n" << "#$ -cwd\n" << "#$ -V\n";
        outf << "#$ -l vf=2G\n" << "#$ -l ironfs\n";
        outf << "#$ -l h_rt=" << hrs << ":00:00\n";
      }
      outf << op.getExecName() << " --pL " << listFile << " --o " << dbFile;
      // keep all other options from the call to self
      for (int j = 0; j < allOpts.size(); j++) {
        if ((allOpts[j].compare("batch") == 0) || (allOpts[j].compare("pL") == 0) || (allOpts[j].compare("o") == 0) || (allOpts[j].compare("sim") == 0)) continue;
        outf << " --" << allOpts[j] << " " << op.getString(allOpts[j]);
      }
      outf << endl;
      outf.close();
      dblf << dbFile << endl;
      toClean.push_back(base + ".sh*"); toClean.push_back(base + ".db"); toClean.push_back(base + ".list");
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
  cout << "Done" << endl; // this makes it easy to check if all the jobs were completed
  return 0; // success
}
