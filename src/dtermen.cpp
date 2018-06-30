#include "dtermen.h"

dTERMen::dTERMen(const string& configFile) {
  init(configFile);
  buildBackgroundPotentials();
}

void dTERMen::init(const string& configFile) {
  vector<string> lines = MstUtils::fileToArray(configFile);
  for (int i = 0; i < lines.size(); i++) {
    string line = MstUtils::trim(MstUtils::removeComment(lines[i], "#"));
    if (line.empty()) continue;
    vector<string> ents = MstUtils::trim(MstUtils::split(line, "="));
    if (ents.size() != 2) MstUtils::error("could not parse parameter line '" + lines[i] + "' from file " + configFile, "dTERMen::dTERMen(const string&)");
    if (ents[0].compare("fasstdb") == 0) {
      fasstdbPath = ents[1];
    } else {
      MstUtils::error("unknown parameter name '" + ents[0] + "'", "dTERMen::dTERMen(const string&)");
    }
  }

  if (fasstdbPath.empty()) MstUtils::error("FASST database not defined in configuration file " + configFile, "dTERMen::dTERMen(const string&)");
  F.readDatabase(fasstdbPath);
}

void dTERMen::buildBackgroundPotentials() {
  vector<string> propNames = {"phi", "psi", "omega", "env"};
  map<string, vector<mstreal> > propVals;
  for (int i = 0; i < propNames.size(); i++) MstUtils::assert(F.isResiduePropertyDefined(propNames[i]), "property " + propNames[i] + " is not defined in the FASST database", "dTERMen::buildBackgroundPotentials()");

  for (int ti = 0; ti < F.numTargets(); ti++) {
    int N = F.getTarget(ti)->residueSize();
    for (int ri = 0; ri < N; ri++) {
      for (int i = 0; i < propNames.size(); i++) {
        propVals[propNames[i]].push_back(F.getResidueProperty(ti, propNames[i], ri));
      }
    }
    map<int, map<int, set<int> > > simsToTarget = F.getResidueRelationships(ti, "sim");
    vector<mstreal> vals(N, 1); // multiplicity of each residue in the database
    for (auto it = simsToTarget.begin(); it != simsToTarget.end(); ++it) {
      map<int, set<int> >& simsInTarget = it->second;
      for (auto jt = simsInTarget.begin(); jt != simsInTarget.end(); ++jt) {
        int ri = jt->first;
        vals[ri] += (jt->second).size();
      }
    }
    propVals["mult"].insert(propVals["mult"].end(), vals.begin(), vals.end());
  }
}


/* ----------- EnergyTable --------------- */
EnergyTable::EnergyTable(const string& tabFile) {
  readFromFile(tabFile);
}

void EnergyTable::readFromFile(const string& tabFile) {
  vector<string> lines = MstUtils::fileToArray(tabFile);
  int i = 0, s, a;

  /* read self energies and site addresses */
  for (; i < lines.size(); i++) {
    vector<string> line = MstUtils::split(lines[i]);
    if (line.size() != 3) break;
    if (siteIndices.find(line[0]) == siteIndices.end()) {
      s = siteIndices.size();
      siteIndices[line[0]] = s;
      aaIndices.push_back(map<string, int>());
      selfE.push_back(vector<mstreal>());
      indToAA.push_back(map<int, string>());
    } else { s = siteIndices[line[0]]; }
    a = aaIndices[s].size();
    aaIndices[s][line[1]] = a;
    indToAA[s][a] = line[1];
    selfE[s].push_back(MstUtils::toReal(line[2]));
  }

  /* read pair energies */
  pairs.clear(); pairs.resize(siteIndices.size());
  pairMaps.clear(); pairMaps.resize(siteIndices.size());
  pairE.clear(); pairE.resize(siteIndices.size());
  for (; i < lines.size(); i++) {
    vector<string> line = MstUtils::split(lines[i]);
    if (line.size() != 5) MstUtils::error("could not parse line '" + lines[i] + "'");
    if ((siteIndices.find(line[0]) == siteIndices.end()) || (siteIndices.find(line[1]) == siteIndices.end())) {
      MstUtils::error("unexpected site names in pair line: " + lines[i]);
    }
    int si = siteIndices[line[0]];
    int sj = siteIndices[line[1]];
    if ((aaIndices[si].find(line[2]) == aaIndices[si].end()) || (aaIndices[sj].find(line[3]) == aaIndices[sj].end())) {
      MstUtils::error("unexpected amino-acid names in pair line: " + lines[i]);
    }
    int aai = aaIndices[si][line[2]];
    int aaj = aaIndices[si][line[3]];
    pairMaps[si][sj] = true; pairMaps[sj][si] = true;
    // the first time we encounter a site pair, fill them up with all-by-all amino-acid zero energies
    if (pairE[si].size() == 0) pairE[si].resize(siteIndices.size());
    if (pairE[si][sj].size() == 0) {
      pairE[si][sj].resize(aaIndices[si].size(), vector<mstreal>(aaIndices[sj].size(), 0));
    }
    pairE[si][sj][aai][aaj] = MstUtils::toReal(line[4]);
  }
  for (i = 0; i < pairMaps.size(); i++) pairs[i] = MstUtils::keys(pairMaps[i]);
}

mstreal EnergyTable::meanEnergy() const {
  mstreal mE = 0;
  for (int i = 0; i < selfE.size(); i++) {
    mstreal m = 0;
    for (int aa = 0; aa < selfE[i].size(); aa++) m += selfE[i][aa];
    mE += m/selfE[i].size();
  }
  for (int i = 0; i < pairE.size(); i++) {
    for (int j = 0; j < pairE[i].size(); j++) {
      mstreal m = 0;
      int nij = 0;
      for (int aai = 0; aai < pairE[i][j].size(); aai++) {
        for (int aaj = 0; aaj < pairE[i][j][aai].size(); aaj++) {
          m += pairE[i][j][aai][aaj];
          nij++;
        }
      }
      if (nij != 0) mE += m/nij;
    }
  }
  return mE;
}

mstreal EnergyTable::selfEnergy(int s, int aa) {
  return selfE[s][aa];
}

mstreal EnergyTable::pairEnergy(int si, int sj, int aai, int aaj) {
  if (si < sj) return pairE[si][sj][aai][aaj];
  return pairE[sj][si][aaj][aai];
}

mstreal EnergyTable::scoreSolution(const vector<int>& sol) {
  if (sol.size() != selfE.size()) MstUtils::error("solution of wrong length for table", "EnergyTable::scoreSolution(const vector<int>&)");
  mstreal ener = 0;
  for (int si = 0; si < selfE.size(); si++) {
    ener += selfEnergy(si, sol[si]);
    const vector<int>& sj = pairs[si];
    for (int j = 0; j < sj.size(); j++) {
      if (sj[j] < si) continue; // do not double-count interactions
      ener += pairEnergy(si, sj[j], sol[si], sol[sj[j]]);
    }
  }
  return ener;
}

mstreal EnergyTable::scoreSequence(const Sequence& seq) {
  return scoreSolution(sequenceToSolution(seq));
}

mstreal EnergyTable::scoreMutation(const vector<int>& sol, int mutSite, int mutAA) {
  if (sol.size() != selfE.size()) MstUtils::error("wild-type solution of wrong length for table", "EnergyTable::scoreMutation(const vector<int>&, int, const string&)");
  if ((mutSite < 0) || (mutSite >= selfE.size())) MstUtils::error("mutation site index out of range for table", "EnergyTable::scoreMutation(const vector<int>&, int, const string&)");

  mstreal dE = selfE[mutSite][mutAA] - selfE[mutSite][sol[mutSite]];
  const vector<int>& intSites = pairs[mutSite];
  for (int j = 0; j < intSites.size(); j++) {
    int intSite = intSites[j];
    dE += pairEnergy(mutSite, intSite, mutAA, sol[intSite]) - pairEnergy(mutSite, intSite, sol[mutSite], sol[intSite]);
  }
  return dE;
}

mstreal EnergyTable::scoreMutation(const vector<int>& sol, const vector<int>& mutSites, const vector<int>& mutAAs) {
  vector<int> mutSol = sol;
  mstreal dE = 0;
  for (int i = 0; i < mutSites.size(); i++) {
    dE += scoreMutation(mutSol, mutSites[i], mutAAs[i]);
    mutSol[mutSites[i]] = mutAAs[i];
  }
  return dE;
}

vector<int> EnergyTable::randomSolution() const {
  vector<int> sol(selfE.size());
  for (int i = 0; i < selfE.size(); i++) {
    sol[i] = MstUtils::randInt(0, selfE[i].size() - 1);
  }
  return sol;
}

int EnergyTable::randomResidue(int si) const {
  return MstUtils::randInt(0, selfE[si].size() - 1);
}

vector<int> EnergyTable::mc(int Nc, int Ni, mstreal kTi, mstreal kTf, int annealType, void* rec, void (*add)(void*, const vector<int>&, mstreal)) {
  if (kTf < 0) kTf = kTi;
  mstreal kT = kTi;
  int L = numSites();
  mstreal bE = 0; vector<int> bS;
  for (int cyc = 0; cyc < Nc; cyc++) {
    // equilibration
    vector<int> seq = randomSolution();
    for (int i = 0; i < int(0.2*Ni + 1); i++) {
      int s = MstUtils::randInt(0, L - 1);
      int aa = randomResidue(s);
      mstreal dE = scoreMutation(seq, s, aa);
      if (MstUtils::randUnit() < exp(-dE/kT)) seq[s] = aa;
    }

    // data collection
    mstreal bcE = scoreSolution(seq); vector<int> bcS = seq;
    if (cyc == 0) { bE = bcE; bS = bcS; }
    mstreal dEb = 0; // current change from the best energy
    for (int i = 0; i < Ni; i++) {
      mstreal f = i*1.0/(Ni-1);
      // every now and again, calculate the total energy to avoid accumulation of addition errors
      if ((i+1) % 10000 == 0) {
        dEb = scoreSolution(seq) - bcE;
      }
      int s = MstUtils::randInt(0, L - 1);
      int aa = randomResidue(s);
      mstreal dE = scoreMutation(seq, s, aa);
      // annealing schedule
      switch (annealType) {
        case 1:
          // linear
          kT = kTf*f + kTi*(1-f);
          break;

        case 2:
          // exponential
          kT = kTi*pow(kTf/kTi, f);
          break;

        default:
          MstUtils::error("unrecognized annealing schedule type '" + MstUtils::toString(annealType) + "'");
      }

      if (MstUtils::randUnit() < exp(-dE/kT)) {
        seq[s] = aa;
        dEb += dE;
        if (dEb < 0) {
          bcE += dEb;
          bcS = seq;
          dEb = 0;
        }
      }
      if ((rec != NULL) && (add != NULL)) (*add)(rec, seq, bcE + dEb);
    }
    if ((cyc == 0) || (bE > bcE)) {
      bE = bcE;
      bS = bcS;
    }
  }

  return bS;
}

Sequence EnergyTable::solutionToSequence(const vector<int>& sol) {
  Sequence seq(sol.size());
  for (int i = 0; i < sol.size(); i++) {
    seq[i] = SeqTools::aaToIdx(getResidueString(i, sol[i]));
  }
  return seq;
}

vector<int> EnergyTable::sequenceToSolution(const Sequence& seq) {
  if (seq.size() != selfE.size()) MstUtils::error("sequence of wrong length for table", "EnergyTable::sequenceToSolution(const Sequence&)");
  vector<int> seqInts(seq.size());
  for (int i = 0; i < seq.size(); i++) {
    string aa3 = seq.getResidue(i, true);
    if (aaIndices[i].find(aa3) == aaIndices[i].end()) MstUtils::error("sequence is not from table alphabet", "EnergyTable::sequenceToSolution(const Sequence&)");
    seqInts[i] = aaIndices[i][aa3];
  }
  return seqInts;
}

string EnergyTable::getResidueString(int si, int ri) {
  return indToAA[si][ri];
}
