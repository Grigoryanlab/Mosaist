#include "dtermen.h"

dTERMen::dTERMen() {
  kT = 1.0;
  aaMapType = 1;
  setAminoAcidMap();
}

dTERMen::dTERMen(const string& configFile) {
  kT = 1.0;
  aaMapType = 1;
  setAminoAcidMap();
  readConfigFile(configFile);
}

void dTERMen::readConfigFile(const string& configFile) {
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
  if (!backPotFile.empty()) {
    readBackgroundPotentials(backPotFile);
  } else {
    buildBackgroundPotentials();
  }
}

void dTERMen::setAminoAcidMap() {
  /* Perfectly corresponding to standard residues. */
  map<string, string> standard = {{"HSD", "HIS"}, {"HSE", "HIS"}, {"HSC", "HIS"}, {"HSP", "HIS"}, {"MSE", "MET"}};

  /* Almost perfectly corresponding to standard residues:
   * MSE -- selenomethyonine; SEC -- selenocysteine */
  map<string, string> almostStandard = {{"MSE", "MET"}, {"SEC", "CYS"}};

  /* A little less perfectly corresponding pairings, but can be acceptable (depends):
   * HIP -- ND1-phosphohistidine; SEP -- phosphoserine; TPO -- phosphothreonine;
   * PTR -- o-phosphotyrosine. */
  map<string, string> lessStandard = {{"HIP", "HIS"}, {"SEP", "SER"}, {"TPO", "THR"}, {"PTR", "TYR"}};

  for (res_t aaIdx = 0; aaIdx <= SeqTools::maxIndex(); aaIdx++) {
    string aa = SeqTools::idxToTriple(aaIdx);
    switch (aaMapType) {
      if (aaIdx > 19) { // natural amino acids go as they are, of course (no mapping needed)
        case 1:
          // map frequent modifications to their closest standard amino acids
          if (standard.find(aa) != standard.end()) {
            aaMap[aa] = standard[aa];
          } else if (almostStandard.find(aa) != almostStandard.end()) {
            aaMap[aa] = almostStandard[aa];
          } else if (lessStandard.find(aa) != lessStandard.end()) {
            aaMap[aa] = lessStandard[aa];
          } else {
            aaMap[aa] = aa;
          }
          break;
        case 2:
          // map only the obvious modifications to their closest standard amino acid,
          // but mostly preserve the modifications (good for designing with modifications)
          if (standard.find(aa) != standard.end()) {
            aaMap[aa] = standard[aa];
          } else if (almostStandard.find(aa) != almostStandard.end()) {
            aaMap[aa] = almostStandard[aa];
          } else {
            aaMap[aa] = aa;
          }
          break;
        case 3:
          // do not do any chemical mapping, interpret everything explicitly
          if (standard.find(aa) != standard.end()) {
            aaMap[aa] = standard[aa];
          } else {
            aaMap[aa] = aa;
          }
          break;
        default:
          MstUtils::error("unrecognized amino-acid mapping type " + MstUtils::toString(aaMapType));
      }
    }
  }
}

res_t dTERMen::aaToIndex(const string& aa) const {
  if (aaMap.find(aa) != aaMap.end()) return SeqTools::aaToIdx(aaMap.at(aa)); // using at() instead of operator[] enables const
  return SeqTools::aaToIdx(aa);
}


void dTERMen::buildBackgroundPotentials() {
  // extract all necessary residue properties
  vector<string> propNames = {"phi", "psi", "omega", "env"};
  map<string, vector<mstreal> > propVals;
  vector<res_t> aa;
  for (int i = 0; i < propNames.size(); i++) MstUtils::assert(F.isResiduePropertyDefined(propNames[i]), "property " + propNames[i] + " is not defined in the FASST database", "dTERMen::buildBackgroundPotentials()");

  for (int ti = 0; ti < F.numTargets(); ti++) {
    Structure* S = F.getTarget(ti);
    int N = S->residueSize();

    // compute multiplicity of each residue
    vector<mstreal> mult(N, 1); // multiplicity of each residue in the structure
    map<int, map<int, set<int> > > simsToTarget = F.getResidueRelationships(ti, "sim");
    for (auto it = simsToTarget.begin(); it != simsToTarget.end(); ++it) {
      map<int, set<int> >& simsInTarget = it->second;
      for (auto jt = simsInTarget.begin(); jt != simsInTarget.end(); ++jt) {
        int ri = jt->first;
        mult[ri] += (jt->second).size();
      }
    }

    // store all properties
    for (int ri = 0; ri < N; ri++) {
      string aaName = (S->getResidue(ri)).getName();
      if (SeqTools::isUnknown(aaName)) continue;
      aa.push_back(aaToIndex(aaName));
      for (int i = 0; i < propNames.size(); i++) {
        propVals[propNames[i]].push_back(F.getResidueProperty(ti, propNames[i], ri));
      }
      propVals["mult"].push_back(mult[ri]);
    }
  }

  vector<mstreal> back(aa.size(), 0); // existing background potential
  buildTwoDimPotential(propVals["phi"], propVals["psi"], {-180, 180, 72}, {-180, 180, 72}, aa, true, back, propVals["mult"]);
  buildOneDimPotential(propVals["omega"], 2, {-180, 180, 1000, 1}, aa, true, back, propVals["mult"]);
  buildOneDimPotential(propVals["env"], 1, {0, 1, 75}, aa, false, back, propVals["mult"]);
}

void dTERMen::buildOneDimPotential(const vector<mstreal>& x, int binSpecType, const vector<mstreal>& binSpec, const vector<res_t>& aa, bool isAngle, const vector<mstreal>& priorPot, const vector<mstreal>& mult) {
  /* bin k (k \in [0; n-1], where n is the number of bins) is:
   * [binEdges[k]; binEdges[k+1]). Except for the last bin, which does include
   * the right limit: [binEdges[k]; binEdges[k+1]]. */
  vector<mstreal> binEdges;
  mstreal minVal, maxVal;
  int nb;
  if (binSpecType == 1) {
    /* Uniform binning */
    if (binSpec.size() != 3) MstUtils::error("expected three values for bin specification type 1", "dTERMen::buildOneDimPotential");
    mstreal minVal = binSpec[0];
    mstreal maxVal = binSpec[1];
    int nb = (int) binSpec[2];
    if ((nb <= 0) || (minVal >= maxVal)) MstUtils::error("inconsistent bin specification, type 1", "dTERMen::buildOneDimPotential");
    mstreal bw = (maxVal - minVal)/nb;
    binEdges.resize(nb + 1, 0);
    for (int i = 0; i < nb; i++) binEdges[i] = minVal + bw*i;
    binEdges[nb] = maxVal;

  } else if (binSpecType == 2) {
// TODO: account for multiplicities when counting even bins!!! To do so, normalize
// counts so that the total stays the same, but relative counts are inversly proportional
// to multiplicities, if they are given. If not, all counts are 1.0.
    /* Non-uniform binning with some minimal number of elements per bin and a
     * minimal bin width. */
    if (binSpec.size() != 4) MstUtils::error("expected four values for bin specification type 2", "dTERMen::buildOneDimPotential");
    mstreal minVal = binSpec[0];
    mstreal maxVal = binSpec[1];
    int minNumPerBin = (int) binSpec[2];
    mstreal minBinWidth = binSpec[3];
    if ((minNumPerBin <= 0) || (minVal >= maxVal) || (minBinWidth < 0)) MstUtils::error("inconsistent bin specification, type 2", "dTERMen::buildOneDimPotential");
    if (minNumPerBin > x.size()) MstUtils::error("requested min number of elements per bin exceeds the total number of elements", "dTERMen::buildOneDimPotential");
    if (minBinWidth > maxVal - minVal) MstUtils::error("requested min bin width exceeds specified range", "dTERMen::buildOneDimPotential");
    vector<int> sortedInds = MstUtils::sortIndices(x);
    vector<mstreal> binEdgesFromLeft, binEdgesFromRight;
    binEdgesFromLeft.push_back(minVal);
    binEdgesFromRight.push_back(maxVal); //
    // leftInd and rightInd will always store the first index (from left or right,
    // respectively) that maps into the next bin. Thus, by setting rightIndex to
    // the last index, the right edge of the last bin will be inclusive
    int leftInd = 0, rightInd = x.size() - 1;
    while (true) {
      int numRemaining = rightInd - leftInd + 1;
      mstreal remWidth = x[sortedInds[rightInd]] - x[sortedInds[leftInd]];
cout << "numRemaining = " << numRemaining << ", leftInd = " << leftInd << ", rightInd = " << rightInd << endl;
      if ((numRemaining >= 2*minNumPerBin) && (remWidth >= 2*minBinWidth)) { // if enough points left for at least two bins
        // add a bin on the left
        for (int k = leftInd + minNumPerBin - 1; k < sortedInds.size(); k++) {
          if (x[sortedInds[k]] - binEdgesFromLeft.back() > minBinWidth) {
            binEdgesFromLeft.push_back(x[sortedInds[k]]);
            leftInd = k; // the left-most element in the next bin (left-to-right)
            break;
          }
        }

        // add a bin on the right
        for (int k = rightInd - minNumPerBin + 1; k >= 0; k--) {
          if (binEdgesFromRight.back() - x[sortedInds[k]] > minBinWidth) {
            binEdgesFromRight.push_back(x[sortedInds[k + 1]]);
            rightInd = k; // the right-most element in the next bin (right-to-left)
            break;
          }
        }

        if ((numRemaining < 3*minNumPerBin) || (remWidth < 3*minBinWidth)) { // if there are enough points left for no more than two bins
          // split the remaining points in half between the two bins we just added
          int k = (leftInd + rightInd)/2;
cout << "splicing: " << endl << MstUtils::vecToString(binEdgesFromLeft) << endl << MstUtils::vecToString(binEdgesFromRight) << endl << "with k = " << x[sortedInds[k]] << endl;
          binEdges = binEdgesFromLeft;
          binEdges.back() = x[sortedInds[k]];
          binEdges.insert(binEdges.end(), binEdgesFromRight.rbegin() + 1, binEdgesFromRight.rend());
cout << "gives" << endl << MstUtils::vecToString(binEdges) << endl;
          break;
        }
      } else if ((numRemaining >= minNumPerBin) && (remWidth >= minBinWidth)) { // if only enough points left for one bin
        // add the last bin to bridge the gap
        binEdges = binEdgesFromLeft;
        binEdges.insert(binEdges.end(), binEdgesFromRight.rbegin(), binEdgesFromRight.rend());
        break;
      } else {
        MstUtils::error("this should not have happened! Ran out of points/width for the middle bin.", "dTERMen::buildOneDimPotential");
      }
    }
  }

cout << "binning of type " << binSpecType << " with parameters: [" << MstUtils::vecToString(binSpec) << "] produces edges: " << endl << MstUtils::vecToString(binEdges) << endl;
}

void dTERMen::buildTwoDimPotential(const vector<mstreal>& x, const vector<mstreal>& y, const vector<mstreal>& xBinSpec, const vector<mstreal>& yBinSpec, const vector<res_t>& aa, bool isAngle, const vector<mstreal>& priorPot, const vector<mstreal>& mult) {

}

void dTERMen::readBackgroundPotentials(const string& file) {

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

vector<int> EnergyTable::mc(int Nc, int Ni, mstreal kTi, mstreal kTf, int annealType, void* rec, void (*add)(void*, const vector<int>&, mstreal), int Ne) {
  if (kTf < 0) kTf = kTi;
  mstreal kT = kTi;
  int L = numSites();
  if (Ne < 0) Ne = int(0.2*Ni + 1);
  mstreal bE = 0; vector<int> bS;
  for (int cyc = 0; cyc < Nc; cyc++) {
    // equilibration
    vector<int> seq = randomSolution();
    for (int i = 0; i < Ne; i++) {
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
