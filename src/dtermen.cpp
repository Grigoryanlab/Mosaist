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
  map<string, string> standard = {{"HSD", "HIS"}, {"HSE", "HIS"}, {"HSC", "HIS"}, {"HSP", "HIS"}};

  /* Almost perfectly corresponding to standard residues:
   * MSE -- selenomethyonine; SEC -- selenocysteine */
  map<string, string> almostStandard = {{"MSE", "MET"}, {"SEC", "CYS"}};

  /* A little less perfectly corresponding pairings, but can be acceptable (depends):
   * HIP -- ND1-phosphohistidine; SEP -- phosphoserine; TPO -- phosphothreonine;
   * PTR -- o-phosphotyrosine. */
  map<string, string> lessStandard = {{"HIP", "HIS"}, {"SEP", "SER"}, {"TPO", "THR"}, {"PTR", "TYR"}};

  for (res_t aai = 0; aai <= SeqTools::maxIndex(); aai++) {
    string aa = SeqTools::idxToTriple(aai);
    if (aai <= 19) {
      aaMap[aai] = aai; // natural amino acids go as they are, of course
    } else {
      switch (aaMapType) {
        case 1:
          // map frequent modifications to their closest standard amino acids
          // discard amino acids not explicitly mentioned in the above sets
          if (standard.find(aa) != standard.end()) {
            aaMap[aai] = SeqTools::aaToIdx(standard[aa]);
          } else if (almostStandard.find(aa) != almostStandard.end()) {
            aaMap[aai] = SeqTools::aaToIdx(almostStandard[aa]);
          } else if (lessStandard.find(aa) != lessStandard.end()) {
            aaMap[aai] = SeqTools::aaToIdx(lessStandard[aa]);
          }
          break;
        case 2:
          // map only the obvious modifications to their closest standard amino acid,
          // but mostly preserve the modifications (good for designing with modifications)
          // discard amino acids not explicitly mentioned in the above sets
          if (standard.find(aa) != standard.end()) {
            aaMap[aai] = SeqTools::aaToIdx(standard[aa]);
          } else if (almostStandard.find(aa) != almostStandard.end()) {
            aaMap[aai] = SeqTools::aaToIdx(almostStandard[aa]);
          } else if (lessStandard.find(aa) != lessStandard.end()) {
            aaMap[aai] = aai;
          }
          break;
        case 3:
          // do no chemical mapping, interpret everything explicitly, BUT still
          // discard amino acids not explicitly mentioned in the above sets
          if (standard.find(aa) != standard.end()) {
            aaMap[aai] = SeqTools::aaToIdx(standard[aa]);
          } else if (almostStandard.find(aa) != almostStandard.end()) {
            aaMap[aai] = aai;
          } else if (lessStandard.find(aa) != lessStandard.end()) {
            aaMap[aai] = aai;
          }
          break;
        case 4:
          // do no chemical mapping, interpret everything explicitly
          if (standard.find(aa) != standard.end()) {
            aaMap[aai] = SeqTools::aaToIdx(standard[aa]);
          } else {
            aaMap[aai] = aai;
          }
          break;
        default:
          MstUtils::error("unrecognized amino-acid mapping type " + MstUtils::toString(aaMapType));
      }
    }
  }
  for (res_t aai = 0; aai <= SeqTools::maxIndex(); aai++) {
    if (aaMap.find(aai) == aaMap.end()) continue; // not an allowed amino acid
    int i = aaMap[aai];
    if (aaIdx.find(i) != aaIdx.end()) { aaIdx[i] = aaIdx[i]; }
    else { aaIdx[i] = aaIdx.size(); }
  }
  globAlph.resize(aaIdx.size());
  for (auto it = aaIdx.begin(); it != aaIdx.end(); ++it) {
    globAlph[it->second] = it->first;
  }
}

void dTERMen::printAminoAcidMap() {
  cout << "---- amino-acid map ----" << endl;
  cout << "universal alphabet has " << globAlph.size() << " amino acids:\n\t";
  for (int i = 0; i < globAlph.size(); i++) cout << SeqTools::idxToTriple(globAlph[i]) << " (" << aaIdx[globAlph[i]] << "); ";
  cout << endl << "allowed mappings are:" << endl;
  for (auto it = aaMap.begin(); it != aaMap.end(); ++it) {
    cout << SeqTools::idxToTriple(it->first) << " -> " << SeqTools::idxToTriple(it->second) << endl;
  }
}

bool dTERMen::isInGlobalAlphabet(const string& aa) const {
  return (aaToIndex(aa) >= 0);
}

bool dTERMen::isInGlobalAlphabet(res_t aa) const {
  return (aaToIndex(aa) >= 0);
}

int dTERMen::aaToIndex(const string& aa) const {
  return aaToIndex(SeqTools::aaToIdx(aa));
}

int dTERMen::aaToIndex(res_t aa) const {
  if (aaIdx.find(aa) != aaIdx.end()) return aaIdx.at(aa);
  return -1;
}

string dTERMen::indexToResName(int idx) const {
  return SeqTools::idxToTriple(globAlph[idx]);
}

res_t dTERMen::indexToAA(int idx) const {
  return globAlph[idx];
}

void dTERMen::buildBackgroundPotentials() {
  // extract all necessary residue properties
  vector<string> propNames = {"phi", "psi", "omega", "env"};
  map<string, vector<mstreal> > propVals;
  vector<int> aa;
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
      if (!isInGlobalAlphabet(aaName)) continue;
      aa.push_back(aaToIndex(aaName));
      for (int i = 0; i < propNames.size(); i++) {
        propVals[propNames[i]].push_back(F.getResidueProperty(ti, propNames[i], ri));
      }
      propVals["mult"].push_back(mult[ri]);
    }
  }

  vector<vector<mstreal> > back(aa.size(), vector<mstreal> (globalAlphabetSize(), 0.0)); // existing background potential
  // TODO: implement potential lookup
  // TODO: implement 2D potential building!
  // ppPot = buildTwoDimPotential(propVals["phi"], propVals["psi"], {-180, 180, 72}, {-180, 180, 72}, aa, true, back);
  omPot = buildOneDimPotential(binData(propVals["omega"], 2, {-180, 180, 1000, 1}, propVals["mult"], true), aa, 1.0, back, true);
cout << "Omega potential:\n"; printOneDimPotential(omPot);
  envPot = buildOneDimPotential(binData(propVals["env"], 1, {0, 1, 75}, propVals["mult"], false), aa, 1.0, back, true);
cout << "Env potential:\n"; printOneDimPotential(envPot);
}

dTERMen::histType dTERMen::binData(const vector<mstreal>& X, int binSpecType, const vector<mstreal>& binSpec, const vector<mstreal>& M, bool isAngle) {
  // transform data if using dihedral angles
  vector<mstreal> X1, M1;
  vector<int> origIndices = MstUtils::range(0, (int) X.size());
  if (isAngle) {
    vector<int> exclude;
    for (int i = 0; i < X.size(); i++) {
      if (Residue::isBadDihedral(X[i])) {
        exclude.push_back(i);
        continue;
      }
      mstreal an = CartesianGeometry::angleDiff(X[i], 0);
      // represent +/- pi as -pi; then, thhe bin boundary definition always works
      if (an >= 180) an = -180;
      X1.push_back(an);
      M1.push_back(M.empty() ? 1.0 : M[i]);
    }
    origIndices = MstUtils::setdiff(origIndices, exclude);
  } else if (M.empty()) {
    for (int i = 0; i < X.size(); i++) M1.push_back(1.0);
  }
  const vector<mstreal>& x = (isAngle ? X1 : X);
  const vector<mstreal>& m = ((isAngle || M.empty()) ? M1 : M);
  vector<int> sortedInds = MstUtils::sortIndices(x);

  histType H;
  mstreal minVal, maxVal;
  if (binSpecType == 1) {
    /* Uniform binning */
    if (binSpec.size() != 3) MstUtils::error("expected three values for bin specification type 1", "dTERMen::binData");
    mstreal minVal = binSpec[0];
    mstreal maxVal = binSpec[1];
    int nb = (int) binSpec[2];
    if ((nb <= 0) || (minVal >= maxVal)) MstUtils::error("inconsistent bin specification, type 1", "dTERMen::binData");
    mstreal bw = (maxVal - minVal)/nb;
    H.binEdges.resize(nb + 1, 0);
    for (int i = 0; i < nb; i++) H.binEdges[i] = minVal + bw*i;
    H.binEdges[nb] = maxVal;
  } else if (binSpecType == 2) {
    /* Non-uniform binning with some minimal number of elements per bin and a
     * minimal bin width. */
    if (binSpec.size() != 4) MstUtils::error("expected four values for bin specification type 2", "dTERMen::binData");
    mstreal minVal = binSpec[0];
    mstreal maxVal = binSpec[1];
    int minNumPerBin = (int) binSpec[2];
    mstreal minBinWidth = binSpec[3];
    if ((minNumPerBin <= 0) || (minVal >= maxVal) || (minBinWidth < 0)) MstUtils::error("inconsistent bin specification, type 2", "dTERMen::binData");
    if (minNumPerBin > x.size()) MstUtils::error("requested min number of elements per bin exceeds the total number of elements", "dTERMen::binData");
    if (minBinWidth > maxVal - minVal) MstUtils::error("requested min bin width exceeds specified range", "dTERMen::binData");
    /* To account for multiplicities when counting bin sizes, create a vector of
     * cumulative data "mass", whereby massI[k] stores the mass for the first k
     * smallest elements (i.e., for indices in range sortedInds[0], sortedInds[1],
     * ..., sortedInds[k]). massE[k] is the same, but it excludes the mass of
     * the final point sortedInds[k] (having both arrays simplifies some counting
     * later). The mass of each data point is 1 over its multiplicity (if it is
     * defined). Thus, if multiplicity is 1 (i.e., the point is "unique"), then
     * its mass is 1. But if multiplicity is 2, meaning there are two roughly
     * equivalent positions in the dataset, then each will together count as 1
     * (each counting as 1/2). If multiplicities are not given, then all counts
     * are set to 1, which means data counts are interpreted explicitly. */
    vector<mstreal> massI(x.size()), massE(x.size());
    if (m.size() != x.size()) MstUtils::error("number of points and multiplicities specified is not consistent!", "dTERMen::binData");
    massI[0] = (1/m[sortedInds[0]]); massE[0] = 0;
    for (int i = 1; i < m.size(); i++) {
      massI[i] = massI[i-1] + 1/m[sortedInds[i]];
      massE[i] = massE[i-1] + 1/m[sortedInds[i-1]];
    }

    vector<mstreal> binEdgesFromLeft, binEdgesFromRight;
    binEdgesFromLeft.push_back(minVal);
    binEdgesFromRight.push_back(maxVal);
    // leftInd and rightInd will always store the first index (from left or right,
    // respectively) that maps into the next bin. Thus, by setting rightIndex to
    // the last index, the right edge of the last bin will be inclusive
    int leftInd = 0, rightInd = x.size() - 1;
    while (true) {
      int numRemaining = massI[rightInd] - massE[leftInd];
      mstreal remWidth = x[sortedInds[rightInd]] - x[sortedInds[leftInd]];
      if ((numRemaining >= 2*minNumPerBin) && (remWidth >= 2*minBinWidth)) { // if enough points left for at least two bins
        // add a bin on the left
        for (int k = leftInd; k < sortedInds.size(); k++) {
          if ((massE[k] - massI[leftInd] >= minNumPerBin) && (x[sortedInds[k]] - binEdgesFromLeft.back() > minBinWidth)) {
            binEdgesFromLeft.push_back(x[sortedInds[k]]);
            leftInd = k; // the left-most element in the next bin (left-to-right)
            break;
          }
        }

        // add a bin on the right
        for (int k = rightInd; k >= 0; k--) {
          if ((massI[rightInd] - massI[k] >= minNumPerBin) && (binEdgesFromRight.back() - x[sortedInds[k]] > minBinWidth)) {
            binEdgesFromRight.push_back(x[sortedInds[k]]);
            rightInd = k - 1; // the right-most element in the next bin (right-to-left)
            break;
          }
        }

        // if there are enough points left for no more than two bins, just split
        // the rest between the two bins
        if ((numRemaining < 3*minNumPerBin) || (remWidth < 3*minBinWidth)) {
          if (rightInd < leftInd) {
            /* Even though in this iteration there were initially enough points to fit three bins
            * (by width and number), it is possible that to add a left and right bin (respecting
            * minimal width and number of points) can move the terminal indices past each other.
            * Example, sorted values between leftInd and rightInd are: [0.5, 0.51, 0.52, 0.53,
            * 0.54, 0.55, 0.6, 0.7] minimum bin width is 0.1, and minimum number of points per bin
            * is 3. Then, to satisfy all conditions going from the left we need to do [0.5 - 0.6),
            * but to satisfy all conditions going from the right would require [0.54 - 0.7). The
            * two segments overlap. If this happens, we will make just one bin between leftInd and
            * rightInd, since it is not possible to make two while respecting all rules. */
            binEdgesFromLeft.pop_back();
            binEdgesFromRight.pop_back();
            H.binEdges = binEdgesFromLeft;
            H.binEdges.insert(H.binEdges.end(), binEdgesFromRight.rbegin(), binEdgesFromRight.rend());
          } else {
            int k = (leftInd + rightInd)/2;
            H.binEdges = binEdgesFromLeft;
            H.binEdges.back() = x[sortedInds[k]];
            H.binEdges.insert(H.binEdges.end(), binEdgesFromRight.rbegin() + 1, binEdgesFromRight.rend());
          }
          break;
        }
      } else if ((numRemaining >= minNumPerBin) && (remWidth >= minBinWidth)) { // if only enough points left for one bin
        // add the last bin to bridge the gap
        H.binEdges = binEdgesFromLeft;
        H.binEdges.insert(H.binEdges.end(), binEdgesFromRight.rbegin(), binEdgesFromRight.rend());
        break;
      } else {
        MstUtils::error("this should not have happened! Ran out of points/width for the middle bin.", "dTERMen::binData");
      }
    }
  }

  // classify points into bins
  H.bins.resize(H.binEdges.size() - 1);
  H.binMasses.resize(H.binEdges.size() - 1);
  H.weights.resize(X.size(), 0.0);
  int bi = 0; mstreal binMass = 0;
  for (int i = 0; i < sortedInds.size(); i++) {
    // the last bin includes the right-most end of the range
    if ((x[sortedInds[i]] >= H.binEdges[bi+1]) && ((bi < H.bins.size() - 1) || (x[sortedInds[i]] > H.binEdges[bi+1]))) {
      H.binMasses[bi] = binMass;
      binMass = 0;
      bi++;
    }
    H.bins[bi].push_back(origIndices[sortedInds[i]]);
    binMass += 1/m[sortedInds[i]];
    H.weights[origIndices[sortedInds[i]]] = 1/m[sortedInds[i]];
  }

  return H;
}

dTERMen::oneDimPotType dTERMen::buildOneDimPotential(const histType& H, const vector<int>& AA, mstreal pc, vector<vector<mstreal> >& backPot, bool updateBackPot) {
  oneDimPotType pot;
  pot.binEdges = H.binEdges;
  int nb = pot.binEdges.size() - 1;
  int naa = globalAlphabetSize();
  pot.aaEnergies.resize(nb, vector<mstreal> (naa, 0.0));

  // overall amino-acid frequencies
  CartesianPoint fAA(naa);
  for (int bi = 0; bi < nb; bi++) { // iterating over bin indices ignores any points excluded from binning
    const vector<int>& binInds = H.bins[bi];
    for (int i = 0; i < binInds.size(); i++) {
      int k = binInds[i];
      fAA[AA[k]] += H.weights[k];
    }
  }
  fAA /= fAA.sum();

  // add up the weighted counts of each amino acid in each bin
  for (int bi = 0; bi < nb; bi++) {
    vector<mstreal> expCounts = vector<mstreal> (naa, 0.0); // expectation of each amino acid in this bin
    const vector<int>& binInds = H.bins[bi];
    for (int i = 0; i < binInds.size(); i++) {
      int k = binInds[i];
      mstreal w = H.weights[k]; // position weight
      int aa = AA[k];
      pot.aaEnergies[bi][aa] += w;

      // compute the expectation of each amino acid in this position
      if (!backPot.empty()) {
        mstreal Z = 0; // partition function for this position (over all amino acids)
        for (int aai = 0; aai < naa; aai++) Z += exp(-backPot[k][aai]);
        for (int aai = 0; aai < naa; aai++) expCounts[aai] += w*exp(-backPot[k][aai])/Z;
      } else {
        // if no prior potential given, assume uniform prior
        for (int aai = 0; aai < naa; aai++) expCounts[aai] += w/naa;
      }
    }

    // compute the potential using frequency-proportional pseudocounts
    for (int aa = 0; aa < naa; aa++) {
      pot.aaEnergies[bi][aa] = -log((pot.aaEnergies[bi][aa] + pc*fAA[aa]*naa)/(expCounts[aa] + pc*fAA[aa]*naa));
    }
  }

  // update the prior potential
  if (!backPot.empty() && updateBackPot) {
    for (int bi = 0; bi < nb; bi++) {
      const vector<int>& binInds = H.bins[bi];
      for (int i = 0; i < binInds.size(); i++) {
        int k = binInds[i];
        for (int aai = 0; aai < naa; aai++) {
          backPot[k][aai] += pot.aaEnergies[bi][aai];
        }
      }
    }
  }

  return pot;
}

dTERMen::oneDimPotType dTERMen::buildOneDimPotential(const histType& H, const vector<int>& AA, mstreal pc) {
  vector<vector<mstreal> > backPot;
  return buildOneDimPotential(H, AA, pc, backPot, false);
}

void dTERMen::printOneDimPotential(const oneDimPotType& P) {
  int nb = P.aaEnergies.size();
  if (nb == 0) return;
  int naa = P.aaEnergies[0].size();
  for (int aa = 0; aa < naa; aa++) {
    cout << indexToResName(aa) << " ";
  }
  cout << endl;
  for (int bi = 0; bi < nb; bi++) {
    printf("%f %f", P.binEdges[bi], P.binEdges[bi+1]);
    for (int aa = 0; aa < naa; aa++) {
      printf(" %5.3f", P.aaEnergies[bi][aa]);
    }
    printf("\n");
  }
}

dTERMen::twoDimPotType dTERMen::buildTwoDimPotential(const vector<mstreal>& x, const vector<mstreal>& y, const vector<mstreal>& xBinSpec, const vector<mstreal>& yBinSpec, const vector<int>& aa, bool isAngle, const vector<mstreal>& priorPot, const vector<mstreal>& mult) {
  twoDimPotType pot;
  return pot;
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

mstreal EnergyTable::energyStdEst(int n) {
  if (n <= 0) MstUtils::error("sample size must be positive", "EnergyTable::energyStdEst");
  CartesianPoint energies(n, 0);
  for (int i = 0; i < n; i++) energies[i] = scoreSolution(randomSolution());
  return energies.stdev();
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
