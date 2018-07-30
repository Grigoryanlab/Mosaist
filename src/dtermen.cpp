#include "dtermen.h"

dTERMen::dTERMen() {
  init();
}

dTERMen::dTERMen(const string& configFile) {
  init();
  readConfigFile(configFile);
}

void dTERMen::init() {
  kT = 1.0;
  aaMapType = 1;
  cdCut = 0.01;
  pmSelf = 1;
  pmPair = 1;
  selfResidualPC = selfCorrPC = 1.0;
  selfResidualMinN = 1000;
  selfResidualMaxN = 5000;
  selfCorrMinN = 200;
  selfCorrMinN = 5000;
  setAminoAcidMap();
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
    } else if (ents[0].compare("rotlib") == 0) {
      rotLibFile = ents[1];
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
  if (rotLibFile.empty()) { MstUtils::error("dTERMen configuration file does not specify a rotamer library, '" + configFile + "'", "dTERMen::readConfigFile"); }
  else {
    RL.readRotamerLibrary(rotLibFile);
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

  // for each position in the database, accumulate total statistical potential,
  // for  all possible amino acids, from all known pseudo-energy types
  vector<vector<mstreal> > back(aa.size(), vector<mstreal> (globalAlphabetSize(), 0.0));
  bkPot = buildZeroDimPotential(aa, back);
  // cout << "Background frequency potential:\n"; printZeroDimPotential(bkPot);
  ppPot = buildTwoDimPotential(binData(propVals["phi"], propVals["psi"], {-180, 180, 36}, {-180, 180, 36}, propVals["mult"], true), aa, 10.0, back, true);
  // cout << "Phi/psi potential:\n"; printTwoDimPotential(ppPot);
  omPot = buildOneDimPotential(binData(propVals["omega"], 2, {-180, 180, 1000, 1}, propVals["mult"], true), aa, 10.0, back, true);
  // cout << "Omega potential:\n"; printOneDimPotential(omPot);
  envPot = buildOneDimPotential(binData(propVals["env"], 1, {0, 1, 75}, propVals["mult"], false), aa, 10.0, back, true);
  // cout << "Env potential:\n"; printOneDimPotential(envPot);
// TODO: could also do an end potential: 1, 2, 3 treated specially and N-2, N-1
// N also.
}

int dTERMen::findBin(const vector<mstreal>& binEdges, mstreal x) {
  if ((x < binEdges.front()) || (x > binEdges.back())) return -1;
  // do a binary search
  int nb = binEdges.size() - 1;
  int left = 0, right = nb - 1, k;
  while (true) {
    k = (left + right)/2;
    if (x < binEdges[k]) {
      right = k - 1;
    } else if ((x > binEdges[k+1]) || ((k < nb - 1) && (x == binEdges[k+1]))) { // the last bin includes the opper limit
      left = k + 1;
    } else {
      break;
    }
  }
  return k;
}

mstreal dTERMen::lookupZeroDimPotential(const zeroDimPotType& P, int aa) {
  return P.aaEnergies[aa];
}

mstreal dTERMen::lookupOneDimPotential(const oneDimPotType& P, mstreal x, int aa) {
  int k = findBin(P.binEdges, x);
  if (k < 0) return 0;
  return P.aaEnergies[k][aa];
}

mstreal dTERMen::lookupTwoDimPotential(const twoDimPotType& P, mstreal x, mstreal y, int aa) {
  int kx = findBin(P.xBinEdges, x);
  if (kx < 0) return 0;
  int ky = findBin(P.yBinEdges, y);
  if (kx < 0) return 0;
  return P.aaEnergies[kx][ky][aa];
}

dTERMen::twoDimHist dTERMen::binData(const vector<mstreal>& X, const vector<mstreal>& Y, const vector<mstreal>& xBinSpec, const vector<mstreal>& yBinSpec, const vector<mstreal>& M, bool isAngle) {
  dTERMen::oneDimHist xH = binData(X, 1, xBinSpec, M, isAngle);
  dTERMen::oneDimHist yH = binData(Y, 1, yBinSpec, M, isAngle);
  dTERMen::twoDimHist xyH;
  xyH.xBinEdges = xH.binEdges;
  xyH.yBinEdges = yH.binEdges;
  xyH.weights = xH.weights; // the weights vector should be the same for X- and Y- histograms, so can copy from either
  xyH.bins.resize(xH.bins.size(), vector<vector<int> >(yH.bins.size()));

  // find intersections between X and Y bins
  map<int, int> yPointsToBins;
  for (int i = 0; i < yH.bins.size(); i++) {
    for (int k = 0; k < yH.bins[i].size(); k++) {
      yPointsToBins[yH.bins[i][k]] = i;
    }
  }
  for (int i = 0; i < xH.bins.size(); i++) {
    for (int k = 0; k < xH.bins[i].size(); k++) {
      int idx = xH.bins[i][k];
      if (yPointsToBins.find(idx) != yPointsToBins.end()) {
        int j = yPointsToBins[idx];
        xyH.bins[i][j].push_back(idx);
// cout << "[" << X[idx] << ", " << Y[idx] << "] -> " << "{" << xyH.xBinEdges[i] << ":" << xyH.xBinEdges[i+1] << ", " << xyH.yBinEdges[j] << ":" << xyH.yBinEdges[j+1] << "}" << endl;
      }
    }
  }

  return xyH;
}

dTERMen::oneDimHist dTERMen::binData(const vector<mstreal>& X, int binSpecType, const vector<mstreal>& binSpec, const vector<mstreal>& M, bool isAngle) {
  // limit to data points within range and transform dihedral angles
  mstreal minVal = binSpec[0];
  mstreal maxVal = binSpec[1];
  vector<mstreal> x = X, m = M;
  vector<int> origIndices = MstUtils::range(0, (int) X.size());
  vector<int> exclude;
  for (int i = 0; i < x.size(); i++) {
    if ((x[i] < minVal) || (x[i] > maxVal)) {
      exclude.push_back(i);
      continue;
    }
    if (isAngle) {
      mstreal an = CartesianGeometry::angleDiff(x[i], 0);
      // represent +/- pi as -pi; then, the bin boundary definition always works
      if (an >= 180) an = -180;
      x[i] = an;
    }
  }
  if (m.empty()) m.resize(x.size(), 1.0);
  if (!exclude.empty()) {
    vector<mstreal> cleanedData(x.size() - exclude.size());
    vector<mstreal> cleanedMult(cleanedData.size());
    vector<int> cleanedIndices(cleanedData.size());
    int k = 0, j = 0;
    for (int i = 0; i < x.size(); i++) {
      if (i == exclude[k]) { k++; continue; }
      cleanedData[j] = x[i];
      cleanedMult[j] = m[i];
      cleanedIndices[j] = origIndices[i];
      j++;
    }
    x = cleanedData;
    origIndices = cleanedIndices;
    m = cleanedMult;
  }

  // actually do binning
  vector<int> sortedInds = MstUtils::sortIndices(x);
  oneDimHist H;
  if (binSpecType == 1) {
    /* Uniform binning */
    if (binSpec.size() != 3) MstUtils::error("expected three values for bin specification type 1", "dTERMen::binData");
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
  for (int i = 0; i < X.size(); i++) { // copy all weights (even for excluded points)
    H.weights[i] = M.empty() ? 1.0 : 1.0/M[i];
  }
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
  }

  return H;
}

dTERMen::zeroDimPotType dTERMen::buildZeroDimPotential(const vector<int>& AA, vector<vector<mstreal> >& backPot) {
  // get the potential
  zeroDimPotType pot;
  int naa = globalAlphabetSize();
  map<int, mstreal> aaCounts;
  for (int i = 0; i < naa; i++) aaCounts[i] = 0.0;
  for (int i = 0; i < AA.size(); i++) aaCounts[AA[i]] += 1;
  pot.aaEnergies = vector<mstreal>(naa, 1.0/0.0);
  for (int i = 0; i < naa; i++) {
    pot.aaEnergies[i] = -kT*log(aaCounts[i]/AA.size());
  }

  // compute the background energy
  for (int i = 0; i < AA.size(); i++) {
    for (int j = 0; j < naa; j++) {
      backPot[i][j] = pot.aaEnergies[j];
    }
  }

  return pot;
}

dTERMen::oneDimPotType dTERMen::buildOneDimPotential(const oneDimHist& H, const vector<int>& AA, mstreal pc, vector<vector<mstreal> >& backPot, bool updateBackPot) {
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
        for (int aai = 0; aai < naa; aai++) Z += exp(-backPot[k][aai]/kT);
        for (int aai = 0; aai < naa; aai++) expCounts[aai] += w*exp(-backPot[k][aai]/kT)/Z;
      } else {
        // if no prior potential given, assume uniform prior
        for (int aai = 0; aai < naa; aai++) expCounts[aai] += w/naa;
      }
    }

    // compute the potential using frequency-proportional pseudocounts
    for (int aa = 0; aa < naa; aa++) {
      pot.aaEnergies[bi][aa] = -kT*log((pot.aaEnergies[bi][aa] + pc*fAA[aa]*naa)/(expCounts[aa] + pc*fAA[aa]*naa));
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


dTERMen::twoDimPotType dTERMen::buildTwoDimPotential(const twoDimHist& H, const vector<int>& AA, mstreal pc, vector<vector<mstreal> >& backPot, bool updateBackPot) {
  twoDimPotType pot;
  pot.xBinEdges = H.xBinEdges;
  pot.yBinEdges = H.yBinEdges;
  int xnb = H.xBinEdges.size() - 1;
  int ynb = H.yBinEdges.size() - 1;
  int naa = globalAlphabetSize();
  pot.aaEnergies.resize(xnb, vector<vector<mstreal> >(ynb, vector<mstreal>(naa, 0.0)));

  // overall amino-acid frequencies
  CartesianPoint fAA(naa);
  for (int bi = 0; bi < xnb; bi++) { // iterating over bin indices ignores any points excluded from binning
    for (int bj = 0; bj < ynb; bj++) {
      const vector<int>& binInds = H.bins[bi][bj];
      for (int i = 0; i < binInds.size(); i++) {
        int k = binInds[i];
        fAA[AA[k]] += H.weights[k];
      }
    }
  }
  fAA /= fAA.sum();

  // add up the weighted counts of each amino acid in each bin
  for (int bi = 0; bi < xnb; bi++) { // iterating over bin indices ignores any points excluded from binning
    for (int bj = 0; bj < ynb; bj++) {
      vector<mstreal> expCounts = vector<mstreal> (naa, 0.0); // expectation of each amino acid in this bin
      const vector<int>& binInds = H.bins[bi][bj];
      for (int i = 0; i < binInds.size(); i++) {
        int k = binInds[i];
        mstreal w = H.weights[k]; // position weight
        int aa = AA[k];
        pot.aaEnergies[bi][bj][aa] += w;

        // compute the expectation of each amino acid in this position
        if (!backPot.empty()) {
          mstreal Z = 0; // partition function for this position (over all amino acids)
          for (int aai = 0; aai < naa; aai++) Z += exp(-backPot[k][aai]/kT);
          for (int aai = 0; aai < naa; aai++) expCounts[aai] += w*exp(-backPot[k][aai]/kT)/Z;
        } else {
          // if no prior potential given, assume uniform prior
          for (int aai = 0; aai < naa; aai++) expCounts[aai] += w/naa;
        }
      }

      // compute the potential using frequency-proportional pseudocounts
      for (int aa = 0; aa < naa; aa++) {
        pot.aaEnergies[bi][bj][aa] = -kT*log((pot.aaEnergies[bi][bj][aa] + pc*fAA[aa]*naa)/(expCounts[aa] + pc*fAA[aa]*naa));
      }
    }
  }

  // update the prior potential
  if (!backPot.empty() && updateBackPot) {
    for (int bi = 0; bi < xnb; bi++) { // iterating over bin indices ignores any points excluded from binning
      for (int bj = 0; bj < ynb; bj++) {
        const vector<int>& binInds = H.bins[bi][bj];
        for (int i = 0; i < binInds.size(); i++) {
          int k = binInds[i];
          for (int aai = 0; aai < naa; aai++) {
            backPot[k][aai] += pot.aaEnergies[bi][bj][aai];
          }
        }
      }
    }
  }

  return pot;
}

dTERMen::oneDimPotType dTERMen::buildOneDimPotential(const oneDimHist& H, const vector<int>& AA, mstreal pc) {
  vector<vector<mstreal> > backPot;
  return buildOneDimPotential(H, AA, pc, backPot, false);
}

void dTERMen::printZeroDimPotential(const zeroDimPotType& P) {
  int naa = P.aaEnergies.size();
  for (int aa = 0; aa < naa; aa++) {
    cout << indexToResName(aa) << "\t" << P.aaEnergies[aa] << endl;
  }
}

void dTERMen::printOneDimPotential(const oneDimPotType& P) {
  int nb = P.aaEnergies.size();
  if (nb == 0) return;
  int naa = globalAlphabetSize();
  for (int aa = 0; aa < naa; aa++) cout << indexToResName(aa) << " ";
  cout << endl;
  for (int bi = 0; bi < nb; bi++) {
    printf("%f %f", P.binEdges[bi], P.binEdges[bi+1]);
    for (int aa = 0; aa < naa; aa++) {
      printf(" %5.3f", P.aaEnergies[bi][aa]);
    }
    printf("\n");
  }
}

void dTERMen::printTwoDimPotential(const twoDimPotType& P) {
  int xnb = P.xBinEdges.size() - 1;
  int ynb = P.yBinEdges.size() - 1;
  int naa = globalAlphabetSize();
  for (int aa = 0; aa < naa; aa++) {
    cout << "---------------> " << indexToResName(aa) << ":" << endl;
    for (int i = 0; i < xnb; i++) {
      for (int j = 0; j < ynb; j++) {
        printf("%5.3f ", P.aaEnergies[i][j][aa]);
      }
      printf("\n");
    }
  }
}

void dTERMen::readBackgroundPotentials(const string& file) {

}

mstreal dTERMen::selfEnergy(Residue* R, const string& aa) {
  return (selfEnergies(R))[dTERMen::aaToIndex(aa.empty() ? R->getName() : aa)];
}

vector<mstreal> dTERMen::selfEnergies(Residue* R, bool verbose) {
  if (R->getStructure() == NULL) MstUtils::error("cannot operate on a disembodied residue!", "dTERMen::selfEnergy");
  if (verbose) cout << "\tdTERMen::selfEnergies -> getting contacts for " << *R << "..." << endl;
  Structure& S = *(R->getStructure());
  auto rmsdCutSelfRes = [](const vector<int>& fragResIdx, const Structure& S) { return RMSDCalculator::rmsdCutoff(fragResIdx, S, 1.0, 20); };
  auto rmsdCutSelfCor = [](const vector<int>& fragResIdx, const Structure& S) { return RMSDCalculator::rmsdCutoff(fragResIdx, S, 1.1, 15); };

  // -- get contacts
  ConFind C(&RL, S);
  contactList cL = C.getContacts(R, cdCut);
  cL.sortByDegree();
  vector<Residue*> contResidues = cL.dstResidues();

  // -- simple environment components
  if (verbose) cout << "\tdTERMen::selfEnergies -> trivial statistical components..." << endl;
  int naa = globalAlphabetSize();
  CartesianPoint selfE(naa, 0.0);
  // vector<mstreal> selfE(naa, 0.0);
  for (int aai = 0; aai < naa; aai++) {
    selfE[aai] = backEner(aai) + bbOmegaEner(R->getOmega(), aai) + bbPhiPsiEner(R->getPhi(), R->getPsi(), aai) + envEner(C.getFreedom(R), aai);
  }

  // -- self residual
  if (verbose) cout << "\tdTERMen::selfEnergies -> self residual..." << endl;
  Structure selfTERM;
  vector<int> fragResIdx;
  int cInd = TERMUtils::selectTERM({R}, selfTERM, pmSelf, &fragResIdx)[0];
  F.setOptions(fasstSearchOptions());
  F.setQuery(selfTERM);
  F.setRMSDCutoff(rmsdCutSelfRes(fragResIdx, S));
  F.setMinNumMatches(selfResidualMinN);
  F.setMaxNumMatches(selfResidualMaxN);
  F.setRedundancyProperty("sim");
  fasstSolutionSet matches = F.search();
  selfE += singleBodyStatEnergy(matches, cInd, selfResidualPC);

  // -- self correction
  if (verbose) cout << "\tdTERMen::selfEnergies -> self correction..." << endl;
  class clique {
    // represents a local interaction "clique" for self-correction calculations.
    // NOTE: the word "clique" is used loosely here, "connected component" would
    // probably be more accurate, because we don't actually look at contacts
    // between the residues that contact the central residue. BUT, the final
    // set of motifs we hope to end up with can still be thought of as cliques,
    // because each motif occurs together frequently enough (so the residues
    // involved in the motif "go together" frequently -- i.e., form a "clique").
    public:
      vector<Residue*> residues;
      int centResIdx;
      fasstSolutionSet matches;
      string toString() {
        stringstream ss;
        ss << "clique with " << residues.size() << " residues, centered on " << centResIdx << ", " << matches.size() << " matches:";
        for (int i = 0; i < residues.size(); i++) ss << " " << *(residues[i]);
        return ss.str();
      }
  };

  // consider each contacting residue with the central one and see if there are
  // enough matches. If so, create a clique to be grown later. If not, still
  // create a clique (with number of matches), which will not be grown.
  vector<clique> finalCliques;
  map<Residue*, clique> cliquesToGrow;
  for (int i = 0; i < contResidues.size(); i++) {
    if (verbose) cout << "\t\tdTERMen::selfEnergies -> seed clique with contact " << *(contResidues[i]) << "..." << endl;
    clique c;
    c.residues = {R, contResidues[i]};
    vector<int> fragResIdx;
    Structure term;
    c.centResIdx = TERMUtils::selectTERM(c.residues, term, pmSelf, &fragResIdx)[0];
    F.setOptions(fasstSearchOptions());
    F.setQuery(term);
    mstreal cut = rmsdCutSelfCor(fragResIdx, S);
    F.setRMSDCutoff(cut);
    F.setMinNumMatches(selfCorrMinN);
    F.setMaxNumMatches(selfCorrMaxN);
    c.matches = F.search();
    if (matches[selfCorrMinN - 1].getRMSD() > cut) { finalCliques.push_back(c); }
    else { cliquesToGrow[contResidues[i]] = c; }
  }

  // for those contacts with sufficient matches, try to combine into larger cliques
  while (!cliquesToGrow.empty()) {
    // sort remaining contacting residues by decreasing number of matches of the
    // clique comprising the residue and the central residue
    vector<Residue*> remConts = MstUtils::keys(cliquesToGrow);
    sort(remConts.begin(), remConts.end(), [&cliquesToGrow](Residue* i, Residue* j) { return cliquesToGrow[i].matches.size() > cliquesToGrow[j].matches.size(); });

    // pick the one with highest number of contacts, and try to grow the clique
    clique parentClique = cliquesToGrow[remConts[0]];
    remConts.erase(remConts.begin());
    clique grownClique;
    if (verbose) cout << "\t\tdTERMen::selfEnergies -> will try to grow [" << parentClique.toString() << "]..." << endl;
    while (!remConts.empty()) {
      // try to add every remaining contact
      for (int j = 0; j < remConts.size(); j++) {
        if (verbose) cout << "\t\t\tdTERMen::selfEnergies -> trying to add " << *(remConts[j]) << "..." << endl;
        clique newClique;
        newClique.residues = parentClique.residues;
        newClique.residues.push_back(remConts[j]);
        Structure term; vector<int> fragResIdx;
        newClique.centResIdx = TERMUtils::selectTERM(newClique.residues, term, pmSelf, &fragResIdx)[0];
        F.setOptions(fasstSearchOptions());
        F.setQuery(term);
        F.setRMSDCutoff(rmsdCutSelfCor(fragResIdx, S));
        F.setMaxNumMatches(selfCorrMaxN);
        newClique.matches = F.search();
        if ((j == 0) || (newClique.matches.size() > grownClique.matches.size())) {
          if (verbose) cout << "\t\t\t\tdTERMen::selfEnergies -> new best" << endl;
          grownClique = newClique;
        }
      }
      // EITHER a sufficient number of matches OR not too many fewer than for the
      // parent clique is the logic we had implemented in original dTERMen
      if ((grownClique.matches.size() >= selfCorrMinN) || (grownClique.matches.size() > 0.8*parentClique.matches.size())) {
        remConts = MstUtils::setdiff(remConts, {grownClique.residues.back()});
        parentClique = grownClique;
        if (verbose) cout << "\t\t\tdTERMen::selfEnergies -> chose to add " << *(grownClique.residues.back()) << "..." << endl;
      } else {
        if (verbose) cout << "\t\t\tdTERMen::selfEnergies -> nothing more worked, sticking with parent clique..." << endl;
        grownClique = parentClique;
        break;
      }
    }

    // add the grown clique to the list of final cliques, and remove used up
    // residues from the list of cliques to grow
    if (verbose) cout << "\t\tdTERMen::selfEnergies -> adding final clique [" << grownClique.toString() << "]..." << endl;
    finalCliques.push_back(grownClique);
    vector<Residue*> deleted = MstUtils::setdiff(MstUtils::keys(cliquesToGrow), remConts);
    if (verbose) cout << "\t\tdTERMen::selfEnergies -> neighboring residues taken care of this cycle:";
    for (int i = 0; i < deleted.size(); i++) {
      if (verbose) cout << " " << *(deleted[i]);
      cliquesToGrow.erase(deleted[i]);
    }
    if (verbose) cout << endl;
  }

  // finally, actually compute the self correction for each final clique
  if (verbose) cout << "\tdTERMen::selfEnergies -> final cliques:" << endl;
  for (int i = 0; i < finalCliques.size(); i++) {
    if (verbose) cout << "\t\t" << finalCliques[i].toString() << endl;
    selfE += singleBodyStatEnergy(finalCliques[i].matches, finalCliques[i].centResIdx, selfCorrPC);
  }

  return selfE;
}

CartesianPoint dTERMen::singleBodyStatEnergy(fasstSolutionSet& matches, int cInd, int pc) {
  int naa = globalAlphabetSize();
  CartesianPoint selfE(naa, 0.0);
  vector<mstreal> No(naa, pc), Ne(naa, pc);
  vector<mstreal> p(naa, 0);
  mstreal phi, psi, omg, env, ener;
  for (int i = 0; i < matches.size(); i++) {
    // see what amino acid is actually observed at the position
    const fasstSolution& m = matches[i];
    res_t aa = F.getMatchSequence(m)[cInd];
    int aaIdx = dTERMen::aaToIndex(aa);
    if (aaIdx < 0) continue; // this match has some amino acid that is not known in the current alphabet
    No[aaIdx] += 1.0;

    // now try out all amino acids at this position and compute expectations
    phi = F.getResidueProperties(m, "phi")[cInd];
    psi = F.getResidueProperties(m, "psi")[cInd];
    omg = F.getResidueProperties(m, "omega")[cInd];
    env = F.getResidueProperties(m, "env")[cInd];
    for (int aai = 0; aai < naa; aai++) {
      p[aai] = backEner(aai) + bbOmegaEner(omg, aai) + bbPhiPsiEner(phi, psi, aai) + envEner(env, aai);
    }
    p = dTERMen::enerToProb(p);
    for (int aai = 0; aai < naa; aai++) Ne[aai] += p[aai];
  }
  for (int aai = 0; aai < naa; aai++) {
    selfE[aai] = -kT*log(No[aai]/Ne[aai]);
  }

  return selfE;
}

vector<mstreal> dTERMen::enerToProb(const vector<mstreal>& ener) {
  mstreal minEner = MstUtils::min(ener);
  mstreal Z = 0;
  vector<mstreal> p(ener.size(), 0);
  for (int i = 0; i < ener.size(); i++) Z += exp(-(ener[i] - minEner)/kT);
  for (int i = 0; i < ener.size(); i++) p[i] = exp(-(ener[i] - minEner)/kT)/Z;
  return p;
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
