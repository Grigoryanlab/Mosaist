#include "mstfuser.h"

vector<string> fusionTopology::bba = {"N", "CA", "C", "O"};

/* --------- fusionParams ------------ */
void fusionParams::setStartingStructure(const Structure& _S) {
  AtomPointerVector atoms;
  vector<Residue*> R = _S.getResidues();
  for (int i = 0; i < R.size(); i++) {
    for (int j = 0; j < fusionTopology::bba.size(); j++) {
      atoms.push_back(R[i]->findAtom(fusionTopology::bba[j]));
    }
  }
  startStruct = Structure(atoms);
}

/* --------- fusionOutput ------------ */
void fusionOutput::addSnapshot(const Structure& snap, mstreal ener) {
  if (snap.atomSize() != fused.atomSize()) MstUtils::error("fused structure and added snapshot have different numbers of atoms", "fusionOutput::addSnapshot");
  AtomPointerVector atoms = snap.getAtoms();
  int na = snap.atomSize();
  int dim = 3;
  for (int k = 0; k < dim; k++) trajSnaps[k].push_back(vector<mstreal>(na, 0.0));
  for (int i = 0; i < na; i++) {
    for (int k = 0; k < dim; k++) trajSnaps[k].back()[i] = (*atoms[i])[k];
  }
  trajScores.push_back(ener);
}

Structure fusionOutput::getSnapshot(int j) const {
  Structure snap = fused;
  AtomPointerVector atoms = snap.getAtoms();
  int dim = 3;
  for (int i = 0; i < atoms.size(); i++) {
    for (int k = 0; k < dim; k++) (*atoms[i])[k] = trajSnaps[k][j][i];
  }
  return snap;
}

/* --------- fusionTopology ------------ */
fusionTopology::fusionTopology(int L) {
  overlappingResidues.resize(L);
  fixed.resize(L, false);
  updated = true;
  verbose = false;
}

fusionTopology::fusionTopology(const vector<vector<Residue*> >& resTopo) {
  overlappingResidues = resTopo;
  fixed.resize(resTopo.size(), false);

  /* find fragments moving together based on the Structure they belong to. */
  map<Structure*, vector<pair<Residue*, int> > > frags;
  for (int i = 0; i < resTopo.size(); i++) {
    for (int j = 0; j < resTopo[i].size(); j++) {
      Structure* S = resTopo[i][j]->getStructure();
      frags[S].push_back(pair<Residue*, int>(resTopo[i][j], i));
    }
  }
  for (auto it = frags.begin(); it != frags.end(); ++it) {
    vector<pair<Residue*, int> >& residues = it->second;
    AtomPointerVector fragAtoms; vector<int> fragResIdx;
    for (int i = 0; i < residues.size(); i++) {
      int ri = residues[i].second;
      Residue& fragRes = *(residues[i].first);
      fragResIdx.push_back(ri);
      for (int j = 0; j < bba.size(); j++) {
        fragAtoms.push_back(fragRes.findAtom(bba[j]));
      }
    }
    fragments.push_back(pair<AtomPointerVector, vector<int> > (fragAtoms, fragResIdx));
    fragWeights.push_back(1.0);
  }
  updated = true;
  verbose = false;
}

fusionTopology::fusionTopology(const fusionTopology& topo) {
  *this = topo;
}

fusionTopology& fusionTopology::operator=(const fusionTopology& topo) {
  this->fragments = topo.fragments;
  this->fragWeights = topo.fragWeights;
  this->overlappingResidues = topo.overlappingResidues;
  this->fixed = topo.fixed;
  this->alignedFrags = topo.alignedFrags;
  this->numMobAtoms = topo.numMobAtoms;
  this->updated = topo.updated;
  this->overlap = topo.overlap;
  this->fragsByOverlap = topo.fragsByOverlap;
  this->chainLengths = topo.chainLengths;
  this->fixedInChain = topo.fixedInChain;
  this->verbose = topo.verbose;
  return *this;
}

void fusionTopology::addFragment(vector<Residue*>& R, const vector<int>& fragResIdx, mstreal weight) {
  if (fragResIdx.size() != 0) MstUtils::assertCond(R.size() == fragResIdx.size(), "fragment residue index vector not the same length as the number of residues in the fragment", "fusionTopology::addFragment(vector<Residue*>&, mstreal, const vector<int>&)");
  AtomPointerVector fragAtoms; vector<int> fragRes = fragResIdx;
  for (int i = 0; i < R.size(); i++) {
    for (int j = 0; j < bba.size(); j++) {
      fragAtoms.push_back(R[i]->findAtom(bba[j]));
    }
    if (fragResIdx.size() == 0) fragRes.push_back(R[i]->getNum());
    overlappingResidues[fragRes[i]].push_back(R[i]);
  }
  fragments.push_back(pair<AtomPointerVector, vector<int> > (fragAtoms, fragRes));
  fragWeights.push_back(weight);
  set<int> fragResSet = MstUtils::contents(fragRes);
  if (overlap.find(fragResSet) == overlap.end()) {
    overlap[fragResSet] = fragsByOverlap.size();
    fragsByOverlap.resize(fragsByOverlap.size() + 1);
  }
  fragsByOverlap[overlap[fragResSet]].push_back(fragments.size() - 1);
  updated = true;
}

void fusionTopology::addFragment(Structure& S, const vector<int>& fragResIdx, mstreal weight) {
  vector<Residue*> residues = S.getResidues();
  addFragment(residues, fragResIdx, weight);
}

void fusionTopology::setAlignedFrags(Structure& fused) {
  alignedFrags.resize(0);
  alignedFrags.resize(fragments.size());
  for (int i = 0; i < fragments.size(); i++) {
    AtomPointerVector& refFrag = fragments[i].first;
    AtomPointerVector fusedFrag;
    vector<int>& fragResIdx = fragments[i].second;
    for (int j = 0; j < fragResIdx.size(); j++) {
      Residue& fusedRes = fused.getResidue(fragResIdx[j]);
      for (int k = 0; k < bba.size(); k++) {
        fusedFrag.push_back(fusedRes.findAtom(bba[k]));
      }
    }
    alignedFrags[i].first = fusedFrag;
    alignedFrags[i].second = refFrag;
  }
}

bool fusionTopology::isConsistent(const Structure& S) {
  if (S.residueSize() != this->length()) return false;
  vector<int> chainLengths = getChainLengths();
  if (S.chainSize() != chainLengths.size()) return false;
  for (int i = 0; i < chainLengths.size(); i++) {
    if (S[i].residueSize() != chainLengths[i]) return false;
  }
  return true;
}

Structure fusionTopology::remapChains(const Structure& S) {
  Structure Sc;
  if (S.residueSize() != this->length()) MstUtils::error("cannot remap chains because number of residues in Structure does not match that in topology");
  vector<int> chainLengths = getChainLengths();
  int k = 0;
  for (int ci = 0; ci < chainLengths.size(); ci++) {
    Chain* C = Sc.appendChain("A", true);
    for (int ri = 0; ri < chainLengths[ci]; ri++) {
      C->appendResidue(new Residue(S.getResidue(k)));
      k++;
    }
  }
  return Sc;
}

void fusionTopology::addFixedPositions(vector<int> fixedInds) {
  for (int i = 0; i < fixedInds.size(); i++) addFixedPosition(fixedInds[i]);
}

vector<int> fusionTopology::getFixedPositions() {
  vector<int> fixedPositions;
  for (int i = 0; i < fixed.size(); i++) {
    if (fixed[i]) fixedPositions.push_back(i);
  }
  return fixedPositions;
}

int fusionTopology::numFixedPositions() {
  int num = 0;
  for (int i = 0; i < fixed.size(); i++) { if (fixed[i]) num++; }
  return num;
}

int fusionTopology::numFixedInChain(int ci) {
  updateConnectivity();
  int num = 0;
  for (int i = 0; i < fixedInChain[ci].size(); i++) { if (fixedInChain[ci][i]) num++; }
  return num;
}

void fusionTopology::updateConnectivity() {
  if (!updated) return;
  int L = overlappingResidues.size();
  int nbba = fusionTopology::bba.size();
  vector<bool> connectedToNext(L, false); // whether each residue is bonded to the next one in the topology

  for (int i = 0; i < fragments.size(); i++) {
    AtomPointerVector& frag = fragments[i].first;
    vector<int> posInds = fragments[i].second;
    // walk over residues of the fragment in the order these residues appear in the topology
    vector<int> sortedInds = MstUtils::sortIndices(posInds);
    for (int k = 0; k < sortedInds.size() - 1; k++) {
      // if this and the next one are consecutive in the topology AND the residues
      // in the fragment are bonded, then mark a connection. NOTE: it is sufficient
      // for the two residues to be bonded in just one example fragment for us to
      // conclude that there is a chain connection (this is by construction).
      int ri = sortedInds[k];
      int rj = sortedInds[k+1];
      if ((posInds[ri] == posInds[rj] - 1) && (Residue::areBonded(frag[nbba*ri]->getResidue(), frag[nbba*rj]->getResidue()))) {
        connectedToNext[posInds[ri]] = true;
      }
    }
  }

  // now walk over the topology and cut chains when there is no connectivity
  int cL = 0; vector<bool> fixedInCurrChain;
  for (int i = 0; i < L; i++) {
    cL++;
    if (fixed[i]) fixedInCurrChain.push_back(true);
    else fixedInCurrChain.push_back(false);
    if (!connectedToNext[i] && (!fixed[i] || (i == L-1) || !fixed[i+1])) {
      if (verbose) cout << "fusionTopology::updateConnectivity -> found a break at topology position " << i << endl;
      chainLengths.push_back(cL);
      fixedInChain.push_back(fixedInCurrChain);
      fixedInCurrChain.clear();
      cL = 0;
    }
  }
  if (verbose) cout << "fusionTopology::updateConnectivity -> in the end, found " << chainLengths.size() << " chains of lengths: " << MstUtils::vecToString(chainLengths) << endl;

  // finally figure out the number of mobile atoms in each chain
  numMobAtoms.resize(chainLengths.size());
  for (int ci = 0; ci < numMobAtoms.size(); ci++) {
    numMobAtoms[ci] = chainLengths[ci];
    for (int i = 0; i < chainLengths[ci]; i++) {
      if (fixedInChain[ci][i]) numMobAtoms[ci]--;
    }
    numMobAtoms[ci] *= nbba;
  }

  updated = false;
}

int fusionTopology::numMobileAtoms() {
  updateConnectivity();
  int n = 0;
  for (int i = 0; i < numMobAtoms.size(); i++) n += numMobAtoms[i];
  return n;
}

/* --------- fusionEvaluator ----------- */
fusionEvaluator::fusionEvaluator(const fusionTopology& _topo, const fusionParams& _params) {
  params = _params;
  topo = _topo;
  init();
}

void fusionEvaluator::init() {
  // create room for fused structure (initialize with average coordinates)
  if (!params.isStartingStructureGiven()) {
    vector<int> chainLengths = topo.getChainLengths();
    fused.reset();
    int i = 0;
    for (int ci = 0; ci < chainLengths.size(); ci++) {
      Chain* C = fused.appendChain("A", true);
      for (int ri = 0; ri < chainLengths[ci]; ri++, i++) {
        if (topo.numOverlappingResidues(i) == 0) MstUtils::error("position index " + MstUtils::toString(i) + " has no overlapping residues in the specified topology, and no starting structure is specified; cannot begin...", "fusionEvaluator::fusionEvaluator");
        Residue* res = new Residue(topo.getOverlappingResidue(i, 0)->getName(), 1);
        C->appendResidue(res);
        for (int j = 0; j < fusionTopology::bba.size(); j++) {
          /* I previously disallowed multiple residues to be aligned in fixed
           * positions in the topology. BUT, it turns out this can be quite useful.
           * For example, if we have residues 3 and 4 fixed, we may still want to
           * have a fragment covering [2, 3, 4] and another one covering [3, 4, 5],
           * because those fragments would consrain geometries between residues 2
           * and 3 and residues 4 and 5. Sure, 3 and 4 could not move and so could
           * not lower the RMSD score themselves. BUT, 2 and 5 can move and would
           * affect the RMSDs arising from the above fragments. So, since allowing
           * multiple residues aligned onto fixed topology positions does make
           * sense, we have to decide how to initialize the coordinates of these
           * fixed positions. If we take the usual average, then it is difficult
           * for the user to specify what exactly these positions should be fixed
           * to (specifying this via the average would be awkward for the user). So
           * here I decided to adopt a convention, where for fixed positions the
           * coordinates would be initialized from the FIRST available residue.
           * The first one always exists, even if only one residue is aligned. And
           * it is also relatively simple to use the first residue at fixed positions
           * as a way of communicating the fixed coordinates. */
          CartesianPoint m = topo.isFixed(i) ? atomInstances(i, fusionTopology::bba[j])[0] : atomInstances(i, fusionTopology::bba[j]).getGeometricCenter();
          res->appendAtom(new Atom(1, fusionTopology::bba[j], m.getX(), m.getY(), m.getZ(), 0, 1.0, false));
        }
      }
    }
    fused.renumber();
  } else {
    fused = params.getStartingStructure();
    // sometimes, the topology may end up lumping several consecutive fixed chain
    // together, so try to redistribute residues to match the topology chains
    if (!topo.isConsistent(fused)) fused = topo.remapChains(fused);
    fused.renumber();
  }
  guess = fused; // save initial "averaged" guess for later alignment (easier to visualize output)
  topo.setAlignedFrags(fused);
  chooseBuildOrigin(); // pick build origins for each chain

  // optionally, initialize the structure using average IC's of overlapping segments
  if (params.getCoorInitType() == fusionParams::coorInitType::meanIC) {
    bool optCart = params.getOptimCartesian();
    params.setOptimCartesian(false); // temporarily switch to IC coordinate representation
    eval(vector<mstreal>(0));        // fill initial point with IC-built structure
    eval(initPoint);                 // fill the Structure with corresponding coordinates
    initPoint.resize(0);
    params.setOptimCartesian(optCart);
  }

  // orient so that the first atom is at the origin, the second is along the X-
  // axis and the third is in the XY plane
  Frame L(CartesianPoint(0, 0, 0), CartesianPoint(1, 0, 0), CartesianPoint(0, 1, 0), CartesianPoint(0, 0, 1));
  CartesianPoint A(fused[0][0][0]), B(fused[0][0][1]), C(fused[0][0][2]);
  CartesianPoint X = (B - A).getUnit();
  CartesianPoint Z = (X.cross(C - B)).getUnit();
  Frame F(A, X, Z.cross(X), Z);
  Transform T = TransformFactory::switchFrames(L, F);
  T.apply(fused);

  // initialize the random number number generator
  long int x = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
  srand(x);
}

fusionEvaluator::fusionEvaluator(const vector<vector<Residue*> >& resTopo, vector<int> _fixedResidues, const fusionParams& _params) {
  params = _params;
  topo = fusionTopology(resTopo);
  topo.addFixedPositions(_fixedResidues);
  init();
}

void fusionEvaluator::chooseBuildOrigin(bool randomize) {
  vector<int> fixedResidues = topo.getFixedPositions();
  if (randomize) {
    buildOriginRes = (fixedResidues.size() > 0) ? fixedResidues[MstUtils::randInt(0, fixedResidues.size() - 1)] : MstUtils::randInt(0, topo.length() - 1);
  } else {
    buildOriginRes = (fixedResidues.size() > 0) ? MstUtils::min(fixedResidues) : 0;
  }
}

Structure fusionEvaluator::getAlignedStructure() {
  Structure aligned = fused;
  RMSDCalculator rc;
  rc.align(aligned.getAtoms(), guess.getAtoms(), aligned);
  return aligned;
}

mstreal fusionEvaluator::eval(const vector<mstreal>& point) {
  bool init = point.empty();
  if (init) {
    initPoint.resize(0);
    gradOfXYZ.clear();
    gradient.resize(numDF());
    bounds.clear();
    masses.clear();
  }
  mstreal bR = 0.01; mstreal aR = 1.0; mstreal dR = 1.0; mstreal xyzR = 0.01; // randomness scale factors
  if (!init && (point.size() != numDF())) MstUtils::error("need to place " + MstUtils::toString(topo.numMobileAtoms()) + " atoms, with " + MstUtils::toString(topo.numFixedPositions()) + " fixed residues, and received " + MstUtils::toString(point.size()) + " parameters: that appears wrong!", "fusionEvaluator::eval");
  int k = 0;
  Atom *pN = NULL, *pCA = NULL, *pC = NULL, *pO = NULL;
  mstreal noise = params.getNoise();

  if (params.getOptimCartesian()) {
    for (int i = 0; i < fused.residueSize(); i++) {
      if (topo.isFixed(i)) continue;
      Residue& res = fused.getResidue(i);
      for (int j = 0; j < res.atomSize(); j++) {
        bool skipDFs = ((i == 0) && !isAnchored());
        for (int dim = 0; dim < 3; dim++) {
          /* If the fused structure is not anchored in space, skip all coordinates
           * of the first atom, the Y and the Z coordinates of the second atom,
           * and the Z coordinate of the third atom. In this case, the constructor
           * would have placed the first atom at the origin, the second atom on
           * the X-axis, and the third atom in the X-Y plane. */
          if (skipDFs && (j < 3) && (dim >= j)) continue;
          if (init) {
            initPoint.push_back(res[j][dim] + xyzR * MstUtils::randUnit() * noise);
            gradOfXYZ.addPartial(res[j], dim, k, 1.0);
            masses.push_back(res[j].getMass());
          } else {
            res[j][dim] = point[k];
          }
          k++;
        }
      }
    }
  } else {
    // compute all reduced masses
    mstreal m_N_CA, m_CA_C, m_N_C, m_C_O, m_N_CA_C, m_CA_C_O, m_N_C_O, m_N_CA_C_N, m_CA_C_N_CA, m_C_N_CA_C, m_N_CA_C_O;
    if (init) { // these "seem" right, but I have not checked (good enough for now)
      mstreal mN = Atom::getMass("N");
      mstreal mCA = Atom::getMass("CA");
      mstreal mC = Atom::getMass("C");
      mstreal mO = Atom::getMass("O");
      m_N_CA = mN*mCA/(mN + mCA);
      m_CA_C = mCA*mC/(mCA + mC);
      m_N_C = mN*mC/(mN + mC);
      m_C_O = mC*mO/(mC + mO);
      m_N_CA_C = mN*mCA*mC/(mN*mCA + mCA*mC + mN*mC);
      m_CA_C_O = mCA*mC*mO/(mCA*mC + mCA*mO + mC*mO);
      m_N_C_O = mN*mC*mO/(mN*mC + mN*mO + mC*mO);
      m_N_CA_C_N = mN*mCA*mC*mN/(mN*mCA*mC + mN*mCA*mN + mN*mC*mN + mCA*mC*mN);
      m_C_N_CA_C = mC*mN*mCA*mC/(mC*mN*mCA + mC*mN*mC + mC*mCA*mC + mN*mCA*mC);
      m_CA_C_N_CA = mCA*mC*mN*mCA/(mCA*mC*mN + mCA*mC*mCA + mCA*mN*mCA + mC*mN*mCA);
      m_N_CA_C_O = mN*mCA*mC*mO/(mN*mCA*mC + mN*mCA*mO + mN*mC*mO + mCA*mC*mO);
    }
    // -- build the fused backbone, atom by atom
    // build forward
    int startIdx = isAnchored() ? buildOriginRes : 0;
    for (int i = startIdx; i < fused.residueSize(); i++) {
      Residue& res = fused.getResidue(i);
      Atom* N = &(res[0]);
      Atom* CA = &(res[1]);
      Atom* C = &(res[2]);
      Atom* O = &(res[3]);
      if (!topo.isFixed(i)) {
        if ((i == 0) && !isAnchored()) {
          if (init) {
            initPoint.push_back(bondInitValue(i, i, "N", "CA") + bR * MstUtils::randUnit() * noise);
            masses.push_back(m_N_CA);
            initPoint.push_back(bondInitValue(i, i, "CA", "C") + bR * MstUtils::randUnit() * noise);
            masses.push_back(m_CA_C);
            initPoint.push_back(angleInitValue(i, i, i, "N", "CA", "C") + aR * MstUtils::randUnit() * noise);
            masses.push_back(m_N_CA_C);
          } else {
            mstreal r2d = M_PI/180;
            mstreal d0 = point[k];
            mstreal d1 = point[k+1];
            mstreal a0 = point[k+2]*r2d;
            N->setCoor(0.0, 0.0, 0.0);
            CA->setCoor(d0, 0.0, 0.0);
            gradOfXYZ.addPartial(CA, 0, k, 1.0);
            C->setCoor(d0 - d1*cos(a0), d1*sin(a0), 0.0);
            gradOfXYZ.addPartial(C, 0, k, 1.0);
            gradOfXYZ.addPartial(C, 0, k+1, -cos(a0));
            gradOfXYZ.addPartial(C, 0, k+2, sin(a0)*r2d);
            gradOfXYZ.addPartial(C, 1, k+1, sin(a0));
            gradOfXYZ.addPartial(C, 1, k+2, cos(a0)*r2d);
            k += 3;
          }
        } else {
          // place N relative to pN, pCA, pC
          if (init) {
            initPoint.push_back(bondInitValue(i-1, i, "C", "N") + bR * MstUtils::randUnit() * noise);
            masses.push_back(m_N_C);
            initPoint.push_back(angleInitValue(i-1, i-1, i, "CA", "C", "N") + aR * MstUtils::randUnit() * noise);
            masses.push_back(m_N_CA_C);
            initPoint.push_back(dihedralInitValue(i-1, i-1, i-1, i, "N", "CA", "C", "N") + dR * MstUtils::randUnit() * noise);
            masses.push_back(m_N_CA_C_N);
          } else {
            N->build(pC, pCA, pN, point[k], point[k+1], point[k+2]);
            // TODO: when we want to enable gradient calculation
            // // partial derivatives of XYZ coordinates of the placed atom with
            // // respect to bond, angle, and dihedral
            // vector<mstreal> icPartials(9, 0.0);
            // // partial derivatives of XYZ coordinates of the placed atom with
            // // respect to XYZ coordinates of each of the anchoring atoms
            // vector<vector<mstreal> > xyzPartials(3, vector<mstreal>(9, 0.0));
            // N->build(pC, pCA, pN, point[k], point[k+1], point[k+2], partials, xyzPartials);
            // gradOfXYZ.addPartials(N, {k, k+1, k+2}, icPartials);
            // gradOfXYZ.addRecursivePartials(N, pC, xyzPartials[0]);
            // gradOfXYZ.addRecursivePartials(N, pCA, xyzPartials[1]);
            // gradOfXYZ.addRecursivePartials(N, pN, xyzPartials[2]);
            k += 3;
          }

          // place CA relative to pCA, pC, N
          if (init) {
            initPoint.push_back(bondInitValue(i, i, "N", "CA") + bR * MstUtils::randUnit() * noise);
            masses.push_back(m_N_CA);
            initPoint.push_back(angleInitValue(i-1, i, i, "C", "N", "CA") + aR * MstUtils::randUnit() * noise);
            masses.push_back(m_N_CA_C);
            initPoint.push_back(dihedralInitValue(i-1, i-1, i, i, "CA", "C", "N", "CA") + dR * MstUtils::randUnit() * noise);
            masses.push_back(m_CA_C_N_CA);
          } else {
            CA->build(N, pC, pCA, point[k], point[k+1], point[k+2]); k += 3;
          }

          // place C relative to pC, N, CA
          if (init) {
            initPoint.push_back(bondInitValue(i, i, "CA", "C") + bR * MstUtils::randUnit() * noise);
            masses.push_back(m_CA_C);
            initPoint.push_back(angleInitValue(i, i, i, "N", "CA", "C") + aR * MstUtils::randUnit() * noise);
            masses.push_back(m_N_CA_C);
            initPoint.push_back(dihedralInitValue(i-1, i, i, i, "C", "N", "CA", "C") + dR * MstUtils::randUnit() * noise);
            masses.push_back(m_C_N_CA_C);
          } else {
            C->build(CA, N, pC, point[k], point[k+1], point[k+2]); k += 3;
          }

          // if this is the last residue, place the O relative to N-CA-C
          if (i == fused.residueSize() - 1) {
            if (init) {
              initPoint.push_back(bondInitValue(i, i, "C", "O") + bR * MstUtils::randUnit() * noise);
              masses.push_back(m_C_O);
              initPoint.push_back(angleInitValue(i, i, i, "CA", "C", "O") + aR * MstUtils::randUnit() * noise);
              masses.push_back(m_CA_C_O);
              initPoint.push_back(dihedralInitValue(i, i, i, i, "N", "CA", "C", "O") + dR * MstUtils::randUnit() * noise);
              masses.push_back(m_N_CA_C_O);
            } else {
              O->build(C, CA, N, point[k], point[k+1], point[k+2]); k += 3;
            }
          }
        }
      }
      // place previous O relative to pCA, N, pC (an improper)
      if ((i > startIdx) && !topo.isFixed(i-1)) {
        if (init) {
          initPoint.push_back(bondInitValue(i-1, i-1, "C", "O") + bR * MstUtils::randUnit() * noise);
          masses.push_back(m_C_O);
          initPoint.push_back(angleInitValue(i, i-1, i-1, "N", "C", "O") + aR * MstUtils::randUnit() * noise);
          masses.push_back(m_CA_C_O);
          initPoint.push_back(dihedralInitValue(i-1, i, i-1, i-1, "CA", "N", "C", "O") + dR * MstUtils::randUnit() * noise);
          masses.push_back(m_N_CA_C_O);
        } else {
          pO->build(pC, N, pCA, point[k], point[k+1], point[k+2]); k += 3;
        }
      }
      pN = N; pCA = CA; pC = C; pO = O;
    }

    // build backwards (only would happens when there is an anchor and it is not the 0-th residue)
    startIdx = isAnchored() ? buildOriginRes : -1;
    for (int i = startIdx; i >= 0; i--) {
      Residue& res = fused.getResidue(i);
      Atom* N = &(res[0]);
      Atom* CA = &(res[1]);
      Atom* C = &(res[2]);
      Atom* O = &(res[3]);
      if (!topo.isFixed(i)) {
        // place C relative to pC, pCA, pN
        if (init) {
          initPoint.push_back(bondInitValue(i+1, i, "N", "C") + bR * MstUtils::randUnit() * noise);
          masses.push_back(m_N_C);
          initPoint.push_back(angleInitValue(i+1, i+1, i, "CA", "N", "C") + aR * MstUtils::randUnit() * noise);
          masses.push_back(m_N_CA_C);
          initPoint.push_back(dihedralInitValue(i+1, i+1, i+1, i, "C", "CA", "N", "C") + dR * MstUtils::randUnit() * noise);
          masses.push_back(m_C_N_CA_C);
        } else {
          C->build(pN, pCA, pC, point[k], point[k+1], point[k+2]); k += 3;
        }

        // place CA relative to pCA, pN, C
        if (init) {
          initPoint.push_back(bondInitValue(i, i, "C", "CA") + bR * MstUtils::randUnit() * noise);
          masses.push_back(m_CA_C);
          initPoint.push_back(angleInitValue(i+1, i, i, "N", "C", "CA") + aR * MstUtils::randUnit() * noise);
          masses.push_back(m_N_CA_C);
          initPoint.push_back(dihedralInitValue(i+1, i+1, i, i, "CA", "N", "C", "CA") + dR * MstUtils::randUnit() * noise);
          masses.push_back(m_CA_C_N_CA);
        } else {
          CA->build(C, pN, pCA, point[k], point[k+1], point[k+2]); k += 3;
        }

        // place N relative to pN, C, CA
        if (init) {
          initPoint.push_back(bondInitValue(i, i, "CA", "N") + bR * MstUtils::randUnit() * noise);
          masses.push_back(m_N_CA);
          initPoint.push_back(angleInitValue(i, i, i, "C", "CA", "N") + aR * MstUtils::randUnit() * noise);
          masses.push_back(m_N_CA_C);
          initPoint.push_back(dihedralInitValue(i+1, i, i, i, "N", "C", "CA", "N") + dR * MstUtils::randUnit() * noise);
          masses.push_back(m_N_CA_C_N);
        } else {
          N->build(CA, C, pN, point[k], point[k+1], point[k+2]); k += 3;
        }

        // place O relative to pCA, pN, C (an improper)
        if (init) {
          initPoint.push_back(bondInitValue(i, i, "C", "O") + bR * MstUtils::randUnit() * noise);
          masses.push_back(m_C_O);
          initPoint.push_back(angleInitValue(i+1, i, i, "N", "C", "O") + aR * MstUtils::randUnit() * noise);
          masses.push_back(m_N_C_O);
          initPoint.push_back(dihedralInitValue(i, i+1, i, i, "CA", "N", "C", "O") + dR * MstUtils::randUnit() * noise);
          masses.push_back(m_N_CA_C_O);
        } else {
          O->build(C, pN, CA, point[k], point[k+1], point[k+2]); k += 3;
        }
      }
      pN = N; pCA = CA; pC = C; pO = O;
    }
  }

  // built up a list of restraining ICs
  if (init) {
    k = 0;
    for (int ci = 0; ci < fused.chainSize(); ci++) {
      Chain& chain = fused[ci];
      for (int i = 0; i < chain.residueSize(); i++, k++) {
        Residue& res = chain[i];
        Atom* N = &(res[0]);
        Atom* CA = &(res[1]);
        Atom* C = &(res[2]);
        Atom* O = &(res[3]);
        // if residue not fixed, evaluate its internal coordinates
        if (!topo.isFixed(ci, i)) {
          bondInstances(k, N, CA);
          bondInstances(k, CA, C);
          bondInstances(k, C, O);
          angleInstances(k, N, CA, C);
          // if last residue, constrain O relative to this residue (as opposed to the next one)
          if (i == chain.residueSize() - 1) {
            angleInstances(k, CA, C, O);
            dihedralInstances(k, N, CA, C, O);
          }
        }

        // if either the previous residue or the current one is not fixed, evaluate
        // the internal coordinates connecting the two
        if ((i > 0) && (!topo.isFixed(ci, i-1) || !topo.isFixed(ci, i))) {
          bondInstances(k, pC, N);
          angleInstances(k, pCA, pC, N);
          dihedralInstances(k, pN, pCA, pC, N);
          angleInstances(k, pC, N, CA);
          dihedralInstances(k, pCA, pC, N, CA);
          dihedralInstances(k, pC, N, CA, C);
          // use this residue to constrain the previous O
          if (i < chain.residueSize() - 1) {
            angleInstances(k, pO, pC, N);
            dihedralInstances(k-1, pCA, N, pC, pO);
          }
        }
        pN = N; pCA = CA; pC = C; pO = O;
      }
    }

    if (params.isCompOn() || params.isRepOn()) {
      for (int i = 0; i < fused.residueSize(); i++) {
        Residue& resi = fused.getResidue(i);
        for (int j = i+2; j < fused.residueSize(); j++) { // no interactions between atoms of adjacent residues
          if (topo.isFixed(i) && topo.isFixed(j)) continue;
          Residue& resj = fused.getResidue(j);
          for (int ai = 0; ai < resi.atomSize(); ai++) {
            for (int aj = 0; aj < resj.atomSize(); aj++) {
              // to save time, make a distance IC for this pair of atoms only if
              // they are somewhat close to prohibited distance ranges
              if (params.isRepOn()) {
                mstreal rmin = 0.9*(atomRadius(resi[ai]) + atomRadius(resj[aj]));
                if (resi[ai].distance(resj[aj]) < 2*rmin) {
                  bounds.push_back(icBound(icType::icDistRep, rmin, 0, {&(resi[ai]), &(resj[aj])}));
                }
              }
              if (params.isCompOn()) {
                mstreal rmax = 2*params.getCompRad();
                if (resi[ai].distance(resj[aj]) > rmax/2) {
                  bounds.push_back(icBound(icType::icDistComp, 0, rmax, {&(resi[ai]), &(resj[aj])}));
                }
              }
            }
          }
        }
      }
    }
  }
  if (init) return 0.0;

  // compute penalty score for out-of-range ICs and best-fit RMSDs
  resetScore();
  for (int k = 0; k < bounds.size(); k++) scoreIC(bounds[k]);
  scoreRMSD();

  return score;
}

mstreal fusionEvaluator::eval(const vector<mstreal>& point, Vector& grad) {
  eval(point);
  grad = gradient;
  // numerical test of the gradient
  if (false) {
    Vector fdGrad = finiteDifferenceGradient(point, vector<mstreal>(point.size(), 10E-5));
    mstreal diff = (grad - fdGrad).norm();
    if (diff > grad.norm()*10E-5) {
      cout << "comparison FAILED:\n";
      for (int i = 0; i < grad.length(); i++) {
        cout << grad[i] << " " << fdGrad[i] << endl;
      }
      cout << "norm difference: " << diff << " (" << 100.0*diff/(0.5*(grad.norm() + fdGrad.norm())) << " %)" << endl;
    } else {
      cout << "comparison SUCCEEDED\n";
    }
  }
  return score;
}

vector<mstreal> fusionEvaluator::guessPoint() {
  if (initPoint.empty()) {
    eval(vector<mstreal>());
  }
  return initPoint;
}


void fusionEvaluator::resetScore() {
  score = bondPenalty = anglPenalty = dihePenalty = rmsdScore = rmsdTot = 0;
  for (int k = 0; k < gradient.size(); k++) gradient[k] = 0;
}

void fusionEvaluator::scoreIC(const icBound& b) {
  mstreal val = b.getCurrentValue();
  mstreal del = 0.0, f, *comp;
  switch (b.type) {
    case icDihedral: {
      mstreal dmin = CartesianGeometry::angleDiffCCW(b.minVal, val);
      mstreal dmax = CartesianGeometry::angleDiffCCW(b.maxVal, val);
      if (dmin < dmax) {
        // outside of the allowed range (before min and after max), so we need
        // the counter-clockwise distance to min and clockwise distance to max
        // (both as positive numbers). If min ends up being closer, the derivative
        // should go with a minus sign, because as the angle becomes larger (and
        // thus moves towards the min), the energy will decrease.
        dmax = 360 - dmax;
        if (dmax < dmin) del = dmax;
        else del = -dmin;
      }
      f = params.getDihedralFC(); comp = &dihePenalty;
      break;
    }
    case icAngle:
      if (val < b.minVal) { del = val - b.minVal; }
      else if (val > b.maxVal) { del = val - b.maxVal; }
      f = params.getAngleFC(); comp = &anglPenalty;
      break;
    case icBond:
      if (val < b.minVal) { del = val - b.minVal; }
      else if (val > b.maxVal) { del = val - b.maxVal; }
      f = params.getBondFC(); comp = &bondPenalty;
      break;
    case icDistRep:
      if (val < b.minVal) { del = val - b.minVal; }
      f = params.getRepFC(); comp = &bondPenalty;
      break;
    case icDistComp:
      if (val > b.maxVal) { del = val - b.maxVal; }
      f = params.getCompFC(); comp = &bondPenalty;
      break;
    case icBrokenDihedral:
    case icBrokenAngle:
    case icBrokenBond:
      return;
    default:
      MstUtils::error("unknown variable type", "fusionEvaluator::scoreIC");
  }
  if (del != 0) {
    mstreal pen = f * del * del;
    score += pen;
    *comp += pen;

    // update the gradient
    vector<mstreal> innerGradient = b.getCurrentGradient();
    int j = 0;
    for (int i = 0; i < b.atoms.size(); i++) {
      Atom* a = b.atoms[i];
      for (int d = 0; d < 3; d++) {
        map<int, mstreal>& parts = gradOfXYZ.getPartials(a, d);
        for (auto it = parts.begin(); it != parts.end(); ++it) {
          gradient[it->first] += 2 * f * del * innerGradient[j] * it->second;
        }
        j++;
      }
    }
  }
}

void fusionEvaluator::scoreRMSD() {
  RMSDCalculator rms;
  rmsdScore = rmsdTot = 0; int N = 0;
  vector<mstreal> innerGradient;
  vector<mstreal> weights = topo.getFragWeights();
  int df = 3;
  if (params.fragRedundancyWeighting()) { // down-weight RMSD contributions from regions with many overlapped fragments
    map<Residue*, int> numOcc; // proportional to the number of times each residue is overlapped
    mstreal Wo = 0, W = 0;
    for (int i = 0; i < weights.size(); i++) Wo += fabs(weights[i]);
    for (int i = 0; i < topo.numAlignedFrags(); i++) {
      AtomPointerVector& atoms = topo.getAlignedFragFused(i);
      // easier to iterate over atoms, so numOcc[res] is not exactly the number
      // of res' occurances, but proportional to it
      for (int j = 0; j < atoms.size(); j++) numOcc[atoms[j]->getResidue()]++;
    }
    for (int i = 0; i < topo.numAlignedFrags(); i++) {
      int T = 0;
      AtomPointerVector& atoms = topo.getAlignedFragFused(i);
      for (int j = 0; j < atoms.size(); j++) T += numOcc[atoms[j]->getResidue()];
      weights[i] *= 1.0/T;
    }
    for (int i = 0; i < weights.size(); i++) W += fabs(weights[i]);
    // the sum of the weights stays the same, but they get re-normalized based on number of occurances
    for (int i = 0; i < topo.numAlignedFrags(); i++) weights[i] *= Wo/W;
  }
  if (params.normalizeRMSD()) { // score per-fragment residual, not total residual
    int nn = 0;
    for (int i = 0; i < topo.numAlignedFrags(); i++) nn += topo.getAlignedFragFused(i).size();
    for (int i = 0; i < topo.numAlignedFrags(); i++) weights[i] *= (1.0*topo.numMobileAtoms())/nn;
  }

  if (params.adaptiveWeighting()) {
    for (int gi = 0; gi < topo.numUniqueOverlaps(); gi++) {
      int n = topo.numFragsOverlapping(gi);
      if (n <= 0) MstUtils::error("empty overlap type!", "fusionEvaluator::scoreRMSD");
      int L = topo.getAlignedFragFused(topo.getFragOverlapping(gi, 0)).size();
      N += L*n;

      mstreal Z = 0;
      vector<mstreal> w(n, 1.0), r(n, 0.0);
      vector<vector<mstreal> > innerGradients(n, vector<mstreal>(L*df, 0.0));
      vector<mstreal> innerGradientZ(L*df, 0.0);

      // compute partition function and collect RMSD gradients
      for (int i = 0; i < n; i++) {
        int fi = topo.getFragOverlapping(gi, i);
        r[i] = rms.qcpRMSDGrad(topo.getAlignedFragFused(fi), topo.getAlignedFragRef(fi), innerGradients[i]);
      }
      mstreal D = MstUtils::min(r); D = D*D;
      for (int i = 0; i < n; i++) {
        w[i] = exp(-params.adaptiveBeta() * (r[i] * r[i] - D));
        Z += w[i];
      }
      // compute weights and gradient of partition function
      for (int j = 0; j < n; j++) {
        w[j] /= Z;
        rmsdScore += w[j] * r[j] * r[j] * L * weights[j];
        rmsdTot += r[j] * r[j] * L;
        for (int k = 0; k < L*df; k++) {
          innerGradientZ[k] -= 2 * w[j] * params.adaptiveBeta() * r[j] * innerGradients[j][k];
        }
      }
      // cout << "group " << gi << ": " << MstUtils::vecToString(w, ", ") << endl;

      // update gradient of score
      for (int i = 0; i < n; i++) {
        int fi = topo.getFragOverlapping(gi, i);
        for (int k = 0; k < L*df; k++) {
          map<int, mstreal>& parts = gradOfXYZ.getPartials(topo.getAlignedFragFused(fi)[k/df], k%df);
          for (auto it = parts.begin(); it != parts.end(); ++it) {
            gradient[it->first] += L*weights[i]*(w[i]*r[i]*(2*(1 - params.adaptiveBeta()*r[i]*r[i])*innerGradients[i][k] - r[i]*innerGradientZ[k])) * it->second;
          }
        }
      }
    }
  } else {
    for (int i = 0; i < topo.numAlignedFrags(); i++) {
      innerGradient.resize(topo.getAlignedFragFused(i).size()*df, 0.0);
      mstreal r = rms.qcpRMSDGrad(topo.getAlignedFragFused(i), topo.getAlignedFragRef(i), innerGradient);
      mstreal r2 = r*r;
      rmsdScore += weights[i] * r2 * topo.getAlignedFragFused(i).size();
      rmsdTot += r2 * topo.getAlignedFragFused(i).size();
      N += topo.getAlignedFragFused(i).size();

      // gradient of RMSD
      int j = 0; // TODO: do as above, where the coordinate index is directly iterated over
      for (int ai = 0; ai < topo.getAlignedFragFused(i).size(); ai++) {
        for (int d = 0; d < df; d++) {
          map<int, mstreal>& parts = gradOfXYZ.getPartials(topo.getAlignedFragFused(i)[ai], d);
          for (auto it = parts.begin(); it != parts.end(); ++it) {
            gradient[it->first] += weights[i] * 2 * r * topo.getAlignedFragFused(i).size() * innerGradient[j] * it->second;
          }
          j++;
        }
      }
    }
  }

  score += rmsdScore;
  rmsdTot = sqrt(rmsdTot/N);
  if (params.isVerbose()) cout << "rmsdScore = " << rmsdScore << " (overall RMSD " << rmsdTot << "), bond penalty = " << bondPenalty << ",  angle penalty = " << anglPenalty << ", dihedral penalty = " << dihePenalty << endl;
}

fusionOutput fusionEvaluator::getScores() {
  fusionOutput scores(bondPenalty, anglPenalty, dihePenalty, rmsdScore, rmsdTot, score);
  return scores;
}

AtomPointerVector fusionEvaluator::atomInstances(int ri, const string& ai) {
  AtomPointerVector atoms;

  // find all instances of the necessary residue and the necessary atom in it
  vector<Residue*>& resI = topo.getOverlappingResidues(ri);
  for (int i = 0; i < resI.size(); i++) {
    atoms.push_back(resI[i]->findAtom(ai));
  }

  return atoms;
}

mstreal fusionEvaluator::bondInitValue(int ri, int rj, const string& ai, const string& aj, bool doNotAverage) {
  if ((params.getCoorInitType() == fusionParams::coorInitType::meanCoor) || doNotAverage) {
    return (fused.getResidue(ri).findAtom(ai))->distance(fused.getResidue(rj).findAtom(aj));
  }
  return bondInstances(rj, fused.getResidue(ri).findAtom(ai), fused.getResidue(rj).findAtom(aj), false).mean();
  // return bondInstances(ri, rj, ai, aj, false).mean();
}

mstreal fusionEvaluator::angleInitValue(int ri, int rj, int rk, const string& ai, const string& aj, const string& ak) {
  if (params.getCoorInitType() == fusionParams::coorInitType::meanCoor) {
    return (fused.getResidue(ri).findAtom(ai))->angle(fused.getResidue(rj).findAtom(aj), fused.getResidue(rk).findAtom(ak));
  }
  return angleInstances(rk, fused.getResidue(ri).findAtom(ai), fused.getResidue(rj).findAtom(aj), fused.getResidue(rk).findAtom(ak), false).mean();
  // return angleInstances(ri, rj, rk, ai, aj, ak, false).mean();
}

mstreal fusionEvaluator::dihedralInitValue(int ri, int rj, int rk, int rl, const string& ai, const string& aj, const string& ak, const string& al) {
  if (params.getCoorInitType() == fusionParams::coorInitType::meanCoor) {
    return (fused.getResidue(ri).findAtom(ai))->dihedral(fused.getResidue(rj).findAtom(aj), fused.getResidue(rk).findAtom(ak), fused.getResidue(rl).findAtom(al));
  }
  return CartesianGeometry::angleMean(dihedralInstances(rl, fused.getResidue(ri).findAtom(ai), fused.getResidue(rj).findAtom(aj), fused.getResidue(rk).findAtom(ak), fused.getResidue(rl).findAtom(al), false));
  // return CartesianGeometry::angleMean(dihedralInstances(ri, rj, rk, rl, ai, aj, ak, al, false));
}

CartesianPoint fusionEvaluator::bondInstances(int rj, Atom* atomI, Atom* atomJ, bool addToCache) {
  CartesianPoint p;
  bool chainBreak = false;
  int ri = rj + atomI->getResidue()->getResidueIndex() - atomJ->getResidue()->getResidueIndex();
  string ai = atomI->getName(); string aj = atomJ->getName();

  // find Chains that contain all necessary residues
  vector<Residue*>& resI = topo.getOverlappingResidues(ri);
  vector<Residue*>& resJ = topo.getOverlappingResidues(rj);
  map<Structure*, vector<Residue*> > S;
  for (int i = 0; i < resI.size(); i++) S[resI[i]->getStructure()].push_back(resI[i]);
  for (int i = 0; i < resJ.size(); i++) S[resJ[i]->getStructure()].push_back(resJ[i]);
  for (auto it = S.begin(); it != S.end(); it++) {
    vector<Residue*>& residues = it->second;
    if (residues.size() != 2) continue;
    Atom* Ai = residues[0]->findAtom(ai);
    Atom* Aj = residues[1]->findAtom(aj);
    p.push_back(Ai->distance(Aj));
  }
  // if no comboes from the same Structure were found, that means this element crosses
  // a chain boundary, so we should try all combinations of atoms
  if (p.size() == 0) {
    chainBreak = true;
    for (int i = 0; i < resI.size(); i++) {
      Atom* Ai = resI[i]->findAtom(ai);
      for (int j = 0; j < resJ.size(); j++) {
        Atom* Aj = resJ[j]->findAtom(aj);
        p.push_back(Ai->distance(Aj));
      }
    }
  }
  if (addToCache) {
    bounds.push_back(icBound(chainBreak ? icBrokenBond : icBond, MstUtils::min(p), MstUtils::max(p), {atomI, atomJ}, MstUtils::toString(ri) + "-" + ai + " : " + MstUtils::toString(rj) + "-" + aj));
  }

  return p;
}

CartesianPoint fusionEvaluator::angleInstances(int rk, Atom* atomI, Atom* atomJ, Atom* atomK, bool addToCache) {
  CartesianPoint p;
  bool chainBreak = false;
  int ri = rk + atomI->getResidue()->getResidueIndex() - atomK->getResidue()->getResidueIndex();
  int rj = rk + atomJ->getResidue()->getResidueIndex() - atomK->getResidue()->getResidueIndex();
  string ai = atomI->getName(); string aj = atomJ->getName(); string ak = atomK->getName();

  // find Chains that contain all necessary residues
  vector<Residue*>& resI = topo.getOverlappingResidues(ri);
  vector<Residue*>& resJ = topo.getOverlappingResidues(rj);
  vector<Residue*>& resK = topo.getOverlappingResidues(rk);
  map<Structure*, vector<Residue*> > S;
  for (int i = 0; i < resI.size(); i++) S[resI[i]->getStructure()].push_back(resI[i]);
  for (int i = 0; i < resJ.size(); i++) S[resJ[i]->getStructure()].push_back(resJ[i]);
  for (int i = 0; i < resK.size(); i++) S[resK[i]->getStructure()].push_back(resK[i]);
  for (auto it = S.begin(); it != S.end(); it++) {
    vector<Residue*>& residues = it->second;
    if (residues.size() != 3) continue;
    Atom* Ai = residues[0]->findAtom(ai);
    Atom* Aj = residues[1]->findAtom(aj);
    Atom* Ak = residues[2]->findAtom(ak);
    p.push_back(Ai->angle(Aj, Ak));
  }
  // if no comboes from the same Structure were found, that means this element crosses
  // a chain boundary, so we should try all combinations of atoms
  if (p.size() == 0) {
    chainBreak = true;
    for (int i = 0; i < resI.size(); i++) {
      Atom* Ai = resI[i]->findAtom(ai);
      for (int j = 0; j < resJ.size(); j++) {
        Atom* Aj = resJ[j]->findAtom(aj);
        for (int k = 0; k < resK.size(); k++) {
          Atom* Ak = resK[k]->findAtom(ak);
          p.push_back(Ai->angle(Aj, Ak));
        }
      }
    }
  }
  if (addToCache) {
    bounds.push_back(icBound(chainBreak ? icBrokenAngle : icAngle, MstUtils::min(p), MstUtils::max(p), {atomI, atomJ, atomK}));
  }

  return p;
}

CartesianPoint fusionEvaluator::dihedralInstances(int rl, Atom* atomI, Atom* atomJ, Atom* atomK, Atom* atomL, bool addToCache) {
  CartesianPoint p;
  bool chainBreak = false;
  int ri = rl + atomI->getResidue()->getResidueIndex() - atomL->getResidue()->getResidueIndex();
  int rj = rl + atomJ->getResidue()->getResidueIndex() - atomL->getResidue()->getResidueIndex();
  int rk = rl + atomK->getResidue()->getResidueIndex() - atomL->getResidue()->getResidueIndex();
  string ai = atomI->getName(); string aj = atomJ->getName(); string ak = atomK->getName(); string al = atomL->getName();

  // find Chains that contain all necessary residues
  vector<Residue*>& resI = topo.getOverlappingResidues(ri);
  vector<Residue*>& resJ = topo.getOverlappingResidues(rj);
  vector<Residue*>& resK = topo.getOverlappingResidues(rk);
  vector<Residue*>& resL = topo.getOverlappingResidues(rl);
  map<Structure*, vector<Residue*> > S;
  for (int i = 0; i < resI.size(); i++) S[resI[i]->getStructure()].push_back(resI[i]);
  for (int i = 0; i < resJ.size(); i++) S[resJ[i]->getStructure()].push_back(resJ[i]);
  for (int i = 0; i < resK.size(); i++) S[resK[i]->getStructure()].push_back(resK[i]);
  for (int i = 0; i < resL.size(); i++) S[resL[i]->getStructure()].push_back(resL[i]);
  for (auto it = S.begin(); it != S.end(); it++) {
    vector<Residue*>& residues = it->second;
    if (residues.size() != 4) continue;
    Atom* Ai = residues[0]->findAtom(ai);
    Atom* Aj = residues[1]->findAtom(aj);
    Atom* Ak = residues[2]->findAtom(ak);
    Atom* Al = residues[3]->findAtom(al);
    p.push_back(Ai->dihedral(Aj, Ak, Al));
  }
  // if no combos from the same Structure were found, that means this element crosses
  // a chain boundary, so we should try all combinations of atoms
  if (p.size() == 0) {
    chainBreak = true;
    for (int i = 0; i < resI.size(); i++) {
      Atom* Ai = resI[i]->findAtom(ai);
      for (int j = 0; j < resJ.size(); j++) {
        Atom* Aj = resJ[j]->findAtom(aj);
        for (int k = 0; k < resK.size(); k++) {
          Atom* Ak = resK[k]->findAtom(ak);
          for (int l = 0; l < resL.size(); l++) {
            Atom* Al = resL[l]->findAtom(al);
            p.push_back(Ai->dihedral(Aj, Ak, Al));
          }
        }
      }
    }
  }
  if (addToCache) {
    bounds.push_back(icBound(chainBreak ? icBrokenDihedral : icDihedral, CartesianGeometry::angleRange(p), {atomI, atomJ, atomK, atomL}));
  }

  return p;
}

mstreal fusionEvaluator::atomRadius(const Atom& a) {
  if (a.isNamed("N")) {
    return 1.6;
  } else if (a.isNamed("CA")) {
    return 2.3; // between CH1E and CH2E in param 19
  } else if (a.isNamed("C")) {
    return 2.1;
  } else if (a.isNamed("O")) {
    return 1.6;
  } else {
    MstUtils::error("do not know radius for atom name " + a.getName(), "fusionEvaluator::atomRadius(const Atom&)");
  }
  return 0;
}

/* --------- Fuser ----------- */

Structure Fuser::fuse(const fusionTopology& topo, fusionOutput& scores, const fusionParams& params) {
  fusionEvaluator E(topo, params); E.setVerbose(false);
  vector<mstreal> bestSolution; mstreal score, bestScore; int bestAnchor;
  vector<vector<mstreal> > trajectory, bestTrajectory; vector<mstreal> trajScores, bestTrajScores;
  for (int i = 0; i < params.numCycles(); i++) {
    if (i > 0) {
      E.noisifyGuessPoint(0.2);
      E.chooseBuildOrigin(true);
    }
    vector<mstreal> solution;
    if (params.getMinimizerType() == fusionParams::gradDescent) {
      score = Optim::gradDescent(E, solution, params.numIters(), params.errTol(), params.isVerbose());
    } else if (params.getMinimizerType() == fusionParams::conjGrad) {
      score = Optim::conjGradMin(E, solution, params.numIters(), params.errTol(), params.isVerbose());
    } else if (params.getMinimizerType() == fusionParams::langevinDyna) {
      trajectory.clear();
      E.guessPoint(); // fills masses (among other things)
      trajScores = Optim::langevinDynamics(E, CartesianPoint(E.getMasses())*params.massFactor(), params.timeStep(), params.viscosity(), params.thermalEnergy(), params.numIters(), trajectory, params.saveInterval(), true);
      score = trajScores.back();
      solution = trajectory.back();
      if (params.logBaseDefined()) {
        fstream lf;
        MstUtils::openFile(lf, params.getLogBase() + ".dyn.pdb", ios::out);
        for (int k = 0; k < trajectory.size(); k++) {
          E.eval(trajectory[k]);
          lf << "MODEL " << k + 1 << endl;
          E.getAlignedStructure().writePDB(lf);
          lf << "ENDMDL" << endl;
        }
        lf.close();
      }
    } else {
      score = Optim::fminsearch(E, params.numIters(), solution, params.isVerbose());
    }
    if ((i == 0) || (score < bestScore)) {
      bestScore = score; bestSolution = solution; bestAnchor = E.getBuildOrigin();
      bestTrajectory = trajectory; bestTrajScores = trajScores;
    }
    if (params.isVerbose()) { E.setVerbose(true); E.eval(solution); E.setVerbose(false); }
  }

  E.setBuildOrigin(bestAnchor);
  if (params.isVerbose()) {
    cout << "best score = " << bestScore << ":" << endl;
    E.setVerbose(true);
  }
  E.eval(bestSolution);
  scores = E.getScores();
  Structure fused = E.getAlignedStructure();
  scores.setFused(fused);
  if (!bestTrajScores.empty()) {
    for (int i = 0; i < bestTrajectory.size(); i++) {
      E.eval(bestTrajectory[i]);
      Structure snap = E.getAlignedStructure();
      scores.addSnapshot(snap, bestTrajScores[i]);
    }
  }
  return fused;
}

Structure Fuser::fuse(const fusionTopology& topo, const fusionParams& params) {
  fusionOutput scores;
  return fuse(topo, scores, params);
}

Structure Fuser::fuse(const vector<vector<Residue*> >& resTopo, fusionOutput& scores, const vector<int>& fixed, const fusionParams& params) {
  fusionTopology topo(resTopo);
  topo.addFixedPositions(fixed);
  return fuse(topo, scores, params);
}

Structure Fuser::fuse(const vector<vector<Residue*> >& resTopo, const vector<int>& fixed, const fusionParams& params) {
  fusionOutput scores;
  return fuse(resTopo, scores, fixed, params);
}

Structure Fuser::autofuse(const vector<Residue*>& residues, int flexOnlyNearOverlaps, const fusionParams& params) {
  // build a proximity search object
  mstreal closeDist = 2.0, pepBondMax = 2.0, pepBondMin = 0.5, pepBondIdeal = 1.3;
  AtomPointerVector CAs(residues.size(), NULL), Ns(residues.size(), NULL), Cs(residues.size(), NULL);;
  for (int i = 0; i < residues.size(); i++) {
    Ns[i] = residues[i]->findAtom("N");
    CAs[i] = residues[i]->findAtom("CA");
    Cs[i] = residues[i]->findAtom("C");
  }
  ProximitySearch psN(Ns, 2*closeDist);
  ProximitySearch psCA(CAs, 2*closeDist);
  ProximitySearch psC(Cs, 2*closeDist);

  // build topology
  vector<vector<Residue*> > resTopo;
  map<Residue*, int> resTopoIdx; // stores indices into the above topology structure that each residue maps into
  vector<Residue*> unsortedResidues = residues;
  while (!unsortedResidues.empty()) {
    Residue* res = unsortedResidues.back();
    unsortedResidues.pop_back();
    vector<int> bucket = psCA.getPointsWithin(res->findAtom("CA"), 0, closeDist);

    // are any close residues already assigned to a position?
    int pos = -1; mstreal dmin;
    for (int i = 0; i < bucket.size(); i++) {
      if (resTopoIdx.find(residues[bucket[i]]) != resTopoIdx.end()) {
        mstreal dist = res->findAtom("CA")->distance(CAs[bucket[i]]);
        if ((pos < 0) || (dist < dmin)) {
          pos = resTopoIdx[residues[bucket[i]]];
          dmin = dist;
        }
      }
    }
    // if so, take their position assignment
    if (pos >= 0) {
      resTopoIdx[res] = pos;
      resTopo[pos].push_back(res);
    } else {
      // otherwise, make a new position for this residue
      resTopo.push_back(vector<Residue*>());
      resTopo.back().push_back(res);
      resTopoIdx[res] = resTopo.size() - 1;
    }
  }

  if (params.isVerbose()) {
    cout << "autofuser generated topology:" << endl;
    for (int i = 0; i < resTopo.size(); i++) {
      cout << i+1 << ": ";
      for (int j = 0; j < resTopo[i].size(); j++) cout << *(resTopo[i][j]) << " ";
      cout << endl;
    }
  }

  // now order residues N-to-C of the final chain
  vector<int> chain;
  vector<bool> unassigned(resTopo.size(), true);
  while (chain.size() != resTopo.size()) {
    // grab a residue that has not yet been assigned to a chain and put it in a new chain
    vector<int> ntoc;
    for (int i = 0; i < unassigned.size(); i++) {
      if (unassigned[i]) { unassigned[i] = false; ntoc.push_back(i); break; }
    }
    if (params.isVerbose()) cout << "growing a new chain starting from residue " << ntoc[0] << "..." << endl;
    // keep growing the chain by finding bonded neighbors on either end
    for (int nc = 0; nc < 2; nc++) {
      while (ntoc.size() != resTopo.size()) {
        int ri = nc ? ntoc.front() : ntoc.back();
        bool found = false;
        for (int i = 0; i < resTopo[ri].size(); i++) {
          Residue* r = resTopo[ri][i];
          Atom* n = nc ? r->findAtom("N") : r->findAtom("C");
          vector<int> neigh = (nc ? psC.getPointsWithin(n->getCoor(), pepBondMin, pepBondMax) : psN.getPointsWithin(n->getCoor(), pepBondMin, pepBondMax));
          if (neigh.empty()) continue;

          // if more than one residue is within a peptide-bond distance, pick the
          // one with the most close-to-ideal bond length
          int next = neigh[0];
          mstreal btol = fabs(n->distance(residues[next]->findAtom(nc ? "C" : "N")) - pepBondIdeal);
          for (int ii = 0; ii < neigh.size(); ii++) {
            mstreal tol = n->distance(residues[neigh[ii]]->findAtom(nc ? "C" : "N")) - pepBondIdeal;
            if (btol > tol) {
              next = neigh[ii];
              btol = tol;
            }
          }

          // now add neighbor on the correct side of the growing chain
          if (params.isVerbose()) cout << "\tN-to-C bond: " << (nc ? resTopoIdx[residues[next]] : ntoc.back()) + 1 << " -- " << (nc ? ntoc.front() : resTopoIdx[residues[next]]) + 1 << endl;
          if (nc) ntoc.insert(ntoc.begin(), resTopoIdx[residues[next]]);
          else ntoc.push_back(resTopoIdx[residues[next]]);
          unassigned[resTopoIdx[residues[next]]] = false;
          found = true;
          break;
        }
        if (!found) break;
      }
    }
    chain.insert(chain.end(), ntoc.begin(), ntoc.end()); // append the last found chain
  }

  // finally, re-order
  vector<vector<Residue*> > oldResTopo = resTopo;
  for (int i = 0; i < chain.size(); i++) {
    resTopo[i] = oldResTopo[chain[i]];
  }

  // if specified, mark as fixed all positions except those that are close enough
  // to positions with multiple overlapping residues
  vector<int> fixedInTopo;
  if (flexOnlyNearOverlaps >= 0) {
    set<int> flexible;
    for (int i = 0; i < resTopo.size(); i++) {
      // mark as flexible if close enough to a position with overlapping residues
      for (int j = max(i - flexOnlyNearOverlaps, 0); j <= min(i + flexOnlyNearOverlaps, (int) resTopo.size() - 1); j++) {
        if (resTopo[j].size() > 1) {
          flexible.insert(i);
          break;
        }
      }
    }
    for (int i = 0; i < resTopo.size(); i++) {
      if (flexible.find(i) == flexible.end()) fixedInTopo.push_back(i);
    }
  }

  // call regular fuser to do all the work
  if (params.isVerbose()) cout << "Fuser::autofuse => fusing a structure of " << resTopo.size() << " residues, with " << fixedInTopo.size() << " positions fixed" << endl;
  if (resTopo.size() <= fixedInTopo.size()) MstUtils::error("number of fixed residues >= than number of residues in topology!");
  return fuse(resTopo, fixedInTopo, params);
}
