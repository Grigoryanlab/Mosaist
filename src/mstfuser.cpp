#include "mstfuser.h"

/* --------- fusionEvaluator ----------- */

fusionEvaluator::fusionEvaluator(const vector<vector<Residue*> >& resTopo, vector<int> _fixedResidues, const fusionParams& _params) {
  params = _params;
  kb = 10;            // internal coordinate force constants
  ka =  0.02;
  kh = 0.001;

  vector<string> bba = {"N", "CA", "C", "O"};
  // copy overlapping residues
  overlappingResidues = resTopo;

  // mark fixed residues
  fixedResidues = _fixedResidues;
  fixed.resize(resTopo.size(), false);
  for (int i = 0; i < fixedResidues.size(); i++) {
    if (fixedResidues[i] >= fixed.size())
      MstUtils::error("expected fused structure with " + MstUtils::toString(resTopo.size()) +
                      " residues, but index " + MstUtils::toString(fixedResidues[i]) +
                      " is out of range", "fusionEvaluator::fusionEvaluator");
      fixed[fixedResidues[i]] = true;
  }
  numMobileAtoms = bba.size()*(resTopo.size() - fixedResidues.size());
  buildOriginRes = (fixedResidues.size() > 0) ? MstUtils::min(fixedResidues) : 0;

  // create room for fused structure (initialize with average coordinates)
  fused.appendChain(new Chain());
  for (int i = 0; i < resTopo.size(); i++) {
    if ((fixed[i]) && (resTopo[i].size() != 1)) MstUtils::error("position index " + MstUtils::toString(i) + " is marked as fixed, but appears to have more than one residue aligned onto it in the topology", "fusionEvaluator::fusionEvaluator");
    if (resTopo[i].size() == 0) MstUtils::error("position index " + MstUtils::toString(i) + " has not overlapping residues in the specified topology", "fusionEvaluator::fusionEvaluator");
    Residue* res = new Residue(resTopo[i][0]->getName(), 1);
    fused[0].appendResidue(res);
    for (int j = 0; j < bba.size(); j++) {
      CartesianPoint m = atomInstances(i, bba[j]).getGeometricCenter();
      res->appendAtom(new Atom(1, bba[j], m.getX(), m.getY(), m.getZ(), 0, 1.0, false));
    }
  }
  fused.renumber();
  guess = fused; // save initial "averaged" guess for later alignment (easier to visualize output)

  /* Initialize alignedFrags. frags[C] designates a segment that moves
   * together (i.e., part of Chain pointed to by C), by storing which
   * Residues the segment has and which residue index, in the full fused
   * structure, they correspond to. */
  map<Structure*, vector<pair<Residue*, int> > > frags;
  for (int i = 0; i < resTopo.size(); i++) {
    for (int j = 0; j < resTopo[i].size(); j++) {
      Structure* S = resTopo[i][j]->getStructure();
      frags[S].push_back(pair<Residue*, int>(resTopo[i][j], i));
    }
  }
  Chain& fusedChain = fused[0];
  for (auto it = frags.begin(); it != frags.end(); ++it) {
    vector<pair<Residue*, int> >& residues = it->second;
    AtomPointerVector fusedAtoms, fragAtoms;
    for (int i = 0; i < residues.size(); i++) {
      int ri = residues[i].second;
      Residue& fragRes = *(residues[i].first);
      Residue& fusedRes = fusedChain[ri];
      for (int j = 0; j < fusedRes.atomSize(); j++) {
        fusedAtoms.push_back(&(fusedRes[j]));
        fragAtoms.push_back(fragRes.findAtom(fusedRes[j].getName()));
      }
    }
    alignedFrags.push_back(pair<AtomPointerVector, AtomPointerVector> (fusedAtoms, fragAtoms));
  }

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

Structure fusionEvaluator::getAlignedStructure() {
  Structure aligned = fused;
  RMSDCalculator rc;
  rc.align(aligned.getAtoms(), guess.getAtoms(), aligned);
  return aligned;
}

mstreal fusionEvaluator::eval(const vector<mstreal>& point) {
  bool init = point.empty();
  if (init) initPoint.resize(0);
  mstreal bR = 0.01; mstreal aR = 1.0; mstreal dR = 1.0; mstreal xyzR = 0.01; // randomness scale factors
  if (!init && (point.size() != numDF())) MstUtils::error("need to place " + MstUtils::toString(numMobileAtoms) + " atoms, " + (isAnchored() ? "with" : "without") + " anchor, and received " + MstUtils::toString(point.size()) + " parameters: that appears wrong!", "fusionEvaluator::eval");
  int k = 0;
  Atom *pN = NULL, *pCA = NULL, *pC = NULL, *pO = NULL;
  mstreal noise = params.getNoise();

  Chain& F = fused[0];
  if (params.getOptimCartesian()) {
    for (int i = 0; i < F.residueSize(); i++) {
      if (fixed[i]) continue;
      Residue& res = F[i];
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
          } else {
            res[j][dim] = point[k]; k++;
          }
        }
      }
    }
  } else {
    // -- build the fused backbone, atom by atom
    // build forward
    int startIdx = isAnchored() ? buildOriginRes : 0;
    for (int i = startIdx; i < F.residueSize(); i++) {
      Residue& res = F[i];
      Atom* N = &(res[0]);
      Atom* CA = &(res[1]);
      Atom* C = &(res[2]);
      Atom* O = &(res[3]);
      if (!fixed[i]) {
        if ((i == 0) && !isAnchored()) {
          if (init) {
            initPoint.push_back(bondInitValue(i, i, "N", "CA") + bR * MstUtils::randUnit() * noise);
            initPoint.push_back(bondInitValue(i, i, "CA", "C") + bR * MstUtils::randUnit() * noise);
            initPoint.push_back(angleInitValue(i, i, i, "N", "CA", "C") + aR * MstUtils::randUnit() * noise);
          } else {
            mstreal d0 = point[k]; k++;
            mstreal d1 = point[k]; k++;
            mstreal a0 = point[k]*M_PI/180; k++;
            N->setCoor(0.0, 0.0, 0.0);
            CA->setCoor(d0, 0.0, 0.0);
            C->setCoor(d0 - d1*cos(a0), d1*sin(a0), 0.0);
          }
        } else {
          // place N relative to pN, pCA, pC
          if (init) {
            initPoint.push_back(bondInitValue(i-1, i, "C", "N") + bR * MstUtils::randUnit() * noise);
            initPoint.push_back(angleInitValue(i-1, i-1, i, "CA", "C", "N") + aR * MstUtils::randUnit() * noise);
            initPoint.push_back(dihedralInitValue(i-1, i-1, i-1, i, "N", "CA", "C", "N") + dR * MstUtils::randUnit() * noise);
          } else {
            N->build(pC, pCA, pN, point[k], point[k+1], point[k+2]); k += 3;
          }

          // place CA relative to pCA, pC, N
          if (init) {
            initPoint.push_back(bondInitValue(i, i, "N", "CA") + bR * MstUtils::randUnit() * noise);
            initPoint.push_back(angleInitValue(i-1, i, i, "C", "N", "CA") + aR * MstUtils::randUnit() * noise);
            initPoint.push_back(dihedralInitValue(i-1, i-1, i, i, "CA", "C", "N", "CA") + dR * MstUtils::randUnit() * noise);
          } else {
            CA->build(N, pC, pCA, point[k], point[k+1], point[k+2]); k += 3;
          }

          // place C relative to pC, N, CA
          if (init) {
            initPoint.push_back(bondInitValue(i, i, "CA", "C") + bR * MstUtils::randUnit() * noise);
            initPoint.push_back(angleInitValue(i, i, i, "N", "CA", "C") + aR * MstUtils::randUnit() * noise);
            initPoint.push_back(dihedralInitValue(i-1, i, i, i, "C", "N", "CA", "C") + dR * MstUtils::randUnit() * noise);
          } else {
            C->build(CA, N, pC, point[k], point[k+1], point[k+2]); k += 3;
          }

          // if this is the last residue, place the O relative to N-CA-C
          if (i == F.residueSize() - 1) {
            if (init) {
              initPoint.push_back(bondInitValue(i, i, "C", "O") + bR * MstUtils::randUnit() * noise);
              initPoint.push_back(angleInitValue(i, i, i, "CA", "C", "O") + aR * MstUtils::randUnit() * noise);
              initPoint.push_back(dihedralInitValue(i, i, i, i, "N", "CA", "C", "O") + dR * MstUtils::randUnit() * noise);
            } else {
              O->build(C, CA, N, point[k], point[k+1], point[k+2]); k += 3;
            }
          }
        }
      }
      // place previous O relative to pCA, N, pC (an improper)
      if ((i > startIdx) && !fixed[i-1]) {
        if (init) {
          initPoint.push_back(bondInitValue(i-1, i-1, "C", "O") + bR * MstUtils::randUnit() * noise);
          initPoint.push_back(angleInitValue(i, i-1, i-1, "N", "C", "O") + aR * MstUtils::randUnit() * noise);
          initPoint.push_back(dihedralInitValue(i-1, i, i-1, i-1, "CA", "N", "C", "O") + dR * MstUtils::randUnit() * noise);
        } else {
          pO->build(pC, N, pCA, point[k], point[k+1], point[k+2]); k += 3;
        }
      }
      pN = N; pCA = CA; pC = C; pO = O;
    }

    // build backwards (only would happens when there is an anchor and it is not the 0-th residue)
    for (int i = buildOriginRes - 1; i >= 0; i--) {
      Residue& res = F[i];
      Atom* N = &(res[0]);
      Atom* CA = &(res[1]);
      Atom* C = &(res[2]);
      Atom* O = &(res[3]);
      if (!fixed[i]) {
        // place C relative to pC, pCA, pN
        if (init) {
          initPoint.push_back(bondInitValue(i+1, i, "N", "C") + bR * MstUtils::randUnit() * noise);
          initPoint.push_back(angleInitValue(i+1, i+1, i, "CA", "N", "C") + aR * MstUtils::randUnit() * noise);
          initPoint.push_back(dihedralInitValue(i+1, i+1, i+1, i, "C", "CA", "N", "C") + dR * MstUtils::randUnit() * noise);
        } else {
          C->build(pN, pCA, pC, point[k], point[k+1], point[k+2]); k += 3;
        }

        // place CA relative to pCA, pN, C
        if (init) {
          initPoint.push_back(bondInitValue(i, i, "C", "CA") + bR * MstUtils::randUnit() * noise);
          initPoint.push_back(angleInitValue(i+1, i, i, "N", "C", "CA") + aR * MstUtils::randUnit() * noise);
          initPoint.push_back(dihedralInitValue(i+1, i+1, i, i, "CA", "N", "C", "CA") + dR * MstUtils::randUnit() * noise);
        } else {
          CA->build(C, pN, pCA, point[k], point[k+1], point[k+2]); k += 3;
        }

        // place N relative to pN, C, CA
        if (init) {
          initPoint.push_back(bondInitValue(i, i, "CA", "N") + bR * MstUtils::randUnit() * noise);
          initPoint.push_back(angleInitValue(i, i, i, "C", "CA", "N") + aR * MstUtils::randUnit() * noise);
          initPoint.push_back(dihedralInitValue(i+1, i, i, i, "N", "C", "CA", "N") + dR * MstUtils::randUnit() * noise);
        } else {
          N->build(CA, C, pN, point[k], point[k+1], point[k+2]); k += 3;
        }

        // place O relative to pCA, pN, C (an improper)
        if (init) {
          initPoint.push_back(bondInitValue(i, i, "C", "O") + bR * MstUtils::randUnit() * noise);
          initPoint.push_back(angleInitValue(i+1, i, i, "N", "C", "O") + aR * MstUtils::randUnit() * noise);
          initPoint.push_back(dihedralInitValue(i, i+1, i, i, "CA", "N", "C", "O") + dR * MstUtils::randUnit() * noise);
        } else {
          O->build(C, pN, CA, point[k], point[k+1], point[k+2]); k += 3;
        }
      }
      pN = N; pCA = CA; pC = C; pO = O;
    }
  }

  // compute penalty for out-of-range parameters (based on the built struct)
  mstreal bondPenalty = 0, anglPenalty = 0, dihePenalty = 0;
  k = 0;
  for (int i = 0; i < F.residueSize(); i++) {
    Residue& res = F[i];
    Atom* N = &(res[0]);
    Atom* CA = &(res[1]);
    Atom* C = &(res[2]);
    Atom* O = &(res[3]);
    // if residue not fixed, evaluate its internal coordinates
    if (!fixed[i]) {
      if (init) {
        bondInstances(i, i, "N", "CA", true);
        bondInstances(i, i, "CA", "C", true);
        bondInstances(i, i, "C", "O", true);
        angleInstances(i, i, i, "N", "CA", "C", true);
        // if last residue, constrain O relative to this residue (as opposed to the next one)
        if (i == F.residueSize() - 1) {
          angleInstances(i, i, i, "CA", "C", "O", true);
          dihedralInstances(i, i, i, i, "N", "CA", "C", "O", true);
        }
      } else {
        bondPenalty += harmonicPenalty(N->distance(CA), bounds[k]); k++;
        bondPenalty += harmonicPenalty(CA->distance(C), bounds[k]); k++;
        bondPenalty += harmonicPenalty(C->distance(O), bounds[k]); k++;
        anglPenalty += harmonicPenalty(N->angle(CA, C), bounds[k]); k++;
        // if last residue, constrain O relative to this residue (as opposed to the next one)
        if (i == F.residueSize() - 1) {
          anglPenalty += harmonicPenalty(CA->angle(C, O), bounds[k]); k++;
          dihePenalty += harmonicPenalty(N->dihedral(CA, C, O), bounds[k]); k++;
        }
      }
    }

    // if either the previous residue or the current one is not fixed, evaluate
    // the internal coordinates connecting the two
    if ((i > 0) && (!fixed[i-1] || !fixed[i])) {
      if (init) {
        bondInstances(i-1, i, "C", "N", true);
        angleInstances(i-1, i-1, i, "CA", "C", "N", true);
        dihedralInstances(i-1, i-1, i-1, i, "N", "CA", "C", "N", true);
        angleInstances(i-1, i, i, "C", "N", "CA", true);
        dihedralInstances(i-1, i-1, i, i, "CA", "C", "N", "CA", true);
        dihedralInstances(i-1, i, i, i, "C", "N", "CA", "C", true);
        // use this residue to constrain the previous O
        if (i < F.residueSize() - 1) {
          angleInstances(i-1, i-1, i, "O", "C", "N", true);
          dihedralInstances(i-1, i, i-1, i-1, "CA", "N", "C", "O", true);
        }
      } else {
        bondPenalty += harmonicPenalty(pC->distance(N), bounds[k]); k++;
        anglPenalty += harmonicPenalty(pCA->angle(pC, N), bounds[k]); k++;
        dihePenalty += harmonicPenalty(pN->dihedral(pCA, pC, N), bounds[k]); k++;
        anglPenalty += harmonicPenalty(pC->angle(N, CA), bounds[k]); k++;
        dihePenalty += harmonicPenalty(pCA->dihedral(pC, N, CA), bounds[k]); k++;
        dihePenalty += harmonicPenalty(pC->dihedral(N, CA, C), bounds[k]); k++;
        // use this residue to constrain the previous O
        if (i < F.residueSize() - 1) {
          anglPenalty += harmonicPenalty(pO->angle(pC, N), bounds[k]); k++;
          dihePenalty += harmonicPenalty(pCA->dihedral(N, pC, pO), bounds[k]); k++;
        }
      }
    }
    pN = N; pCA = CA; pC = C; pO = O;
  }

  // finally, compute best-fit RMSD of individual fragments onto the built structure
  mstreal rmsdScore = 0, rmsdTot = 0;
  if (!init) {
    RMSDCalculator rms;
    for (int i = 0; i < alignedFrags.size(); i++) {
      mstreal r = rms.bestRMSD(alignedFrags[i].second, alignedFrags[i].first);
      // rmsdScore += r;
      rmsdScore += r * r * alignedFrags[i].first.size();
      rmsdTot += r;
    }
  }
  if (params.isVerbose()) cout << "rmsdScore = " << rmsdScore << " (total RMSD " << rmsdTot << "), bond penalty = " << bondPenalty << ",  angle penalty = " << anglPenalty << ", dihedral penalty = " << dihePenalty << endl;

  return rmsdScore + bondPenalty + anglPenalty + dihePenalty;
}

AtomPointerVector fusionEvaluator::atomInstances(int ri, const string& ai) {
  AtomPointerVector atoms;

  // find all instances of the necessary residue and the necessary atom in it
  vector<Residue*>& resI = overlappingResidues[ri];
  for (int i = 0; i < resI.size(); i++) {
    atoms.push_back(resI[i]->findAtom(ai));
  }

  return atoms;
}

mstreal fusionEvaluator::bondInitValue(int ri, int rj, const string& ai, const string& aj, bool doNotAverage) {
  if ((params.getCoorInitType() == fusionParams::coorInitType::meanCoor) || doNotAverage) {
    return (fused.getResidue(ri).findAtom(ai))->distance(fused.getResidue(rj).findAtom(aj));
  }
  return bondInstances(ri, rj, ai, aj).mean();
}

mstreal fusionEvaluator::angleInitValue(int ri, int rj, int rk, const string& ai, const string& aj, const string& ak) {
  if (params.getCoorInitType() == fusionParams::coorInitType::meanCoor) {
    return (fused.getResidue(ri).findAtom(ai))->angle(fused.getResidue(rj).findAtom(aj), fused.getResidue(rk).findAtom(ak));
  }
  return angleInstances(ri, rj, rk, ai, aj, ak).mean();
}

mstreal fusionEvaluator::dihedralInitValue(int ri, int rj, int rk, int rl, const string& ai, const string& aj, const string& ak, const string& al) {
  if (params.getCoorInitType() == fusionParams::coorInitType::meanCoor) {
    return (fused.getResidue(ri).findAtom(ai))->dihedral(fused.getResidue(rj).findAtom(aj), fused.getResidue(rk).findAtom(ak), fused.getResidue(rl).findAtom(al));
  }
  return CartesianGeometry::angleMean(dihedralInstances(ri, rj, rk, rl, ai, aj, ak, al));
}

CartesianPoint fusionEvaluator::bondInstances(int ri, int rj, const string& ai, const string& aj, bool addToCache) {
  CartesianPoint p;
  bool chainBreak = false;

  // find Chains that contain all necessary residues
  vector<Residue*>& resI = overlappingResidues[ri];
  vector<Residue*>& resJ = overlappingResidues[rj];
  map<Chain*, vector<Residue*> > S;
  for (int i = 0; i < resI.size(); i++) S[resI[i]->getChain()].push_back(resI[i]);
  for (int i = 0; i < resJ.size(); i++) S[resJ[i]->getChain()].push_back(resJ[i]);
  for (auto it = S.begin(); it != S.end(); it++) {
    vector<Residue*>& residues = it->second;
    if (residues.size() != 2) continue;
    Atom* Ai = residues[0]->findAtom(ai);
    Atom* Aj = residues[1]->findAtom(aj);
    p.push_back(Ai->distance(Aj));
  }
  // if no comboes from the same chain were found, that means this element crosses
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
    bounds.push_back(icBound(chainBreak ? icBrokenBond : icBond, MstUtils::min(p), MstUtils::max(p), MstUtils::toString(ri) + "-" + ai + " : " + MstUtils::toString(rj) + "-" + aj));
  }

  return p;
}

CartesianPoint fusionEvaluator::angleInstances(int ri, int rj, int rk, const string& ai, const string& aj, const string& ak, bool addToCache) {
  CartesianPoint p;
  bool chainBreak = false;

  // find Chains that contain all necessary residues
  vector<Residue*>& resI = overlappingResidues[ri];
  vector<Residue*>& resJ = overlappingResidues[rj];
  vector<Residue*>& resK = overlappingResidues[rk];
  map<Chain*, vector<Residue*> > S;
  for (int i = 0; i < resI.size(); i++) S[resI[i]->getChain()].push_back(resI[i]);
  for (int i = 0; i < resJ.size(); i++) S[resJ[i]->getChain()].push_back(resJ[i]);
  for (int i = 0; i < resK.size(); i++) S[resK[i]->getChain()].push_back(resK[i]);
  for (auto it = S.begin(); it != S.end(); it++) {
    vector<Residue*>& residues = it->second;
    if (residues.size() != 3) continue;
    Atom* Ai = residues[0]->findAtom(ai);
    Atom* Aj = residues[1]->findAtom(aj);
    Atom* Ak = residues[2]->findAtom(ak);
    p.push_back(Ai->angle(Aj, Ak));
  }
  // if no comboes from the same chain were found, that means this element crosses
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
    bounds.push_back(icBound(chainBreak ? icBrokenAngle : icAngle, MstUtils::min(p), MstUtils::max(p)));
  }

  return p;
}

CartesianPoint fusionEvaluator::dihedralInstances(int ri, int rj, int rk, int rl, const string& ai, const string& aj, const string& ak, const string& al, bool addToCache) {
  CartesianPoint p;
  bool chainBreak = false;

  // find Chains that contain all necessary residues
  vector<Residue*>& resI = overlappingResidues[ri];
  vector<Residue*>& resJ = overlappingResidues[rj];
  vector<Residue*>& resK = overlappingResidues[rk];
  vector<Residue*>& resL = overlappingResidues[rl];
  map<Chain*, vector<Residue*> > S;
  for (int i = 0; i < resI.size(); i++) S[resI[i]->getChain()].push_back(resI[i]);
  for (int i = 0; i < resJ.size(); i++) S[resJ[i]->getChain()].push_back(resJ[i]);
  for (int i = 0; i < resK.size(); i++) S[resK[i]->getChain()].push_back(resK[i]);
  for (int i = 0; i < resL.size(); i++) S[resL[i]->getChain()].push_back(resL[i]);
  for (auto it = S.begin(); it != S.end(); it++) {
    vector<Residue*>& residues = it->second;
    if (residues.size() != 4) continue;
    Atom* Ai = residues[0]->findAtom(ai);
    Atom* Aj = residues[1]->findAtom(aj);
    Atom* Ak = residues[2]->findAtom(ak);
    Atom* Al = residues[3]->findAtom(al);
    p.push_back(Ai->dihedral(Aj, Ak, Al));
  }
  // if no combos from the same chain were found, that means this element crosses
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
    bounds.push_back(icBound(chainBreak ? icBrokenDihedral : icDihedral, CartesianGeometry::angleRange(p)));
  }

  return p;
}

mstreal fusionEvaluator::harmonicPenalty(mstreal val, const icBound& b) {
  switch (b.type) {
    case icDihedral: {
      if (CartesianGeometry::angleDiffCCW(b.minVal, val) > CartesianGeometry::angleDiffCCW(b.maxVal, val)) return 0;
      mstreal dx2 = MstUtils::min(pow(CartesianGeometry::angleDiff(b.minVal, val), 2), pow(CartesianGeometry::angleDiff(b.maxVal, val), 2));
      return kh * dx2;
    }
    case icAngle:
    case icBond: {
      mstreal pen = 0;
      mstreal K = (b.type == icBond) ? kb : ka;
      if (val < b.minVal) { pen = K * (val - b.minVal) * (val - b.minVal); }
      else if (val > b.maxVal) { pen = K * (val - b.maxVal) * (val - b.maxVal); }
      return pen;
    }
    case icBrokenDihedral:
    case icBrokenAngle:
    case icBrokenBond:
      return 0;
    default:
      MstUtils::error("uknown variable type", "fusionEvaluator::harmonicPenalty");
  }
}


/* --------- Fuser ----------- */

Structure Fuser::fuse(const vector<vector<Residue*> >& resTopo, const vector<int>& fixed, const fusionParams& params) {
  bool useGradDescent = true;
  fusionEvaluator E(resTopo, fixed, params); E.setVerbose(false);
  vector<mstreal> bestSolution; mstreal score, bestScore;
  if (useGradDescent) {
    bestScore = Optim::gradDescent(E, bestSolution, params.numIters(), params.errTol(), params.isVerbose());
  } else {
    mstreal bestScore = Optim::fminsearch(E, params.numIters(), bestSolution);
  }
  int bestAnchor = E.getBuildOrigin();
  if (params.isVerbose()) { E.setVerbose(true); E.eval(bestSolution); E.setVerbose(false); }
  for (int i = 0; i < params.numCycles() - 1; i++) {
    E.noisifyGuessPoint(0.2);
    vector<mstreal> solution;
    int anchor = E.randomizeBuildOrigin();
    if (useGradDescent) {
      score = Optim::gradDescent(E, solution, params.numIters(), params.errTol(), params.isVerbose());
    } else {
      score = Optim::fminsearch(E, params.numIters(), solution);
    }
    if (score < bestScore) { bestScore = score; bestSolution = solution; bestAnchor = anchor; }
    if (params.isVerbose()) { E.setVerbose(true); E.eval(solution); E.setVerbose(false); }
  }
  E.setBuildOrigin(bestAnchor);
  if (params.isVerbose()) {
    cout << "best score = " << bestScore << ":" << endl;
    E.setVerbose(true);
  }
  E.eval(bestSolution);
  return E.getAlignedStructure();
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

  // vector<bool> counted(residues.size(), false);
  // for (int i = 0; i < residues.size(); i++) {
  //   if (counted[i]) continue;
  //   resTopo.push_back(vector<Residue*>());
  //   vector<int> bucket = psCA.getPointsWithin(CAs[i]->getCoor(), 0, closeDist);
  //   for (int j = 0; j < bucket.size(); j++) {
  //     resTopo.back().push_back(residues[bucket[j]]);
  //     resTopoIdx[residues[bucket[j]]] = resTopo.size() - 1;
  //     counted[bucket[j]] = true;
  //   }
  // }

  // now order residues N-to-C of the final chain
  vector<int> ntoc(1, 0); // put first residue in chain
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
        if (nc) ntoc.insert(ntoc.begin(), resTopoIdx[residues[next]]);
        else ntoc.push_back(resTopoIdx[residues[next]]);
        found = true;
        break;
      }
      if (!found) break;
    }
  }
  MstUtils::assert(ntoc.size() == resTopo.size(), "could not deduce chain connectivity automatically", "Fuser::autofuse");

  // finally, re-order
  vector<vector<Residue*> > oldResTopo = resTopo;
  for (int i = 0; i < ntoc.size(); i++) {
    resTopo[i] = oldResTopo[ntoc[i]];
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
  return fuse(resTopo, fixedInTopo, params);
}
