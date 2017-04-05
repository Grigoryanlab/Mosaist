#include "mstfuser.h"

fusionEvaluator::fusionEvaluator(const vector<vector<Residue*> >& resTopo, vector<int> _fixedResidues, bool _verbose) {
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
  numMobileAtoms = 4*(resTopo.size() - fixedResidues.size());
  anchorRes = (fixedResidues.size() > 0) ? MstUtils::min(fixedResidues) : -1;

  // create room for fused structure
  fused.appendChain(new Chain());
  for (int i = 0; i < resTopo.size(); i++) {
    if (fixed[i]) {
      if (resTopo[i].size() != 1) MstUtils::error("residue index " + MstUtils::toString(i) + " is marked as fixed, but appears to have more than one residue aligned onto it in the topology", "fusionEvaluator::fusionEvaluator");
      fused[0].appendResidue(new Residue(*(resTopo[i][0])));
    } else {
      Residue* res = new Residue();
      fused[0].appendResidue(res);
      res->appendAtom(new Atom(1, "N", 0.0, 0.0, 0.0, 0, 1.0, false));
      res->appendAtom(new Atom(1, "CA", 0.0, 0.0, 0.0, 0, 1.0, false));
      res->appendAtom(new Atom(1, "C", 0.0, 0.0, 0.0, 0, 1.0, false));
      res->appendAtom(new Atom(1, "O", 0.0, 0.0, 0.0, 0, 1.0, false));
    }
  }
  fused.renumber();

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

  // set internal coordinate force constants
  kb = 100;
  ka =  0.02;
  kh = 0.01;
  verbose = _verbose;
  noise = 0; // the first time we provide a guess point, just take means (later can try noisy start points)
}

double fusionEvaluator::eval(const vector<double>& point) {
  bool init = point.empty();
  if (init) initPoint.resize(0);
  real bR = 0.01; real aR = 1.0; real dR = 1.0; // randomness scale factors

  if (!init && (point.size() != numDF())) MstUtils::error("need to place " + MstUtils::toString(numMobileAtoms) + " atoms, " + (isAnchored() ? "with" : "without") + " anchor, and received " + MstUtils::toString(point.size()) + " parameters -- these appear to disagree", "fusionEvaluator::eval");

  // TODO: take care of the case, where there are fixed residues, but the N-terminus is not fixed!!!!!!!!!!!!!!!

  // -- build the fused backbone, atom by atom
  Chain& F = fused[0];
  int k = 0;                 // parameter index
  Atom *pN = NULL, *pCA = NULL, *pC = NULL, *pO = NULL;  // the previous three placed atoms
  // build forward
  for (int i = anchorRes; i < F.residueSize(); i++) {
    Residue& res = F[i];
    Atom* N = res.findAtom("N");
    Atom* CA = res.findAtom("CA");
    Atom* C = res.findAtom("C");
    Atom* O = res.findAtom("O");
    if (!fixed[i]) {
      if ((i == 0) && !isAnchored()) {
        if (init) {
          initPoint.push_back(bondInstances(i, i, "N", "CA").mean() + bR * MstUtils::randUnit() * noise);
          initPoint.push_back(bondInstances(i, i, "CA", "C").mean() + bR * MstUtils::randUnit() * noise);
          initPoint.push_back(angleInstances(i, i, i, "N", "CA", "C").mean() + aR * MstUtils::randUnit() * noise);
        } else {
          double d0 = point[k]; k++;
          double d1 = point[k]; k++;
          double a0 = point[k]; k++;
          N->setCoor(0.0, 0.0, 0.0);
          CA->setCoor(d0, 0.0, 0.0);
          C->setCoor(d0 - d1*cos(a0), d1*sin(a0), 0.0);
        }
      } else {
        // place N relative to pN, pCA, pC
        if (init) {
          initPoint.push_back(bondInstances(i-1, i, "C", "N").mean() + bR * MstUtils::randUnit() * noise);
          initPoint.push_back(angleInstances(i-1, i-1, i, "CA", "C", "N").mean() + aR * MstUtils::randUnit() * noise);
          initPoint.push_back(CartesianGeometry::angleMean(dihedralInstances(i-1, i-1, i-1, i, "N", "CA", "C", "N")) + dR * MstUtils::randUnit() * noise);
        } else {
          N->build(pC, pCA, pN, point[k], point[k+1], point[k+2]); k += 3;
        }

        // place CA relative to pCA, pC, N
        if (init) {
          initPoint.push_back(bondInstances(i, i, "N", "CA").mean() + bR * MstUtils::randUnit() * noise);
          initPoint.push_back(angleInstances(i-1, i, i, "C", "N", "CA").mean() + aR * MstUtils::randUnit() * noise);
          initPoint.push_back(CartesianGeometry::angleMean(dihedralInstances(i-1, i-1, i, i, "CA", "C", "N", "CA")) + dR * MstUtils::randUnit() * noise);
        } else {
          CA->build(N, pC, pCA, point[k], point[k+1], point[k+2]); k += 3;
        }

        // place C relative to pC, N, CA
        if (init) {
          initPoint.push_back(bondInstances(i, i, "CA", "C").mean() + bR * MstUtils::randUnit() * noise);
          initPoint.push_back(angleInstances(i, i, i, "N", "CA", "C").mean() + aR * MstUtils::randUnit() * noise);
          initPoint.push_back(CartesianGeometry::angleMean(dihedralInstances(i-1, i, i, i, "C", "N", "CA", "C")) + dR * MstUtils::randUnit() * noise);
        } else {
          C->build(CA, N, pC, point[k], point[k+1], point[k+2]); k += 3;
        }

        // if this is the last residue, place the O relative to N-CA-C
        if (i == F.residueSize() - 1) {
          if (init) {
            initPoint.push_back(bondInstances(i, i, "C", "O").mean() + bR * MstUtils::randUnit() * noise);
            initPoint.push_back(angleInstances(i, i, i, "CA", "C", "O").mean() + aR * MstUtils::randUnit() * noise);
            initPoint.push_back(CartesianGeometry::angleMean(dihedralInstances(i, i, i, i, "N", "CA", "C", "O")) + dR * MstUtils::randUnit() * noise);
          } else {
            O->build(C, CA, N, point[k], point[k+1], point[k+2]); k += 3;
          }
        }
      }
    }
    // place previous O relative to pCA, N, pC (an improper)
    if ((i > 0) && !fixed[i-1]) {
      if (init) {
        initPoint.push_back(bondInstances(i-1, i-1, "C", "O").mean() + bR * MstUtils::randUnit() * noise);
        initPoint.push_back(angleInstances(i, i-1, i-1, "N", "C", "O").mean() + aR * MstUtils::randUnit() * noise);
        initPoint.push_back(CartesianGeometry::angleMean(dihedralInstances(i-1, i, i-1, i-1, "CA", "N", "C", "O")) + dR * MstUtils::randUnit() * noise);
      } else {
        pO->build(pC, N, pCA, point[k], point[k+1], point[k+2]); k += 3;
      }
    }
    pN = N; pCA = CA; pC = C; pO = O;
  }

  // build backwards (only would happens when there is an anchor and it is not the 0-th residue)
  for (int i = anchorRes; i >= 0; i--) {
    Residue& res = F[i];
    Atom* N = res.findAtom("N");
    Atom* CA = res.findAtom("CA");
    Atom* C = res.findAtom("C");
    Atom* O = res.findAtom("O");
    if (!fixed[i]) {
      // place C relative to pC, pCA, pN
      if (init) {
        initPoint.push_back(bondInstances(i+1, i, "N", "C").mean() + bR * MstUtils::randUnit() * noise);
        initPoint.push_back(angleInstances(i+1, i+1, i, "CA", "N", "C").mean() + aR * MstUtils::randUnit() * noise);
        initPoint.push_back(CartesianGeometry::angleMean(dihedralInstances(i+1, i+1, i+1, i, "C", "CA", "N", "C")) + dR * MstUtils::randUnit() * noise);
      } else {
        C->build(pN, pCA, pC, point[k], point[k+1], point[k+2]); k += 3;
      }

      // place CA relative to pCA, pN, C
      if (init) {
        initPoint.push_back(bondInstances(i, i, "C", "CA").mean() + bR * MstUtils::randUnit() * noise);
        initPoint.push_back(angleInstances(i+1, i, i, "N", "C", "CA").mean() + aR * MstUtils::randUnit() * noise);
        initPoint.push_back(CartesianGeometry::angleMean(dihedralInstances(i+1, i+1, i, i, "CA", "N", "C", "CA")) + dR * MstUtils::randUnit() * noise);
      } else {
        CA->build(C, pN, pCA, point[k], point[k+1], point[k+2]); k += 3;
      }

      // place N relative to pN, C, CA
      if (init) {
        initPoint.push_back(bondInstances(i, i, "CA", "N").mean() + bR * MstUtils::randUnit() * noise);
        initPoint.push_back(angleInstances(i, i, i, "C", "CA", "N").mean() + aR * MstUtils::randUnit() * noise);
        initPoint.push_back(CartesianGeometry::angleMean(dihedralInstances(i+1, i, i, i, "N", "C", "CA", "N")) + dR * MstUtils::randUnit() * noise);
      } else {
        N->build(CA, C, pN, point[k], point[k+1], point[k+2]); k += 3;
      }

      // place O relative to pCA, pN, C (an improper), except for the residue
      // proceeding the anchor, since that residue's O was placed in forward building
      if (i != anchorRes - 1) {
        if (init) {
          initPoint.push_back(bondInstances(i, i, "C", "O").mean() + bR * MstUtils::randUnit() * noise);
          initPoint.push_back(angleInstances(i+1, i, i, "N", "C", "O").mean() + aR * MstUtils::randUnit() * noise);
          initPoint.push_back(CartesianGeometry::angleMean(dihedralInstances(i, i+1, i, i, "CA", "N", "C", "O")) + dR * MstUtils::randUnit() * noise);
        } else {
          O->build(C, pN, CA, point[k], point[k+1], point[k+2]); k += 3;
        }
      }

    }
    pN = N; pCA = CA; pC = C; pO = O;
  }

  // compute penalty for out-of-range parameters (based on the built struct)
  if (verbose) cout << endl << "evaluating:" << endl;
  double penalty = 0;
  k = 0;
  for (int i = 0; i < F.residueSize(); i++) {
    Residue& res = F[i];
    Atom* N = res.findAtom("N");
    Atom* CA = res.findAtom("CA");
    Atom* C = res.findAtom("C");
    Atom* O = res.findAtom("O");
    // if residue not fixed, evaluate its internal coordinates
    if (!fixed[i]) {
      if (init) {
        bondInstances(i, i, "N", "CA", true);
        bondInstances(i, i, "CA", "C", true);
        bondInstances(i, i, "C", "O", true);
        angleInstances(i, i, i, "N", "CA", "C", true);
        // if O was placed relative to this residue (as opposed to relative to needing the next residue)
        if (i == F.residueSize() - 1) {
          angleInstances(i, i, i, "CA", "C", "O", true);
          dihedralInstances(i, i, i, i, "N", "CA", "C", "O", true);
        }
      } else {
        penalty += harmonicPenalty(N->distance(CA), bounds[k], kb); k++;
        penalty += harmonicPenalty(CA->distance(C), bounds[k], kb); k++;
        penalty += harmonicPenalty(C->distance(O), bounds[k], kb); k++;
        penalty += harmonicPenalty(N->angle(CA, C), bounds[k], ka); k++;
        // if O was placed relative to this residue (as opposed to relative to needing the next residue)
        if (i == F.residueSize() - 1) {
          penalty += harmonicPenalty(CA->angle(C, O), bounds[k], ka); k++;
          penalty += harmonicPenalty(N->dihedral(CA, C, O), bounds[k], kh, true); k++;
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
        // if previous O was placed relative to this resudue
        if (i < F.residueSize() - 1) {
          angleInstances(i-1, i-1, i, "O", "C", "N", true);
          dihedralInstances(i-1, i, i-1, i-1, "CA", "N", "C", "O", true);
        }
      } else {
        penalty += harmonicPenalty(pC->distance(N), bounds[k], kb); k++;
        penalty += harmonicPenalty(pCA->angle(pC, N), bounds[k], ka); k++;
        penalty += harmonicPenalty(pN->dihedral(pCA, pC, N), bounds[k], kh, true); k++;
        penalty += harmonicPenalty(pC->angle(N, CA), bounds[k], ka); k++;
        penalty += harmonicPenalty(pCA->dihedral(pC, N, CA), bounds[k], kh, true); k++;
        penalty += harmonicPenalty(pC->dihedral(N, CA, C), bounds[k], kh, true); k++;
        // if previous O was placed relative to this resudue
        if (i < F.residueSize() - 1) {
          penalty += harmonicPenalty(pO->angle(pC, N), bounds[k], ka); k++;
          penalty += harmonicPenalty(pCA->dihedral(N, pC, pO), bounds[k], kh, true); k++;
        }
      }
    }
    pN = N; pCA = CA; pC = C; pO = O;
  }

  // finally, compute best-fit RMSD of individual fragments onto the built structure
  double rmsdScore = 0;
  if (!init) {
    RMSDCalculator rms;
    for (int i = 0; i < alignedFrags.size(); i++) {
      rmsdScore += rms.bestRMSD(alignedFrags[i].second, alignedFrags[i].first);
    }
  }
  if (verbose) cout << "rmsdScore = " << rmsdScore << ", penalty = " << penalty << ", EVAL DONE" << endl;

  return rmsdScore + penalty;
}

CartesianPoint fusionEvaluator::bondInstances(int ri, int rj, const string& ai, const string& aj, bool addToCache) {
  CartesianPoint p;
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
  if (addToCache) {
    pair<double, double> b(MstUtils::min(p), MstUtils::max(p));
    bounds.push_back(b);
  }

  return p;
}

CartesianPoint fusionEvaluator::angleInstances(int ri, int rj, int rk, const string& ai, const string& aj, const string& ak, bool addToCache) {
  CartesianPoint p;
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
  if (addToCache) {
    pair<double, double> b(MstUtils::min(p), MstUtils::max(p));
    bounds.push_back(b);
  }

  return p;
}

CartesianPoint fusionEvaluator::dihedralInstances(int ri, int rj, int rk, int rl, const string& ai, const string& aj, const string& ak, const string& al, bool addToCache) {
  CartesianPoint p;
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
  if (addToCache) {
    bounds.push_back(CartesianGeometry::angleRange(p));
  }

  return p;
}

double fusionEvaluator::harmonicPenalty(double val, const pair<double, double>& minmax, double K, bool angular) {
  if (angular) {
    if (CartesianGeometry::angleDiffCCW(minmax.first, val) > CartesianGeometry::angleDiffCCW(minmax.second, val)) return 0;
    double dx2 = MstUtils::min(pow(CartesianGeometry::angleDiff(minmax.first, val), 2), pow(CartesianGeometry::angleDiff(minmax.second, val), 2));
    if (verbose) cout << "dihe value = " << val << ", while interval [" << minmax.first << ", " << minmax.second << "], penalty = " << K * dx2 << endl;
    return K * dx2;
  } else {
    double pen = 0;
    if (val < minmax.first) { pen = K * (val - minmax.first) * (val - minmax.first); }
    else if (val > minmax.second) { pen = K * (val - minmax.second) * (val - minmax.second); }
    if (verbose) cout << "value = " << val << ", while interval [" << minmax.first << ", " << minmax.second << "], penalty = " << pen << endl;
    return pen;
  }
}
