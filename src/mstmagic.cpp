#include "mstmagic.h"

vector<AtomPointerVector> TERMUtils::mostDesignableFragments(Structure& C, vector<Structure*>& TERMs, int n, CartesianPoint* cen, CartesianPoint* ext, string outBase, bool verb) {
  int N = 1000;      // examine top this many most promissing points to expand into regions
  mstreal dcut = 2.0;   // inter-atomic distance cutoff for counting neighbors

  // this ProximitySearch object is for checking which residues are already contained within the central TERM (we want to bring new segments)
  AtomPointerVector centAtoms = C.getAtoms();
  ProximitySearch centTERM(centAtoms, (mstreal) 1.0, true, NULL);

  // this ProximitySearch will store all residues brought in by the central TERM
  ProximitySearch* PS;
  if (cen == NULL) {
    // by default, look in a pad x pad x pad box around the center of the central TERM
    mstreal pad = 100;
    PS = new ProximitySearch(centAtoms, (mstreal) 1.0, false, NULL, pad/2);
  } else {
    PS = new ProximitySearch((*cen)[0] - (*ext)[0], (*cen)[1] - (*ext)[1], (*cen)[2] - (*ext)[2], (*cen)[0] + (*ext)[0], (*cen)[1] + (*ext)[1], (*cen)[2] + (*ext)[2]);
  }

  // now add backbone atoms in all overlapping motifs to this ProximitySearch object
  vector<vector<int> > pointSource;
  for (int k = 1; k < TERMs.size(); k++) {
    Structure& S = *(TERMs[k]);
    for (int i = 0; i < S.chainSize(); i++) {
      Chain& chain = S[i];
      for (int j = 0; j < chain.residueSize(); j++) {
        Residue& res = chain[j];
        Atom* a = res.findAtom("CA");
        if (PS->isPointWithinGrid(a)) {
          if (centTERM.pointsWithin(a, 0, dcut)) continue; // don't want to find commonalities among residues in the overlap regions
          PS->addPoint(a, pointSource.size());
          vector<int> source(3, 0); source[0] = k; source[1] = i; source[2] = j;
          pointSource.push_back(source);
        }
      }
    }
  }
  if (verb) cout << "added a total of " << PS->pointSize() << " points to ProximitySearch\n";

  // now find the number of close neighbors for each point
  vector<int> nn(PS->pointSize(), 0);
  for (int i = 0; i < PS->pointSize(); i++) {
    vector<int> close;
    PS->pointsWithin(PS->getPoint(i), 0.0, dcut, &close);
    nn[i] = close.size();
  }

  // sort points by the number of neighbors
  vector<int> sortedIndex = MstUtils::sortIndices(nn, true);

  // for the most promising points, try to expand regions around them
  N = min(N, (int) sortedIndex.size());
  if (verb) cout << "looking at the top " << N << " most promising points:\n";
  vector<double> regScores(N, 0); vector<AtomPointerVector> reg(N);
  for (int i = 0; i < N; i++) {
    vector<int> close;
    PS->pointsWithin(PS->getPoint(sortedIndex[i]), 0.0, dcut, &close);
    vector<int> source = pointSource[sortedIndex[i]];

    // walk in both directions from the central point
    Structure& S = *(TERMs[source[0]]);
    Chain& chain = S[source[1]];
    int L, U, ri;
    for (int s = -1; s <= 1; s += 2) {
      vector<int> inRegister = close; // the set of points, whose corresponding structures continue to be in register with the structure of the current promising point
      for (ri = source[2] + s; (ri >= 0) && (ri < chain.residueSize()); ri += s) {
        Residue& res = chain[ri];
        Atom* a = res.findAtom("CA");
        vector<int> closeCur;
        PS->pointsWithin(a, 0.0, dcut, &closeCur);
        map<int, vector<int> > closeCurTERM;
        for (int ii = 0; ii < closeCur.size(); ii++) {
          int t = pointSource[closeCur[ii]][0]; // which TERM this point belongs to
          closeCurTERM[t].push_back(closeCur[ii]);
        }

        // figure out how many of the previous are still in register at this point
        vector<int> surviving;
        for (int ii = 0; ii < inRegister.size(); ii++) {
          int p = inRegister[ii];
          int t = pointSource[p][0];
          if (closeCurTERM.find(t) == closeCurTERM.end()) continue;
          for (int k = 0; k < closeCurTERM[t].size(); k++) {
            int p1 = closeCurTERM[t][k];
            if ((pointSource[p][1] == pointSource[p1][1]) && (pointSource[p][2] + s == pointSource[p1][2])) {
              surviving.push_back(p);
            }
          }
        }
        if (surviving.size() * 1.0 / close.size() < 0.1) break;
        inRegister = surviving;
      }
      if (s < 0) L = ri + 1;
      else U = ri - 1;
      regScores[i] += fabs(1.0*(ri - source[2]))*inRegister.size()*1.0;
    }
    if (verb) cout << "promising point " << i << ", with " << close.size() << " close neighbors expanded to a region of " << U - L + 1 << " residues, score = " << regScores[i] << "..." << endl;

    // cut out and output promising region
    for (ri = L; ri <= U; ri++) {
      Residue& res = chain[ri];
      for (int ai = 0; ai < res.atomSize(); ai++) {
        reg[i].push_back(&(res[ai]));
      }
    }
  }

  // select the best of the promising regions, forcing diversity
  vector<int> bestRegs = MstUtils::sortIndices(regScores);
  map<int, bool> chosen;
  vector<AtomPointerVector> ret;
  int k = 0;
  for (int i = 0; i < N; i++) {
    int pi = sortedIndex[bestRegs[i]]; // index of the seed point
    // make sure the seed of this region is not too close to that of a previously chosen region
    vector<int> close;
    PS->pointsWithin(PS->getPoint(pi), 0.0, dcut, &close);
    bool redundant = false;
    for (int ii = 0; ii < close.size(); ii++) {
      if (chosen.find(close[ii]) != chosen.end()) {
        redundant = true; break;
      }
    }
    if (redundant) continue;
    chosen[pi] = true;
    ret.push_back(reg[bestRegs[i]]);

    if (!outBase.empty()) {
      Structure regS; regS.addAtoms(reg[bestRegs[i]]);
      regS.writePDB(outBase + MstUtils::toString(i) + ".pdb");
    }
    k++;
    if (k >= n) break;
  }
  delete PS;

  return ret;
}
