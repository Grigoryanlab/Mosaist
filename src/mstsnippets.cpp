class BallTree {
  public:
    BallTree() { root = NULL; }
    BallTree(const AtomPointerVector& _atoms, bool add = true);
    ~BallTree() { if (root != NULL) delete root; }
    bool addPoint(const CartesianPoint& p);
    CartesianPoint getPoint(int i) { return points[i]; }
    void pointsWithin(const CartesianPoint& p, mstreal dmin, mstreal dmax, vector<int>& within) { root->pointsWithin(p, dmin, dmax, within); }
    friend ostream & operator<<(ostream &_os, const BallTree& tree) {
      _os << "(" << *(tree.root) << ")"; return _os;
    }

  private:
    class BallTreeNode {
      public:
        BallTreeNode() { rad = rad2 = 0.0; hasChildren = false; }
        BallTreeNode(const vector<int>& indices, bool add, BallTree* parent);
        ~BallTreeNode() {
          for (int i = 0; i < subtrees.size(); i++) delete(subtrees[i]);
        }
        bool isPointWithin(const CartesianPoint& p) { return (cen.distance2(p) <= rad2); }
        void pointsWithin(const CartesianPoint& p, mstreal dmin, mstreal dmax, vector<int>& within);
        bool addPoint(const CartesianPoint& p, int tag);
        friend ostream & operator<<(ostream &_os, const BallTreeNode& node) {
          _os << "{" << node.rad << "}, [" << node.points.size() << ", " << node.subpoints.size() << "]: (";
          if (node.hasChildren) {
            for (int i = 0; i < node.subtrees.size(); i++) {
              _os << *(node.subtrees[i]);
              if (i != node.subtrees.size() - 1) _os << " | ";
            }
          }
          _os << ")";
          return _os;
        }

      private:
        vector<BallTreeNode*> subtrees;
        BallTree* tree;
        bool hasChildren;
        mstreal rad, rad2;      // a sphere centered at point cen and with
        CartesianPoint cen;     // radius rad contains all points in the sub-tree
        vector<int> points;     // all points inside the node (but not its sub-trees)
        vector<int> subpoints;  // all points inside the sub-tree of the node
    };

    vector<CartesianPoint> points;
    BallTreeNode* root;
};

BallTree::BallTreeNode::BallTreeNode(const vector<int>& indices, bool add, BallTree* parent) {
  tree = parent;
  if (indices.empty()) { // empty leaf node
    hasChildren = false; rad = rad2 = 0.0; return;
  }

  // compute centroid and bounding sphere
  vector<CartesianPoint> P(indices.size());
  for (int i = 0; i < indices.size(); i++) P[i] = tree->getPoint(indices[i]);
  cen.resize(P[0].size(), 0.0);
  for (int i = 0; i < indices.size(); i++) cen += P[i];
  cen /= indices.size();
  rad = 0.0; mstreal dist;
  for (int i = 0; i < indices.size(); i++) {
    dist = cen.distancenc(P[i]);
    if (dist > rad) rad = dist;
  }
  rad2 = rad*rad;

  if (add) subpoints.insert(subpoints.end(), indices.begin(), indices.end());
  if (indices.size() < 5) { // leaf node
    if (add) points.insert(points.end(), indices.begin(), indices.end());
    hasChildren = false;
    return;
  }

  // split points into two balanced groups via K-means
  vector<vector<int> > clusts = Clusterer::kmeans(P, MstUtils::min(3, (int) P.size()));
  for (int i = 0; i < clusts.size(); i++) {
    vector<int>& clust = clusts[i];
    vector<int> subIndices(clust.size());
    for (int j = 0; j < clust.size(); j++) subIndices[j] = indices[clust[j]];
   subtrees.push_back(new BallTreeNode(subIndices, add, parent));
  }
  hasChildren = add;
}

BallTree::BallTree(const AtomPointerVector& _atoms, bool add) {
  points.resize(_atoms.size());
  vector<int> indices(_atoms.size());
  for (int i = 0; i < _atoms.size(); i++) { points[i] = _atoms[i]; indices[i] = i; }
  root = new BallTreeNode(indices, add, this);
}

void BallTree::BallTreeNode::pointsWithin(const CartesianPoint& p, mstreal dmin, mstreal dmax, vector<int>& within) {
auto begin = chrono::high_resolution_clock::now();
  // unravel recursion by keeping a list of sub-trees to visit
  list<BallTree::BallTreeNode*> toVisit(1, this);
  BallTree::BallTreeNode* subtree;
  int i;
// int c = 0, pc = 0;
  while (!toVisit.empty()) {
    subtree = toVisit.front(); toVisit.pop_front();
// c++;
    // is the entire sphere out of consideration?
    mstreal dist = subtree->cen.distancenc(p);
// pc++;
    if ((dist - subtree->rad > dmax) || (dist + subtree->rad < dmin)) continue;

    // is the entire sphere within consideration?
    if ((dist - subtree->rad >= dmin) && (dist + subtree->rad <= dmax)) {
      within.insert(within.end(), subtree->subpoints.begin(), subtree->subpoints.end());
      continue;
    }

    // if not check points within the sphere
    for (i = 0; i < subtree->points.size(); i++) {
// pc++;
      dist = p.distancenc(tree->getPoint(subtree->points[i]));
      if ((dist <= dmax) && (dist >= dmin)) within.push_back(subtree->points[i]);
    }

    // next check points within children spheres
    if (subtree->hasChildren) {
      toVisit.insert(toVisit.end(), subtree->subtrees.begin(), subtree->subtrees.end());
    }
  }
// auto end = chrono::high_resolution_clock::now();
// int diff = chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
// cout << "\t(num points, num subtrees, num dist calcs, time) " << within.size() << " " << c << " " << pc << " " << diff << endl;
}

bool BallTree::addPoint(const CartesianPoint& p) {
  if (root->addPoint(p, points.size())) {
    points.push_back(p);
    return true;
  }
  return false;
}

bool BallTree::BallTreeNode::addPoint(const CartesianPoint& p, int tag) {
  if (!isPointWithin(p)) return false;
  subpoints.push_back(tag);
  for (int i = 0; i < subtrees.size(); i++) {
    if (subtrees[i]->isPointWithin(p)) {
      subtrees[i]->addPoint(p, tag);
      subpoints.push_back(tag);
      return true;
    }
  }
  points.push_back(tag);
  return true;
}

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Implements the FASSA (FAst Structure Search Algorithm). Options:");
  op.addOption("q", "query PDB file.", true);
  op.addOption("d", "a database file with a list of PDB files.", true);
  op.setOptions(argc, argv);

  if (true) {
    FASST S;
    cout << "Building the database..." << endl;
    S.setQuery(op.getString("q"));
    S.addTargets(MstUtils::fileToArray(op.getString("d")));
    S.setRMSDCutoff(2.0);
    cout << "Searching..." << endl;
    S.search();
    cout << "found " << S.numMatches() << " matches:" << endl;
    set<fasstSolution> matches = S.getMatches();
    for (auto it = matches.begin(); it != matches.end(); ++it) {
      cout << *it << endl;
    }
  } else {
    vector<string> db = MstUtils::fileToArray(op.getString("d"));
    db.push_back(op.getString("q"));
    int N = 100; // number of searches per structure
    vector<vector<mstreal> > dmin(db.size(), vector<mstreal>(N)), dmax(db.size(), vector<mstreal>(N));
    vector<vector<CartesianPoint> > points(db.size(), vector<CartesianPoint>(N));
    vector<Structure> S(db.size());
    for (int i = 0; i < S.size(); i++) {
      S[i].readPDB(db[i], "QUIET");
      mstreal xlo, ylo, zlo, xhi, yhi, zhi;
      ProximitySearch::calculateExtent(S[i], xlo, ylo, zlo, xhi, yhi, zhi);
      mstreal xw = xhi - xlo; mstreal yw = yhi - ylo; mstreal zw = zhi - zlo;
      for (int j = 0; j < N; j++) {
        points[i][j] = CartesianPoint(xlo + MstUtils::randUnit()*xw, ylo + MstUtils::randUnit()*yw, zlo + MstUtils::randUnit()*zw);
        dmin[i][j] = MstUtils::randUnit()*2.0;
        dmax[i][j] = dmin[i][j] + MstUtils::randUnit()*10.0;
      }
    }

    // first try with ProximitySearch
    int psPrep = 0, psSearch = 0;
    cout << "doing with ProximitySearch..." << endl;
    vector<vector<vector<int> > > psResults(S.size(), vector<vector<int> >(N));
    for (int i = 0; i < S.size(); i++) {
      auto begin = chrono::high_resolution_clock::now();
      ProximitySearch ps(S[i].getAtoms(), 6.0);
      auto end = chrono::high_resolution_clock::now();
      psPrep += chrono::duration_cast<std::chrono::microseconds>(end-begin).count();
      begin = chrono::high_resolution_clock::now();
      for (int j = 0; j < N; j++) {
        psResults[i][j] = ps.getPointsWithin(points[i][j], dmin[i][j], dmax[i][j]);
      }
      end = chrono::high_resolution_clock::now();
      psSearch += chrono::duration_cast<std::chrono::microseconds>(end-begin).count();
    }
    cout << "done with ProximitySearch: " << psPrep << " us of prep, " << psSearch << " us of search" << endl;

    // then with BallTree
    int btPrep = 0, btSearch = 0;
    cout << "doing with BallTree..." << endl;
    vector<vector<vector<int> > > btResults(S.size(), vector<vector<int> >(N));
    for (int i = 0; i < S.size(); i++) {
      auto begin = chrono::high_resolution_clock::now();
      BallTree bt(S[i].getAtoms());
      auto end = chrono::high_resolution_clock::now();
      btPrep += chrono::duration_cast<std::chrono::microseconds>(end-begin).count();
      begin = chrono::high_resolution_clock::now();
      for (int j = 0; j < N; j++) {
        bt.pointsWithin(points[i][j], dmin[i][j], dmax[i][j], btResults[i][j]);
      }
      end = chrono::high_resolution_clock::now();
      btSearch += chrono::duration_cast<std::chrono::microseconds>(end-begin).count();
    }
    cout << "done with BallTree: " << btPrep << " us of prep, " << btSearch << " us of search" << endl;
    cout << "ball-tree search / ps search = " << (1.0*btSearch)/psSearch << endl;

    // make sure results are the same
    bool same = true;
    for (int i = 0; i < S.size(); i++) {
      for (int j = 0; j < N; j++) {
        if (psResults[i][j].size() != btResults[i][j].size()) { same = false; break; }
        set<int> psSet = MstUtils::contents(psResults[i][j]);
        for (int k = 0; k < btResults[i][j].size(); k++) {
          if (psSet.find(btResults[i][j][k]) == psSet.end()) { same = false; break; }
        }
        if (!same) break;
      }
      if (!same) break;
    }
    if (!same) cout << endl << "Results are NOT the same!" << endl;
    else cout << endl << "Results are the same!" << endl;
  }
}
