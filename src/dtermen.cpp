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
  intCut = 0.01; // set to a value over 1.0 to not count sidechain-backbone contacts via interference
  pmSelf = 1;
  pmPair = 1;
  selfResidualPC = selfCorrPC = 1.0;
  selfResidualMinN = 1000;
  selfResidualMaxN = 5000;
  selfCorrMaxCliqueSize = -1;
  selfCorrMinN = 200;
  selfCorrMaxN = 5000;
  pairMinN = 1000;
  pairMaxN = 5000;
  recordData = false;
  homCut = 0.6;
  setAminoAcidMap();
  setEnergyFunction("35");

  // set up FASST base options
  F.setOptions(fasstSearchOptions());
  F.setRedundancyProperty("sim");
  foptsBase = F.options();
}

void dTERMen::setEnergyFunction(const string& ver) {
  efunVer = ver;
  if (efunVer.compare("35") == 0) {
    // this is the default
  } else if (efunVer.compare("35.f") == 0) {
    // do not use interference-based contacts
    intCut = 1.1;
  } else {
    MstUtils::error("unknown energy function version '" + efunVer + "'");
  }
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
    } else if (ents[0].compare("efun") == 0) {
      setEnergyFunction(ents[1]);
    } else if (ents[0].compare("selfCorrMaxCliqueSize") == 0) {
      selfCorrMaxCliqueSize = MstUtils::toInt(ents[1]);
    } else if (ents[0].compare("selfResidualLims") == 0) {
      vector<int> lims = MstUtils::splitToInt(ents[1]);
      if (lims.size() != 2) MstUtils::error("expected two integers in selfResidualLims field", "dTERMen::readConfigFile");
      selfResidualMinN = lims[0];
      selfResidualMaxN = lims[1];
    } else if (ents[0].compare("pairLims") == 0) {
      vector<int> lims = MstUtils::splitToInt(ents[1]);
      if (lims.size() != 2) MstUtils::error("expected two integers in pairLims field", "dTERMen::readConfigFile");
      pairMinN = lims[0];
      pairMaxN = lims[1];
    } else if (ents[0].compare("selfCorrLims") == 0) {
      vector<int> lims = MstUtils::splitToInt(ents[1]);
      if (lims.size() != 2) MstUtils::error("expected two integers in selfCorrLims field", "dTERMen::readConfigFile");
      selfCorrMinN = lims[0];
      selfCorrMaxN = lims[1];
    } else if (ents[0].compare("homCut") == 0) {
      homCut = MstUtils::toReal(ents[1]);
    } else {
      MstUtils::error("unknown parameter name '" + ents[0] + "'", "dTERMen::dTERMen(const string&)");
    }
  }

  if (fasstdbPath.empty()) MstUtils::error("FASST database not defined in configuration file " + configFile, "dTERMen::dTERMen(const string&)");
  F.readDatabase(fasstdbPath, 2);
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

vector<pair<Residue*, Residue*>> dTERMen::getContactsWith(const vector<Residue*>& source, ConFind& C, int type, bool verbose) {
  set<Residue*> sourceSet = MstUtils::contents(source);
  vector<pair<Residue*, Residue*>> conts;
  set<pair<Residue*, Residue*>> contsSet;
  pair<Residue*, Residue*> c;
  contactList contList;
  if (verbose) {
    cout << "identifying contacts with:";
    for (int i = 0; i < source.size(); i++) cout << " " << *(source[i]);
    cout << endl;
  }

  // get all contacts involving the source residues
  for (int cType = 0; cType < 2; cType++) {
    if (cType == 0) contList = C.getContacts(source, cdCut);
    else {
      if (intCut > 1.0) continue;
      contList = C.getInterfering(source, intCut);
    }

    // go through each and insert into list, in the right order, if it qualifies
    contList.sortByDegree();
    for (int i = 0; i < contList.size(); i++) {
      bool isInA = (sourceSet.find(contList.residueA(i)) != sourceSet.end());
      bool isInB = (sourceSet.find(contList.residueB(i)) != sourceSet.end());
      if (((type == 0) && (isInA == isInB)) || ((type == 1) && !(isInA && isInB)) || ((type == 2) && !(isInA || isInB))) continue;

      if (!isInA) c = pair<Residue*, Residue*>(contList.residueB(i), contList.residueA(i));
      else c = pair<Residue*, Residue*>(contList.residueA(i), contList.residueB(i));

      if (verbose) cout << "\t" << (cType ? "interference " : "contact-degree ") << "contact with " << *(c.second) << " (from " << *(c.first) << "); " << contList.degree(i) << endl;
      if (contsSet.find(c) == contsSet.end()) {
        contsSet.insert(c);
        conts.push_back(c);
      }
    }
  }
  if (verbose) cout << "in the end, found " << conts.size() << " contacts" << endl;

  return conts;
}


EnergyTable dTERMen::buildEnergyTable(const vector<Residue*>& variable, const vector<vector<string>>& allowed, const vector<vector<Residue*>>& images, EnergyTable* specTable, const vector<Residue*>& specContext) {
  EnergyTable E;
  if (variable.empty()) return E;
  Structure* S = variable[0]->getStructure();
  set<Residue*> variableSet = MstUtils::contents(variable);
  map<Residue*, string> siteNames;
  for (int i = 0; i < variable.size(); i++) {
    if (variable[i]->getStructure() != S) MstUtils::error("mutable positions belong to different Structure objects", "dTERMen::buildEnergyTable");
    string siteName = variable[i]->getChain()->getID() + "," + MstUtils::toString(variable[i]->getNum());
    if (variable[i]->getIcode() != ' ') siteName += variable[i]->getIcode();
    siteNames[variable[i]] = siteName;
    E.addSite(siteName);
    if (allowed.empty()) {
      for (int k = 0; k < globAlph.size(); k++) E.addToSiteAlphabet(i, SeqTools::idxToTriple(globAlph[k]));
    } else {
      for (int k = 0; k < allowed[i].size(); k++) {
        int aai = aaToIndex(allowed[i][k]);
        if (!aaIndexKnown(aai)) MstUtils::error("error in alphabet specification, amino acid '" + allowed[i][k] + "' at site " + siteName + " not known", "dTERMen::buildEnergyTable");
        E.addToSiteAlphabet(i, indexToResName(aai));
      }
    }
  }
  if (specTable != NULL) *specTable = E;
  ConFind C(&RL, *S); // make a ConFind object that will keep getting reused in energy calculations

  /* If dealing with crystal symmetry, create map:
   * imgToCen[Ri] is the residue in the central unit cell corresponding to residue
   *              Ri in some image.
   * cenToImg[Ri] is the list of image residues, corresponding to the residue Ri
   *              from the central unit cell. Residues are listed in the order
   *              in which images are specified in the images vector.
   * */
   map<Residue*, Residue*> imgToCen;
   map<Residue*, vector<Residue*>> cenToImg;
   if (!images.empty()) {
     for (int mi = 1; mi < images.size(); mi++) {
       if (images[mi].size() != images[0].size()) MstUtils::error("central unit cell has " + MstUtils::toString(images[0].size()) + " residues, while image " + MstUtils::toString(mi) + " has " + MstUtils::toString(images[mi].size()) + " residues", "dTERMen::buildEnergyTable");
       for (int ri = 0; ri < images[mi].size(); ri++) imgToCen[images[mi][ri]] = images[0][ri];
     }
     for (int ri = 0; ri < images[0].size(); ri++) {
       for (int mi = 1; mi < images.size(); mi++) {
         cenToImg[images[0][ri]].push_back(images[mi][ri]);
       }
     }
   }

  // all residues that contact variable positions, but are not variable, are fixed
  set<Residue*> fixedSet; // this set will include fixed residues in image unit cells
  vector<pair<Residue*, Residue*>> conts = getContactsWith(variable, C, 2);
  for (int i = 0; i < conts.size(); i++) {
    Residue* res = conts[i].second;
    if (variableSet.find(res) == variableSet.end()) {
      if ((imgToCen.find(res) == imgToCen.end()) || (variableSet.find(imgToCen[res]) == variableSet.end())) {
        fixedSet.insert(res);
      } else {
        // add to variable set images of variable positions, with which variable
        // positions in the central unit cell interact
        variableSet.insert(res);
      }
    }
  }

  // make sure amino acids at fixed positions are legal
  for (auto it = fixedSet.begin(); it != fixedSet.end(); ++it) {
    Residue* res = *it;
    if (!aaIndexKnown(aaToIndex(res->getName()))) MstUtils::error("fixed position " + MstUtils::toString(*(*it)) + " occupied with unknown amino acid", "dTERMen::buildEnergyTable");
    if (imgToCen.find(res) != imgToCen.end()) {
      Residue* cres = imgToCen[res];
      if (!(res->isNamed(cres->getName()))) MstUtils::error("fixed residue " + MstUtils::toString(*cres) + " and its corresponding image residue " + MstUtils::toString(*res) + " are occupied with different amino acids", "dTERMen::buildEnergyTable");
    }
  }

  // specificity context is a sub-set of the fixed positions
  set<Residue*> specContextSet;
  if (specTable != NULL) {
    if (specContext.empty()) specContextSet = fixedSet;
    else {
      specContextSet = MstUtils::contents(specContext);
      for (Residue* contRes : specContext) {
        if (cenToImg.find(contRes) == cenToImg.end()) continue;
        vector<Residue*> imgResis = cenToImg[contRes];
        specContextSet.insert(imgResis.begin(), imgResis.end());
      }
      specContextSet = MstUtils::contents(MstUtils::setintersect(MstUtils::keys(fixedSet), MstUtils::keys(specContextSet)));
    }
    cout << "---> " << specContextSet.size() << " residues are considered as the fixed context" << endl;
  }

  // record target and design problem info
  targetOrigSeq = Sequence(*S);
  variableResidues.resize(variable.size());
  for (int i = 0; i < variable.size(); i++) variableResidues[i] = variable[i]->getResidueIndex();
  vector<Residue*> allResidues = S->getResidues();
  for (Residue* R : allResidues) {
    targetResidueProperties["phi"].push_back(R->getPhi());
    targetResidueProperties["psi"].push_back(R->getPsi());
    targetResidueProperties["omega"].push_back(R->getOmega());
    targetResidueProperties["env"].push_back(C.getFreedom(R));
  }

  // compute self energies
  for (int i = 0; i < variable.size(); i++) {
    cout << "computing self energy for position " << *(variable[i]) << ", " << i+1 << "/" << variable.size() << endl;
    vector<mstreal> selfE = selfEnergies(variable[i], C, true);
    vector<string> alpha = E.getSiteAlphabet(i);
    for (int k = 0; k < alpha.size(); k++) {
      int aai = aaToIndex(alpha[k]);
      E.setSelfEnergy(i, k, selfE[aai]);
      if (specTable != NULL) specTable->setSelfEnergy(i, k, 0.0);
    }
  }

  // compute pair energies
  for (int i = 0; i < conts.size(); i++) {
    Residue* resA = conts[i].first;
    Residue* resB = conts[i].second;
    cout << "computing pair energy for positions " << *resA << " x " << *resB << ", " << i+1 << "/" << conts.size() << endl;
    vector<vector<mstreal>> pairE = pairEnergies(resA, resB, true);
    int si = E.siteIndex(siteNames[resA]); // residue A will always be a variable one (that's how we constructed conts)
    vector<string> alphaA = E.getSiteAlphabet(si);

    // differentiate based on whether resB is variable or not
    if (variableSet.find(resB) != variableSet.end()) {
      if (imgToCen.find(resB) == imgToCen.end()) {
        // bona fide pair interaction within the central unit cell
        int sj = E.siteIndex(siteNames[resB]);
        vector<string> alphaB = E.getSiteAlphabet(sj);
        for (int a = 0; a < alphaA.size(); a++) {
          int ai = aaToIndex(alphaA[a]);
          for (int b = 0; b < alphaB.size(); b++) {
            int bi = aaToIndex(alphaB[b]);
            E.setPairEnergy(si, sj, a, b, pairE[ai][bi]);
          }
        }
      } else {
        // pair interaction between variable positions across unit cells
        if (imgToCen[resB] == resA) {
          // interaction of a residue with its own image (goes into self)
          for (int a = 0; a < alphaA.size(); a++) {
            int ai = aaToIndex(alphaA[a]);
            E.setSelfEnergy(si, a, E.selfEnergy(si, a) + 0.5*pairE[ai][ai]);
          }
        } else {
          // interaction of a residue with an image of another residues (goes into pair)
          Residue* cenB = imgToCen[resB];
          int sj = E.siteIndex(siteNames[cenB]);
          vector<string> alphaB = E.getSiteAlphabet(sj);
          for (int a = 0; a < alphaA.size(); a++) {
            int ai = aaToIndex(alphaA[a]);
            for (int b = 0; b < alphaB.size(); b++) {
              int bi = aaToIndex(alphaB[b]);
              E.setPairEnergy(si, sj, a, b, 0.5*pairE[ai][bi]);
            }
          }
        }
      }
    } else {
      if (fixedSet.find(resB) == fixedSet.end()) MstUtils::error("expected only variable-variable and variable-fixed interactions", "dTERMen::buildEnergyTable");
      // interaction between a variable and a fixed position (goes into self)
      mstreal sf = (imgToCen.find(resB) == imgToCen.end()) ? 1.0 : 0.5; // if across images, divide by two
      int bi = aaToIndex(resB->getName());
      for (int a = 0; a < alphaA.size(); a++) {
        int ai = aaToIndex(alphaA[a]);
        E.setSelfEnergy(si, a, E.selfEnergy(si, a) + sf*pairE[ai][bi]);
      }

      // if this is an interaction with a fixed-context residue, need to update spec gap table
      if (specContextSet.find(resB) != specContextSet.end()) {
        mstreal sf = (imgToCen.find(resB) == imgToCen.end()) ? 1.0 : 0.5; // if across images, divide by two
        vector<mstreal> contextProbs(globalAlphabetSize(), 0.0);
        for (int biAlt = 0; biAlt < globalAlphabetSize(); biAlt++) contextProbs[biAlt] = backEner(biAlt);
        contextProbs[bi] = INFINITY; // exclude the identity currently present at this position
        dTERMen::enerToProb(contextProbs);
        for (int a = 0; a < alphaA.size(); a++) {
          int ai = aaToIndex(alphaA[a]);
          mstreal meanAltPairEner = 0;
          for (int biAlt = 0; biAlt < globalAlphabetSize(); biAlt++) {
            meanAltPairEner += contextProbs[biAlt] * pairE[ai][biAlt];
          }
          mstreal pairEner = pairE[ai][bi];
          mstreal gap = pairEner - meanAltPairEner;
          specTable->setSelfEnergy(si, a, specTable->selfEnergy(si, a) + sf*gap);
        }
      }
    }
  }

  return E;
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
    if (aaIdx.find(i) == aaIdx.end()) { int idx = aaIdx.size(); aaIdx[i] = idx; }
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

bool dTERMen::aaIndexKnown(int aaIdx) const {
  return aaIdx >= 0;
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
    Sequence S = F.getTargetSequence(ti);
    int N = S.length();

    // compute multiplicity of each residue
    vector<mstreal> mult(N, 1); // multiplicity of each residue in the structure
    map<int, vector<FASST::resAddress>> simsToTarget = F.getResidueRelationships(ti, "sim");
    for (auto it = simsToTarget.begin(); it != simsToTarget.end(); ++it) {
      mult[it->first] += (it->second).size();
    }

    // store all properties
    for (int ri = 0; ri < N; ri++) {
      string aaName = S.getResidue(ri, true);
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
  if (ky < 0) return 0;
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

void dTERMen::writeRecordedData(const string& file) {
  // find any duplicate TERMs
  map<set<int>, int> byResidueSet;
  for (int i = 0; i < data.size(); i++) {
    set<int> resIdxSet = MstUtils::contents(data[i].getResidueIndices());
    if (byResidueSet.find(resIdxSet) == byResidueSet.end()) {
      byResidueSet[resIdxSet] = i;
    } else {
      // if a TERM with the same residues was seen before, keep the one with more matches
      int pi = byResidueSet[resIdxSet];
      if (data[i].numMatches() > data[pi].numMatches()) byResidueSet[resIdxSet] = i;
    }
  }
  vector<int> unique = MstUtils::values(byResidueSet);
  sort(unique.begin(), unique.end());

  // write data about the target
  fstream ofs; MstUtils::openFile(ofs, file, ios::out);
  ofs << targetOrigSeq.toString() << endl;
  ofs << MstUtils::vecToString(variableResidues) << endl;
  for (int i = 0; i < targetOrigSeq.length(); i++) {
    ofs << targetResidueProperties["phi"][i] << " ";
    ofs << targetResidueProperties["psi"][i] << " ";
    ofs << targetResidueProperties["omega"][i] << " ";
    ofs << targetResidueProperties["env"][i] << endl;
  }

  // write data about unique TERMs
  for (int i : unique) {
    ofs << "* TERM " << i << endl;
    ofs << MstUtils::vecToString(data[i].getResidueIndices()) << endl;
    for (int j = 0; j < data[i].numMatches(); j++) {
      const fasstSolution& m = data[i].getMatch(j);
      ofs << F.getMatchSequence(m).toString() << " " << m.getRMSD() << " ";
      ofs << MstUtils::vecToString(F.getResidueProperties(m, "phi")) << " ";
      ofs << MstUtils::vecToString(F.getResidueProperties(m, "psi")) << " ";
      ofs << MstUtils::vecToString(F.getResidueProperties(m, "omega")) << " ";
      ofs << MstUtils::vecToString(F.getResidueProperties(m, "env")) << endl;
    }
  }
  ofs.close();
}

mstreal dTERMen::selfEnergy(Residue* R, const string& aa) {
  return (selfEnergies(R))[dTERMen::aaToIndex(aa.empty() ? R->getName() : aa)];
}

vector<mstreal> dTERMen::selfEnergies(Residue* R, bool verbose) {
  if (R->getStructure() == NULL) MstUtils::error("cannot operate on a disembodied residue!", "dTERMen::selfEnergies(Residue*, bool)");
  ConFind C(&RL, *(R->getStructure()));
  return selfEnergies(R, C, verbose);
}

vector<mstreal> dTERMen::selfEnergies(Residue* R, ConFind& C, bool verbose) {
  auto rmsdCutSelfRes = [](const vector<int>& fragResIdx, const Structure& S) { return RMSDCalculator::rmsdCutoff(fragResIdx, S, 1.0, 20); };
  auto rmsdCutSelfCor = [](const vector<int>& fragResIdx, const Structure& S) { return RMSDCalculator::rmsdCutoff(fragResIdx, S, 1.1, 15); };
  if (R->getStructure() == NULL) MstUtils::error("cannot operate on a disembodied residue!", "dTERMen::selfEnergies(Residue*, ConFind&, bool)");
  Structure& S = *(R->getStructure());

  // -- simple environment components
  if (verbose) cout << "\tdTERMen::selfEnergies -> trivial statistical components..." << endl;
  int naa = globalAlphabetSize();
  CartesianPoint selfE(naa, 0.0);
  for (int aai = 0; aai < naa; aai++) {
    selfE[aai] = backEner(aai) + bbOmegaEner(R->getOmega(), aai) + bbPhiPsiEner(R->getPhi(), R->getPsi(), aai) + envEner(C.getFreedom(R), aai);
  }
  if (verbose) printSelfComponent(selfE, "\t");

  // -- self residual
  if (verbose) cout << "\tdTERMen::selfEnergies -> self residual..." << endl;
  termData sT({R}, pmSelf);
  F.setOptions(foptsBase);
  F.setQuery(sT.getTERM());
  F.setRMSDCutoff(rmsdCutSelfRes(sT.getResidueIndices(), S));
  F.setMinNumMatches(selfResidualMinN);
  F.setMaxNumMatches(selfResidualMaxN);
  sT.setMatches(F.search(), homCut, &F);
  CartesianPoint selfResidual = singleBodyStatEnergy(sT.getMatches(), sT.getCentralResidueIndices()[0], selfResidualPC);
  if (recordData) data.push_back(sT);
  if (verbose) printSelfComponent(selfResidual, "\t");
  selfE += selfResidual;

  // -- self correction
  if ((selfCorrMaxCliqueSize >= 0) && (selfCorrMaxCliqueSize < 2)) return selfE; // if max clique size is less than 2, then there is effectively no self residual
  if (verbose) cout << "\tdTERMen::selfEnergies -> self correction..." << endl;

  // -- get contacts
  vector<pair<Residue*, Residue*>> conts = getContactsWith({R}, C, 0, verbose);
  vector<Residue*> contResidues(conts.size(), NULL);
  for (int i = 0; i < conts.size(); i++) contResidues[i] = conts[i].second;

  // consider each contacting residue with the central one and see if there are
  // enough matches. If so, create a clique to be grown later. If not, still
  // create a clique (with number of matches), which will not be grown.
  vector<termData> finalCliques;
  map<Residue*, termData> cliquesToGrow;
  for (int i = 0; i < contResidues.size(); i++) {
    if (verbose) cout << "\t\tdTERMen::selfEnergies -> seed clique with contact " << *(contResidues[i]) << "..." << endl;
    termData c({R, contResidues[i]}, pmSelf);
    F.setOptions(foptsBase);
    F.setQuery(c.getTERM());
    F.setRMSDCutoff(rmsdCutSelfCor(c.getResidueIndices(), S));
    F.setMinNumMatches(selfCorrMinN);
    F.setMaxNumMatches(selfCorrMaxN);
    c.setMatches(F.search(), homCut, &F);
    if ((c.numMatches() < selfCorrMinN) || (c.getMatch(selfCorrMinN - 1).getRMSD()) > F.getRMSDCutoff()) { finalCliques.push_back(c); }
    else { cliquesToGrow[contResidues[i]] = c; }
  }

  // for those contacts with sufficient matches, try to combine into larger cliques
  while (!cliquesToGrow.empty()) {
    // sort remaining contacting residues by decreasing number of matches of the
    // clique comprising the residue and the central residue
    vector<Residue*> remConts = MstUtils::keys(cliquesToGrow);
    sort(remConts.begin(), remConts.end(), [&cliquesToGrow](Residue* i, Residue* j) { return cliquesToGrow[i].numMatches() > cliquesToGrow[j].numMatches(); });

    // pick the one with highest number of contacts, and try to grow the clique
    termData parentClique = cliquesToGrow[remConts[0]];
    remConts.erase(remConts.begin());
    termData grownClique = parentClique;
    if (verbose) cout << "\t\tdTERMen::selfEnergies -> will try to grow [" << parentClique.toString() << "]..." << endl;
    while (!remConts.empty()) {
      if ((selfCorrMaxCliqueSize >= 0) && (parentClique.numCentralResidues() >= selfCorrMaxCliqueSize)) break;
      // try to add every remaining contact
      for (int j = 0; j < remConts.size(); j++) {
        if (verbose) cout << "\t\t\tdTERMen::selfEnergies -> trying to add " << *(remConts[j]) << "..." << endl;
        termData newClique = parentClique;
        newClique.addCentralResidue(remConts[j], pmSelf);
        F.setOptions(foptsBase);
        F.setQuery(newClique.getTERM());
        F.setRMSDCutoff(rmsdCutSelfCor(newClique.getResidueIndices(), S));
        F.setMaxNumMatches(selfCorrMaxN);
        newClique.setMatches(F.search(), homCut, &F);
        if ((j == 0) || (newClique.numMatches() > grownClique.numMatches())) {
          if (verbose) cout << "\t\t\t\tdTERMen::selfEnergies -> new best" << endl;
          grownClique = newClique;
        }
      }
      // EITHER a sufficient number of matches OR not too many fewer than for the
      // parent clique is the logic we had implemented in original dTERMen
      if ((grownClique.numMatches() >= selfCorrMinN) || (grownClique.numMatches() > 0.8*parentClique.numMatches())) {
        remConts = MstUtils::setdiff(remConts, {grownClique.getCentralResidues().back()});
        parentClique = grownClique;
        if (verbose) cout << "\t\t\tdTERMen::selfEnergies -> chose to add " << *(grownClique.getCentralResidues().back()) << "..." << endl;
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
    if (recordData) data.push_back(finalCliques[i]);
    CartesianPoint cliqueDelta = singleBodyStatEnergy(finalCliques[i].getMatches(), finalCliques[i].getCentralResidueIndices()[0], selfCorrPC);
    if (verbose) printSelfComponent(cliqueDelta, "\t\t\t");
    selfE += cliqueDelta;
  }

  return selfE;
}

mstreal dTERMen::pairEnergy(Residue* Ri, Residue* Rj, const string& aai, const string& aaj) {
  return (pairEnergies(Ri, Rj))[dTERMen::aaToIndex(aai.empty() ? Ri->getName() : aai)]
                               [dTERMen::aaToIndex(aaj.empty() ? Rj->getName() : aaj)];
}

vector<vector<mstreal>> dTERMen::pairEnergies(Residue* Ri, Residue* Rj, bool verbose) {
  auto rmsdCutPair = [](const vector<int>& fragResIdx, const Structure& S) { return RMSDCalculator::rmsdCutoff(fragResIdx, S, 1.0, 20); };
  if ((Ri->getStructure() == NULL) || (Rj->getStructure() == NULL)) MstUtils::error("cannot operate on a disembodied residues!", "dTERMen::pairEnergies(Residue*, Residue*, bool)");
  if (Ri->getStructure() != Rj->getStructure()) MstUtils::error("specified residues belong to different structures!", "dTERMen::pairEnergies(Residue*, Residue*, bool)");
  Structure& S = *(Ri->getStructure());

  // isolate TERM and get matches
  int naa = globalAlphabetSize();
  termData pT({Ri, Rj}, pmPair);
  F.setOptions(foptsBase);
  F.setQuery(pT.getTERM());
  F.setRMSDCutoff(rmsdCutPair(pT.getResidueIndices(), S));
  F.setMinNumMatches(pairMinN);
  F.setMaxNumMatches(pairMaxN);
  pT.setMatches(F.search(), homCut, &F);
  if (recordData) data.push_back(pT);

  // for each of the two positions, compute expectation of every amino acid in
  // the context of every match, based on background "trivial" energies and
  // normalizing the totals of each amino acid across all matches to be equal to
  // the number observed
  vector<CartesianPoint> Pi, Pj;
  vector<vector<mstreal> > pairE(naa, vector<mstreal>(naa, 0.0));
  int cIndI = pT.getCentralResidueIndices()[0]; int cIndJ = pT.getCentralResidueIndices()[1];
  CartesianPoint NoI = dTERMen::singleBodyObservations(pT.getMatches(), cIndI);
  CartesianPoint NeI = dTERMen::singleBodyExpectations(pT.getMatches(), cIndI, &Pi);
  CartesianPoint NoJ = dTERMen::singleBodyObservations(pT.getMatches(), cIndJ);
  CartesianPoint NeJ = dTERMen::singleBodyExpectations(pT.getMatches(), cIndJ, &Pj);
  CartesianPoint NoIJ = dTERMen::twoBodyObservations(pT.getMatches(), cIndI, cIndJ);
  for (int aai = 0; aai < naa; aai++) {
    for (int aaj = 0; aaj < naa; aaj++) {
      mstreal Ne = 0;
      for (int k = 0; k < pT.numMatches(); k++) {
        if (Pi[k].empty() || Pj[k].empty()) continue; // e.g., relevant positions in this match had unknown identities
        Ne += Pi[k][aai] * Pj[k][aaj];
      }
      Ne *= (NoI[aai]/NeI[aai])*(NoJ[aaj]/NeJ[aaj]); // normalize to preserve marginals
      mstreal No = NoIJ[dTERMen::pairToIdx(aai, aaj)];
      mstreal pc = naa/MstUtils::max((vector<mstreal>) {No, Ne, 1.0});
      pairE[aai][aaj] = -kT*log((No + pc)/(Ne + pc));
    }
  }

  return pairE;
}

vector<vector<mstreal>> dTERMen::pairEnergiesNew(Residue* Ri, Residue* Rj, bool verbose) {
  auto rmsdCutPair = [](const vector<int>& fragResIdx, const Structure& S) { return RMSDCalculator::rmsdCutoff(fragResIdx, S, 1.0, 20); };
  if ((Ri->getStructure() == NULL) || (Rj->getStructure() == NULL)) MstUtils::error("cannot operate on a disembodied residues!", "dTERMen::pairEnergies(Residue*, Residue*, bool)");
  if (Ri->getStructure() != Rj->getStructure()) MstUtils::error("specified residues belong to different structures!", "dTERMen::pairEnergies(Residue*, Residue*, bool)");
  Structure& S = *(Ri->getStructure());

  // isolate TERM and get matches
  int naa = globalAlphabetSize();
  termData pT({Ri, Rj}, pmPair);
  F.setOptions(foptsBase);
  F.setQuery(pT.getTERM(), false);
  F.setRMSDCutoff(rmsdCutPair(pT.getResidueIndices(), S));
  F.setMinNumMatches(pairMinN);
  F.setMaxNumMatches(pairMaxN);
  pT.setMatches(F.search(), homCut, &F);
  if (recordData) data.push_back(pT);

  // figure out which matches are from "homo-dimers"
  vector<bool> isHomoDimer(pT.getMatches().size(), false);
  for (int i = 0; i < pT.getMatches().size(); i++) {
    // in order to count as a homo-dimer match, the central residues of the two
    // segments must be the same, their sequence context must match closely
    const fasstSolution& match = pT.getMatches()[i];
    Sequence seq = F.getMatchSequence(match);
    if (seq[pT.getCentralResidueIndices()[0]] != seq[pT.getCentralResidueIndices()[1]]) continue;
    fasstSolution invertedMatch(match.getAlignment(), match.getRMSD(), match.getTargetIndex(), match.getTransform(), match.getSegLengths(), {1, 0});
    fasstSolutionSet setOrig(match);
    fasstSolutionSet setInverted(invertedMatch);
    getFASST()->addSequenceContext(setOrig);
    getFASST()->addSequenceContext(setInverted);
    isHomoDimer[i] = !(setOrig.insert(setInverted[0], 0.9));
  }

  // for each of the two positions, compute expectation of every amino acid in
  // the context of every match, based on background "trivial" energies and
  // normalizing the totals of each amino acid across all matches to be equal to
  // the number observed
  vector<vector<mstreal> > pairE(naa, vector<mstreal>(naa, 0.0));
  int cIndI = pT.getCentralResidueIndices()[0]; int cIndJ = pT.getCentralResidueIndices()[1];
  vector<CartesianPoint> Pi = dTERMen::singleBodyExpectationsMatchedMarginals(pT.getMatches(), cIndI);
  vector<CartesianPoint> Pj = dTERMen::singleBodyExpectationsMatchedMarginals(pT.getMatches(), cIndJ);
  CartesianPoint NoIJ = dTERMen::twoBodyObservations(pT.getMatches(), cIndI, cIndJ);
  for (int aai = 0; aai < naa; aai++) {
    for (int aaj = 0; aaj < naa; aaj++) {
      mstreal Ne = 0;
      for (int k = 0; k < pT.numMatches(); k++) {
        if (Pi[k].empty() || Pj[k].empty()) continue; // e.g., relevant positions in this match had unknown identities
        if (isHomoDimer[k]) Ne += (Pi[k][aai] + Pj[k][aaj])/2;
        else Ne += Pi[k][aai] * Pj[k][aaj];
      }
      mstreal No = NoIJ[dTERMen::pairToIdx(aai, aaj)];
      mstreal pc = naa/MstUtils::max((vector<mstreal>) {No, Ne, 1.0});
      pairE[aai][aaj] = -kT*log((No + pc)/(Ne + pc));
    }
  }

  return pairE;
}

vector<vector<mstreal>> dTERMen::pairEnergiesNew2(Residue* Ri, Residue* Rj, bool verbose) {
  auto rmsdCutPair = [](const vector<int>& fragResIdx, const Structure& S) { return RMSDCalculator::rmsdCutoff(fragResIdx, S, 1.0, 20); };
  if ((Ri->getStructure() == NULL) || (Rj->getStructure() == NULL)) MstUtils::error("cannot operate on a disembodied residues!", "dTERMen::pairEnergies(Residue*, Residue*, bool)");
  if (Ri->getStructure() != Rj->getStructure()) MstUtils::error("specified residues belong to different structures!", "dTERMen::pairEnergies(Residue*, Residue*, bool)");
  Structure& S = *(Ri->getStructure());

  // isolate TERM and get matches
  int naa = globalAlphabetSize();
  termData pT({Ri, Rj}, pmPair);
  F.setOptions(foptsBase);
  F.setQuery(pT.getTERM());
  F.setRMSDCutoff(rmsdCutPair(pT.getResidueIndices(), S));
  F.setMinNumMatches(pairMinN);
  F.setMaxNumMatches(pairMaxN);
  pT.setMatches(F.search(), homCut, &F);
  if (recordData) data.push_back(pT);

  // figure out which matches are from "homo-dimers"
  vector<bool> isHomoDimer(pT.getMatches().size(), false);
  for (int i = 0; i < pT.getMatches().size(); i++) {
    // in order to count as a homo-dimer match, the central residues of the two
    // segments must be the same, their sequence context must match closely
    const fasstSolution& match = pT.getMatches()[i];
    Sequence seq = F.getMatchSequence(match);
    if (seq[pT.getCentralResidueIndices()[0]] != seq[pT.getCentralResidueIndices()[1]]) continue;
    fasstSolution invertedMatch(match.getAlignment(), match.getRMSD(), match.getTargetIndex(), match.getTransform(), match.getSegLengths(), {1, 0});
    fasstSolutionSet setOrig(match);
    fasstSolutionSet setInverted(invertedMatch);
    getFASST()->addSequenceContext(setOrig);
    getFASST()->addSequenceContext(setInverted);
    isHomoDimer[i] = !(setOrig.insert(setInverted[0], 0.9));
  }

  // for each of the two positions, compute expectation of every amino acid in
  // the context of every match, based on background "trivial" energies and
  // normalizing the totals of each amino acid across all matches to be equal to
  // the number observed
  vector<CartesianPoint> Pi, Pj;
  vector<vector<mstreal> > pairE(naa, vector<mstreal>(naa, 0.0));
  int cIndI = pT.getCentralResidueIndices()[0]; int cIndJ = pT.getCentralResidueIndices()[1];
  CartesianPoint NoI = dTERMen::singleBodyObservations(pT.getMatches(), cIndI);
  CartesianPoint NeI = dTERMen::singleBodyExpectations(pT.getMatches(), cIndI, &Pi);
  CartesianPoint NoJ = dTERMen::singleBodyObservations(pT.getMatches(), cIndJ);
  CartesianPoint NeJ = dTERMen::singleBodyExpectations(pT.getMatches(), cIndJ, &Pj);
  CartesianPoint NoIJ = dTERMen::twoBodyObservations(pT.getMatches(), cIndI, cIndJ);
  for (int aai = 0; aai < naa; aai++) {
    for (int aaj = 0; aaj < naa; aaj++) {
      mstreal Ne = 0;
      for (int k = 0; k < pT.numMatches(); k++) {
        if (Pi[k].empty() || Pj[k].empty()) continue; // e.g., relevant positions in this match had unknown identities
        if (isHomoDimer[k]) Ne += ((Pi[k][aai] * (NoI[aai]/NeI[aai])) + (Pj[k][aaj] * (NoJ[aaj]/NeJ[aaj])))/2;
        else Ne += (Pi[k][aai] * (NoI[aai]/NeI[aai])) * (Pj[k][aaj] * (NoJ[aaj]/NeJ[aaj]));
      }
      mstreal No = NoIJ[dTERMen::pairToIdx(aai, aaj)];
      mstreal pc = naa/MstUtils::max((vector<mstreal>) {No, Ne, 1.0});
      pairE[aai][aaj] = -kT*log((No + pc)/(Ne + pc));
    }
  }

  return pairE;
}

void dTERMen::printSelfComponent(const CartesianPoint& ener, const string& prefix) {
  cout << prefix;
  for (int i = 0; i < globAlph.size(); i++) printf("%8s", indexToResName(i).c_str());
  cout << endl << prefix;
  for (int i = 0; i < ener.size(); i++) printf("%8.3f", ener[i]);
  cout << endl;
}

CartesianPoint dTERMen::singleBodyStatEnergy(fasstSolutionSet& matches, int cInd, int pc) {
  CartesianPoint selfE(globalAlphabetSize(), 0.0);
  CartesianPoint Ne = dTERMen::singleBodyExpectations(matches, cInd);
  CartesianPoint No = dTERMen::singleBodyObservations(matches, cInd);
  for (int aai = 0; aai < selfE.size(); aai++) {
    selfE[aai] = -kT*log((No[aai] + pc)/(Ne[aai] + pc));
  }

  return selfE;
}

CartesianPoint dTERMen::singleBodyObservations(fasstSolutionSet& matches, int cInd) {
  CartesianPoint No(globalAlphabetSize(), 0.0);
  for (int i = 0; i < matches.size(); i++) {
    int aaIdx = dTERMen::aaToIndex(F.getMatchSequence(matches[i])[cInd]);
    if (!aaIndexKnown(aaIdx)) continue; // this match has some amino acid that is not known in the current alphabet
    No[aaIdx] += 1.0;
  }
  return No;
}

CartesianPoint dTERMen::twoBodyObservations(fasstSolutionSet& matches, int cIndI, int cIndJ) {
  CartesianPoint No(globalAlphabetSize()*globalAlphabetSize(), 0.0);
  for (int i = 0; i < matches.size(); i++) {
    int aaiIdx = dTERMen::aaToIndex(F.getMatchSequence(matches[i])[cIndI]);
    int aajIdx = dTERMen::aaToIndex(F.getMatchSequence(matches[i])[cIndJ]);
    if (!aaIndexKnown(aaiIdx) || !aaIndexKnown(aajIdx)) continue;
    No[dTERMen::pairToIdx(aaiIdx, aajIdx)] += 1.0;
  }
  return No;
}

CartesianPoint dTERMen::singleBodyExpectations(fasstSolutionSet& matches, int cInd, vector<CartesianPoint>* breakDown) {
  CartesianPoint Ne(globalAlphabetSize(), 0.0);
  if (breakDown != NULL) { breakDown->clear(); breakDown->resize(matches.size()); }
  for (int i = 0; i < matches.size(); i++) {
    int aaIdx = dTERMen::aaToIndex(F.getMatchSequence(matches[i])[cInd]);
    if (!aaIndexKnown(aaIdx)) continue; // this match has some amino acid that is not known in the current alphabet
    CartesianPoint inMatchExp = backExpectation(matches[i], cInd);
    if (breakDown != NULL) (*breakDown)[i] = inMatchExp;
    Ne += inMatchExp;
  }
  return Ne;
}

vector<CartesianPoint> dTERMen::singleBodyExpectationsMatchedMarginals(fasstSolutionSet& matches, int cInd) {
  // compute initial expected amino-acid distributions in each match
  vector<CartesianPoint> mlogP;
  CartesianPoint Ne(globalAlphabetSize(), 0.0);
  vector<int> validMatches;
  for (int i = 0; i < matches.size(); i++) {
    int aaIdx = dTERMen::aaToIndex(F.getMatchSequence(matches[i])[cInd]);
    if (!aaIndexKnown(aaIdx)) continue; // this match has some amino acid that is not known in the current alphabet
    validMatches.push_back(i);
    CartesianPoint inMatchExp = backExpectation(matches[i], cInd);
    Ne += inMatchExp;
    mlogP.push_back(inMatchExp);
    for (int j = 0; j < mlogP.back().size(); j++) mlogP.back()[j] = -log(mlogP.back()[j]);
  }
  CartesianPoint No = dTERMen::singleBodyObservations(matches, cInd);

  // next, iterate to find optimal amino-acid bias energies to get marginals to agree
  CartesianPoint delE(globalAlphabetSize(), 0.0);
  vector<CartesianPoint> breakDown(matches.size());
  int Ncyc = 10;
  for (int c = 0; c < Ncyc; c++) {
    for (int aai = 0; aai < delE.size(); aai++) delE[aai] = delE[aai] + log(Ne[aai]/No[aai]);
    Ne.clear(); Ne.resize(globalAlphabetSize(), 0.0);
    for (int i = 0; i < mlogP.size(); i++) {
      CartesianPoint inMatchExp = mlogP[i] + delE;
      dTERMen::enerToProb(inMatchExp);
      Ne += inMatchExp;
      if (c == Ncyc - 1) {
        breakDown[validMatches[i]] = inMatchExp;
      }
    }
  }

  return breakDown;
}

CartesianPoint dTERMen::backExpectation(const fasstSolution& m, int cInd) {
  vector<mstreal> p(globalAlphabetSize(), 0);
  mstreal phi = F.getResidueProperties(m, "phi")[cInd];
  mstreal psi = F.getResidueProperties(m, "psi")[cInd];
  mstreal omg = F.getResidueProperties(m, "omega")[cInd];
  mstreal env = F.getResidueProperties(m, "env")[cInd];
  for (int aai = 0; aai < globalAlphabetSize(); aai++) {
    p[aai] = backEner(aai) + bbOmegaEner(omg, aai) + bbPhiPsiEner(phi, psi, aai) + envEner(env, aai);
  }
  dTERMen::enerToProb(p);
  return p;
}

mstreal dTERMen::enerToProb(vector<mstreal>& ener) {
  mstreal minEner = MstUtils::min(ener);
  mstreal Z = 0;
  for (int i = 0; i < ener.size(); i++) { ener[i] = exp(-(ener[i] - minEner)/kT); Z += ener[i]; }
  for (int i = 0; i < ener.size(); i++) ener[i] /= Z;
  return Z;
}

int dTERMen::pairToIdx(int aai, int aaj) const {
  return aai*globalAlphabetSize() + aaj;
}

pair<int, int> dTERMen::idxToPair(int idx) const {
  return pair<int, int>(idx / globalAlphabetSize(), idx % globalAlphabetSize());
}

/* ----------- EnergyTable --------------- */
int termData::setMatches(const fasstSolutionSet& _matches, mstreal homologyCutoff, FASST* F) {
  if (F == NULL) { matches = _matches; return matches.size(); }

  // need at least one match (to structure the dummy match around) and at least
  // one central residue for the underlying TERM (to infer the parent Structure)
  if (_matches.size() == 0) return 0;
  if (centResidues.empty()) return 0;
  Structure* querySource = centResidues[0]->getStructure();
  if (querySource == NULL) return 0;

  // Create a dummy solution corresponding to the query, and fill its sequence
  // context based on the template it came from. Use an example match to figure
  // out segment boundaries, so that this decision is made only once (in FASST)
  const fasstSolution& ex = *(_matches.begin());
  vector<int> algn(ex.numSegments()), segLen(ex.numSegments());
  int off = 0;
  for (int i = 0; i < ex.numSegments(); i++) {
    algn[i] = fragResIdx[off];
    segLen[i] = ex.segLength(i);
    off += ex.segLength(i);
  }
  fasstSolution orig(fasstSolutionAddress(0, algn), segLen);
  orig.addSequenceContext(*querySource, F->options().getContextLength());

  // fill the sequence context of all solutions
  fasstSolutionSet matchesWithContext = _matches;
  F->addSequenceContext(matchesWithContext);

  // create a solution set with just the dummy solution and add each of the real
  // solutions into it, one by one, with redundancy turned on. If they are not
  // inserted, then they are redundant with the query and should be removed
  orig.setRMSD(-1.0); // make sure the dummy solution always stays on top!
  fasstSolutionSet origSet(orig);
  matches.clear();
  for (int i = 0; i < matchesWithContext.size(); i++) {
    if (origSet.insert(matchesWithContext[i], homologyCutoff)) {
      matches.insert(matchesWithContext[i]);
      origSet.erase(matchesWithContext[i]);
      if (origSet.size() != 1) MstUtils::error("it appears the original dummy solution has been lost (this should not happen)", "dTERMen::termData::setMatches"); // sanity check
    }
  }

  return matches.size();
}

/* ----------- EnergyTable --------------- */
EnergyTable::EnergyTable(const string& tabFile) {
  readFromFile(tabFile);
}

EnergyTable EnergyTable::restrictSiteAlphabet(const vector<vector<string>>& restricted_siteAlphabets, bool constraint_table) {
  EnergyTable restricted_etab;

  if (restricted_siteAlphabets.size() != sites.size()) MstUtils::error("The length of restricted_siteAlphabets and the number of positions in the energy table do not match ("+MstUtils::toString(restricted_siteAlphabets.size())+") and ("+MstUtils::toString(sites.size())+")");

  /* copy site addresses and set alphabets */
  for (int s = 0; s < sites.size(); s++) {
    if ((constraint_table) && (restricted_siteAlphabets[s].size() > 1)) MstUtils::error("In constraint_table mode only a single residue, or none at all, may be specified per position");
    restricted_etab.addSite(sites[s]);
    if (restricted_siteAlphabets[s].empty() || ((constraint_table && restricted_siteAlphabets[s].size() == 1))) {
      // If no residue is specified, copy all residue types
      vector<string> old_siteAlphabet = getSiteAlphabet(s);
      restricted_etab.setSiteAlphabet(s,old_siteAlphabet);
    } else {
      restricted_etab.setSiteAlphabet(s,restricted_siteAlphabets[s]);
    }
  }
  /* copy self energies */
  for (int s = 0; s < sites.size(); s++) {
    vector<string> siteAlphabet = restricted_etab.getSiteAlphabet(s);
    for (int site_aa = 0; site_aa < siteAlphabet.size(); site_aa++) {
      string aa = siteAlphabet[site_aa];
      if (!inSiteAlphabet(s,aa)) MstUtils::error("unexpected amino-acid name in restricted alphabet: " + aa);
      /*
       In the case that this is a constraint table AND this position is specified, the energy of the
       sequence provided in restricted_siteAlphabets[s][0] should be copied to all of the other
       residue types.
       */
      int old_aa_idx = ((constraint_table) && (restricted_siteAlphabets[s].size() == 1)) ? indexInSiteAlphabet(s,restricted_siteAlphabets[s][0]) : indexInSiteAlphabet(s,aa);
      int new_aa_idx = restricted_etab.indexInSiteAlphabet(s,aa);
      mstreal self_e = selfEnergy(s,old_aa_idx);
      restricted_etab.setSelfEnergy(s, new_aa_idx, self_e);
    }
  }
  /* copy pair energies */
  for (int si = 0; si < sites.size(); si++) {
    for (auto it = pairMaps[si].begin(); it != pairMaps[si].end(); ++it) {
      int sj = it->first;
      if (sj < si) continue;
      int k = it->second;
      vector<string> si_siteAlphabet = restricted_etab.getSiteAlphabet(si);
      vector<string> sj_siteAlphabet = restricted_etab.getSiteAlphabet(sj);
      for (int aai = 0; aai < si_siteAlphabet.size(); aai++) {
        for (int aaj = 0; aaj < sj_siteAlphabet.size(); aaj++) {
          /*
           If this is a constraint table:
           In the case that both positions are specified, the pair energy of the residue types
           specified by restricted_siteAlphabets[s][0] is copied to all pairs of residue types.

           In the case that one of the positions is specified and the other is not, there will be 20
           unique pair energies (1 specified residue x 20 residues). Each of these pair energies will
           be copied 19 times, for each of the other residue types at the specified position.
           */
          int aai_old = ((constraint_table) && (restricted_siteAlphabets[si].size() == 1)) ? indexInSiteAlphabet(si,restricted_siteAlphabets[si][0]) : indexInSiteAlphabet(si,si_siteAlphabet[aai]);
          int aaj_old = ((constraint_table) && (restricted_siteAlphabets[sj].size() == 1)) ? indexInSiteAlphabet(sj,restricted_siteAlphabets[sj][0]) : indexInSiteAlphabet(sj,sj_siteAlphabet[aaj]);
          mstreal ener = pairE[si][k][aai_old][aaj_old];
          restricted_etab.setPairEnergy(si, sj, aai, aaj, ener);
        }
      }
    }
  }
  return restricted_etab;
};

EnergyTable EnergyTable::restrictSiteAlphabet(const Structure& S, bool constraint_table) {
  //Use the sequence of the provided structure to restrict the site alphabet at each position of the energy table
  //  Positions with a residue type are restricted to that type
  vector<Residue*> S_seq = S.getResidues();
  if (S_seq.size() != sites.size()) MstUtils::error("The protein size must be the same as the provided energy table");
  vector<vector<string>> all_siteNames; all_siteNames.resize(S_seq.size());
  for (int R_idx = 0; R_idx < S_seq.size(); R_idx++) {
    string R_name = S_seq[R_idx]->getName();
    all_siteNames[R_idx].push_back(R_name);
  }
  return restrictSiteAlphabet(all_siteNames, constraint_table);
};

void EnergyTable::addSite(const string& siteName) {
  if (siteIndices.find(siteName) != siteIndices.end()) MstUtils::error("site '" + siteName + "' is already present!", "EnergyTable::addSite(const string&)");
  siteIndices[siteName] = sites.size();
  sites.push_back(siteName);
  aaAlpha.resize(aaAlpha.size() + 1);
  aaIndices.resize(aaIndices.size() + 1);
  selfE.resize(selfE.size() + 1);
  pairE.resize(pairE.size() + 1);
  pairMaps.resize(pairMaps.size() + 1);
}

void EnergyTable::addSites(const vector<string>& siteNames) {
  for (int i = 0; i < siteNames.size(); i++) addSite(siteNames[i]);
}

void EnergyTable::setSiteAlphabet(int siteIdx, const vector<string>& alpha) {
//  if (!empty()) MstUtils::error("site alphabets must be set before populating energies", "EnergyTable::setSiteAlphabet(int, const vector<string>&)");
  if (aaAlpha.size() < siteIdx + 1) MstUtils::error("site index out of range", "EnergyTable::setSiteAlphabet(int, const vector<string>&)");
  aaAlpha[siteIdx] = alpha;
  for (int i = 0; i < alpha.size(); i++) aaIndices[siteIdx][alpha[i]] = i;
  selfE[siteIdx].clear(); selfE[siteIdx].resize(alpha.size(), 0.0);
  for (auto it = pairMaps[siteIdx].begin(); it != pairMaps[siteIdx].end(); ++it) {
    int sj = it->first;
    int k = it->second;
    pairE[siteIdx][k].clear();
    pairE[siteIdx][k].resize(aaIndices[siteIdx].size(), vector<mstreal>(aaIndices[sj].size()));
  }
}

int EnergyTable::addToSiteAlphabet(int siteIdx, const string& aa) {
  if (aaIndices[siteIdx].find(aa) != aaIndices[siteIdx].end()) MstUtils::error("tried to add an amino acid that already exists at the site!", "EnergyTable::addToSiteAlphabet(int, const string&)");
  int a = aaIndices[siteIdx].size();
  aaIndices[siteIdx][aa] = a;
  aaAlpha[siteIdx].push_back(aa);
  selfE[siteIdx].push_back(0.0);
  return a;
}

void EnergyTable::renameSiteResidue(int si, int ai, const string& aa) {
  string oldAA = aaAlpha[si][ai];
  aaAlpha[si][ai] = aa;
  aaIndices[si].erase(oldAA);
  aaIndices[si][aa] = ai;
}

void EnergyTable::clear() {
  siteIndices.clear();
  sites.clear();
  aaIndices.clear();
  aaAlpha.clear();
  selfE.clear();
  pairMaps.clear();
  pairE.clear();
}

void EnergyTable::readFromFile(const string& tabFile) {
  clear();
  vector<string> lines = MstUtils::fileToArray(tabFile);
  int i = 0, a;

  /* read self energies and site addresses */
  for (; i < lines.size(); i++) {
    vector<string> line = MstUtils::split(lines[i]);
    if (line.size() != 3) break;
    if (!siteExists(line[0])) addSite(line[0]);
    int siteIdx = siteIndex(line[0]);
    int aaIdx = addToSiteAlphabet(siteIdx, line[1]);
    setSelfEnergy(siteIdx, aaIdx, MstUtils::toReal(line[2]));
  }

  /* read pair energies */
  for (; i < lines.size(); i++) {
    vector<string> line = MstUtils::split(lines[i]);
    if (line.size() != 5) MstUtils::error("could not parse line '" + lines[i] + "'");
    if (!siteExists(line[0]) || !siteExists(line[1])) MstUtils::error("unexpected site names in pair line: " + lines[i]);
    int si = siteIndex(line[0]);
    int sj = siteIndex(line[1]);
    if (!inSiteAlphabet(si, line[2]) || !inSiteAlphabet(si, line[3])) MstUtils::error("unexpected amino-acid names in pair line: " + lines[i]);
    int aai = indexInSiteAlphabet(si, line[2]);
    int aaj = indexInSiteAlphabet(si, line[3]);
    setPairEnergy(si, sj, aai, aaj, MstUtils::toReal(line[4]));
  }
}

void EnergyTable::writeToFile(const string& tabFile) {
  fstream of;
  MstUtils::openFile(of, tabFile, ios::out);
  for (int si = 0; si < sites.size(); si++) {
    //if there is only one possible amino acid, this position cannot be optimized and should not be written
    if (aaAlpha[si].size() == 1) continue;
    for (int aai = 0; aai < aaAlpha[si].size(); aai++) {
      of << sites[si] << " " << aaAlpha[si][aai] << " " << selfE[si][aai] << endl;
    }
  }
  for (int si = 0; si < sites.size(); si++) {
    for (auto it = pairMaps[si].begin(); it != pairMaps[si].end(); ++it) {
      int sj = it->first;
      if (sj < si) continue;
      int k = it->second;
      if ((aaAlpha[si].size() == 1) && (aaAlpha[sj].size() == 1)) continue; //don't write if there is one possible pair combination
      for (int aai = 0; aai < aaAlpha[si].size(); aai++) {
        for (int aaj = 0; aaj < aaAlpha[sj].size(); aaj++) {
          mstreal ener = pairE[si][k][aai][aaj];
          if (ener == 0.0) continue;
          of << sites[si] << " " << sites[sj] << " " << aaAlpha[si][aai] << " " << aaAlpha[sj][aaj] << " " << ener << endl;
        }
      }
    }
  }
  of.close();
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
  if (sj < si) {
    int tmp = si; si = sj; sj = tmp;
    tmp = aai; aai = aaj; aaj = tmp;
  }
  if (pairMaps[si].find(sj) == pairMaps[si].end()) return 0.0;
  return pairE[si][pairMaps[si][sj]][aai][aaj];
}

void EnergyTable::setSelfEnergy(int s, int aa, mstreal ener) {
  selfE[s][aa] = ener;
}

void EnergyTable::setPairEnergy(int si, int sj, int aai, int aaj, mstreal ener) {
  // store each interaction energy in one order only
  if (sj < si) {
    int tmp = si; si = sj; sj = tmp;
    tmp = aai; aai = aaj; aaj = tmp;
  }
  if (pairMaps[si].find(sj) == pairMaps[si].end()) {
    pairMaps[si][sj] = pairE[si].size();
    // we store index pairs in both directions, to make it easy to look up all
    // interactors of a given site, but we access in only one direction
    pairMaps[sj][si] = pairE[si].size();
    pairE[si].resize(pairE[si].size() + 1);
  }
  int k = pairMaps[si][sj];

  // the first time we encounter a site pair, fill them up with all-by-all amino-acid zero energies
  if (pairE[si][k].size() == 0) {
    pairE[si][k].resize(aaIndices[si].size(), vector<mstreal>(aaIndices[sj].size(), 0));
  }

  // finally set the energy
  pairE[si][k][aai][aaj] = ener;
}

mstreal EnergyTable::scoreSolution(const vector<int>& sol) {
  if (sol.size() != selfE.size()) MstUtils::error("solution of wrong length for table", "EnergyTable::scoreSolution(const vector<int>&)");
  mstreal ener = 0;
  for (int si = 0; si < selfE.size(); si++) {
    ener += selfE[si][sol[si]];
    for (auto it = pairMaps[si].begin(); it != pairMaps[si].end(); ++it) {
      if (it->first < si) continue; // do not overcount pairs
      ener += pairE[si][it->second][sol[si]][sol[it->first]];
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
  if (!pairMaps[mutSite].empty()) {
    for (auto it = pairMaps[mutSite].begin(); it != pairMaps[mutSite].end(); ++it) {
      int intSite = it->first; // interacting site index
      int k = it->second;      // index at which this site's interactions will be stored
      if (intSite > mutSite) {
        dE += pairE[mutSite][k][mutAA][sol[intSite]] - pairE[mutSite][k][sol[mutSite]][sol[intSite]];
      } else {
        dE += pairE[intSite][k][sol[intSite]][mutAA] - pairE[intSite][k][sol[intSite]][sol[mutSite]];
      }
    }
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

vector<int> EnergyTable::mc(int Nc, int Ni, mstreal kTi, mstreal kTf, int annealType, void* rec, void (*add)(void*, const vector<int>&, mstreal), int Ne, void* extra, mstreal (*additionalScore)(void*, const vector<int>&, EnergyTable&, int mutSite, int mutAA)) {
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
      if (additionalScore != NULL) dE = (*additionalScore)(extra, seq, *this, s, aa);
      if (MstUtils::randUnit() < exp(-dE/kT)) seq[s] = aa;
    }
    // data collection
    mstreal bcE = scoreSolution(seq); vector<int> bcS = seq;
    if (additionalScore != NULL) bcE += (*additionalScore)(extra, seq, *this, -1, -1);
    if (cyc == 0) { bE = bcE; bS = bcS; }
    mstreal dEb = 0; // current change from the best energy
    for (int i = 0; i < Ni; i++) {
      mstreal f = i*1.0/(Ni-1);
      // every now and again, calculate the total energy to avoid accumulation of addition errors
      if ((i+1) % 10000 == 0) {
        dEb = scoreSolution(seq) - bcE;
        if (additionalScore != NULL) dEb += (*additionalScore)(extra, seq, *this, -1, -1);
      }
      int s = MstUtils::randInt(0, L - 1);
      int aa = randomResidue(s);
      mstreal dE = scoreMutation(seq, s, aa);
      if (additionalScore != NULL) dE += (*additionalScore)(extra, seq, *this, s, aa);
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

vector<int> EnergyTable::sequenceToSolution(const Sequence& seq, bool strict) {
  if (seq.size() != selfE.size()) MstUtils::error("sequence of wrong length for table", "EnergyTable::sequenceToSolution(const Sequence&)");
  vector<int> seqInts(seq.size());
  for (int i = 0; i < seq.size(); i++) {
    string aa3 = seq.getResidue(i, true);
    if (aaIndices[i].find(aa3) == aaIndices[i].end()) {
      if (strict) MstUtils::error("sequence is not from table alphabet", "EnergyTable::sequenceToSolution(const Sequence&)");
      seqInts[i] = -1;
    } else {
      seqInts[i] = aaIndices[i][aa3];
    }
  }
  return seqInts;
}

string EnergyTable::getResidueString(int si, int ri) {
  return aaAlpha[si][ri];
}
