#include "mstsequence.h"
using namespace MST;

/* ------------ SeqTools ------------ */

// statics must be defined in the cpp file
vector<string> SeqTools::aa1, SeqTools::aa3;
map<string, res_t> SeqTools::aa3ToIdx, SeqTools::aa1ToIdx;
vector<string> SeqTools::idxToAA3, SeqTools::idxToAA1;
res_t SeqTools::_unkIdx, SeqTools::_gapIdx;
bool initialized = SeqTools::initConstants();

res_t SeqTools::aaToIdx(const string& aa) {
  if (aa.length() == 1) {
    if (aa1ToIdx.find(aa) == aa1ToIdx.end()) return SeqTools::unknownIdx();
    return aa1ToIdx[aa];
  } else if (aa.length() == 3) {
    if (aa3ToIdx.find(aa) == aa3ToIdx.end()) return SeqTools::unknownIdx();
    return aa3ToIdx[aa];
  }
  MstUtils::error("uknown amino acid '" + aa + "'", "SeqTools::aaToIdx");
  return 0; // for the compiler to be happy
}

string SeqTools::idxToTriple(res_t idx) {
  MstUtils::assertCond((idx >= 0) && (idx < idxToAA3.size()), "unknown amino-acid index '" + MstUtils::toString(idx) + "'", "SeqTools::idxToTriple");
  return idxToAA3[idx];
}

string SeqTools::idxToSingle(res_t idx) {
  MstUtils::assertCond((idx >= 0) && (idx < idxToAA1.size()), "unknown amino-acid index '" + MstUtils::toString(idx) + "'", "SeqTools::idxToSingle");
  return idxToAA1[idx];
}

string SeqTools::toTriple(const string& aa) {
  return idxToTriple(aaToIdx(aa));
}

string SeqTools::toSingle(const string& aa) {
  return SeqTools::idxToSingle(SeqTools::aaToIdx(aa));
}

string SeqTools::tripleToSingle(const string& triple, const string& del) {
  vector<string> seq = MstUtils::split(triple, del);
  string single;
  for (int i = 0; i < seq.size(); i++) {
    single += SeqTools::toSingle(seq[i]);
  }
  return single;
}

string SeqTools::singleToTriple(const string& single, const string& del) {
  vector<string> seq = MstUtils::split(single, del);
  string triple;
  for (int i = 0; i < seq.size(); i++) {
    triple += SeqTools::toTriple(seq[i]);
    if (i != seq.size() - 1) triple += " ";
  }
  return triple;
}

vector<res_t> SeqTools::seqToIdx(const string& seqstr, const string& del) {
  vector<string> seq = MstUtils::split(seqstr, del);
  vector<res_t> indices(seq.size());
  for (int i = 0; i < seq.size(); i++) {
    indices[i] = SeqTools::aaToIdx(seq[i]);
  }
  return indices;
}

bool SeqTools::initConstants() {
  aa3.push_back("ALA"); aa1.push_back("A");
  aa3.push_back("CYS"); aa1.push_back("C");
  aa3.push_back("ASP"); aa1.push_back("D");
  aa3.push_back("GLU"); aa1.push_back("E");
  aa3.push_back("PHE"); aa1.push_back("F");
  aa3.push_back("GLY"); aa1.push_back("G");
  aa3.push_back("HIS"); aa1.push_back("H");
  aa3.push_back("ILE"); aa1.push_back("I");
  aa3.push_back("LYS"); aa1.push_back("K");
  aa3.push_back("LEU"); aa1.push_back("L");
  aa3.push_back("MET"); aa1.push_back("M");
  aa3.push_back("ASN"); aa1.push_back("N");
  aa3.push_back("PRO"); aa1.push_back("P");
  aa3.push_back("GLN"); aa1.push_back("Q");
  aa3.push_back("ARG"); aa1.push_back("R");
  aa3.push_back("SER"); aa1.push_back("S");
  aa3.push_back("THR"); aa1.push_back("T");
  aa3.push_back("VAL"); aa1.push_back("V");
  aa3.push_back("TRP"); aa1.push_back("W");
  aa3.push_back("TYR"); aa1.push_back("Y");

  aa3.push_back("HSD"); aa1.push_back("H");
  aa3.push_back("HSE"); aa1.push_back("H");
  aa3.push_back("HSC"); aa1.push_back("H");
  aa3.push_back("HSP"); aa1.push_back("H");
  aa3.push_back("MSE"); aa1.push_back("M");

  aa3.push_back("CSO"); aa1.push_back("X"); // S-hydroxycysteine
  aa3.push_back("HIP"); aa1.push_back("H"); // ND1-phosphohistidine
  aa3.push_back("SEC"); aa1.push_back("C"); // selenocysteine
  aa3.push_back("SEP"); aa1.push_back("S"); // phosphoserine
  aa3.push_back("TPO"); aa1.push_back("T"); // phosphothreonine
  aa3.push_back("PTR"); aa1.push_back("Y"); // o-phosphotyrosine
  aa3.push_back("UNK"); aa1.push_back("?"); // unknown residue
  aa3.push_back("---"); aa1.push_back("-"); // gap

  idxToAA3.resize(aa3.size());
  idxToAA1.resize(aa1.size());
  for (int i = 0; i < aa3.size(); i++) {
    aa3ToIdx[aa3[i]] = (res_t) i;
    // single-letter code is ambiguous, take the first amino acid
    // when going from single-letter code to index
    if (aa1ToIdx.find(aa1[i]) == aa1ToIdx.end()) aa1ToIdx[aa1[i]] = i;
    idxToAA3[i] = aa3[i];
    idxToAA1[i] = aa1[i];
  }
  _unkIdx = aa1ToIdx["?"];
  _gapIdx = aa1ToIdx["-"];
  return true;
}

void SeqTools::readFasta(const string& fastaFile, vector<Sequence>& seqs) {
	fstream file;
	MstUtils::openFile(file, fastaFile, fstream::in, "SeqTools::readFasta " + fastaFile);

	string id, seq, line;
	while (getline(file, line)) {
    line = MstUtils::trim(line);
    if (line.empty()) continue;
		if (line[0] == '>') { // identifier lines should start with '>'
			if (id.length() > 0) { // add previous (id, sequence) if it is not the first identifier (i.e. id.length > 0)
				MstUtils::assertCond((seq.length() > 0), "Sequence " + MstUtils::toString(seqs.size() + 1) + " appears to be missing", "SeqTools::readFast");
				seqs.push_back(Sequence(seq, id));
			}
			id = MstUtils::trim(line.substr(1, line.length() - 1)); // skip ">"
			MstUtils::assertCond((id.length() > 0), "Identifier for sequence " + MstUtils::toString(seqs.size() + 1) + " appears to be missing", "SeqTools::readFast");
			seq = ""; // clear for next sequence
		} else {
			seq += line; // seqeunces can be multi-line
		}
	}

	if (seq.length() > 0) {
		seqs.push_back(Sequence(seq, id));
	}
	file.close();
}

vector<Sequence> SeqTools::readFasta(const string& fastaFile) {
  vector<Sequence> seqs;
  SeqTools::readFasta(fastaFile, seqs);
  return seqs;
}

void SeqTools::readSequences(const string& seqsFile, vector<Sequence>& seqs) {
  fstream file;
	MstUtils::openFile(file, seqsFile, fstream::in, "SeqTools::readSequences " + seqsFile);
	string line;
	while (getline(file, line)) {
    line = MstUtils::trim(line);
    if (line.empty()) continue;
    seqs.push_back(Sequence(line, ""));
	}
	file.close();
}

vector<Sequence> SeqTools::readSequences(const string& seqsFile) {
  vector<Sequence> seqs;
  SeqTools::readSequences(seqsFile, seqs);
  return seqs;
}

vector<vector<int> > SeqTools::rSearch(const vector<Sequence>& seqs, mstreal idCut, mstreal a, bool verb) {
  MstUtils::assertCond((idCut >= 0) && (idCut <= 1.0), "ID cutoff value must be [0; 1]", "SeqTools::rSearch()");
  int N = seqs.size();
  vector<vector<int> > result(N);
  if (N == 0) return result;
  int L = seqs[0].length();
  int S = (int) ceil(1.0*L*idCut); // number of identities necessary to pass the cutoff
  MstTimer tim;

  // --- pick word length based on the problem
  // get bias in the sequence set to know how to compute word match expectations
  vector<int> hist(SeqTools::maxIndex(), 0);
  for (int i = 0; i < N; i++) {
    const Sequence& S = seqs[i];
    for (int j = 0; j < S.length(); j++) hist[S[j]]++;
  }
  mstreal pe = 0; // expected probability of two randomly picked amino acids from this set matching
  for (int i = 0; i < hist.size(); i++) pe += pow(hist[i]*1.0/(N*L), 2);
  mstreal d = 1; /* this parameter controls the cost of finding all sequences with
                  * a common word with a given sequence relative to the cost of
                  * comparing two sequences. Note, that because we do the former in
                  * one go for all sequences, this cost is amortized. But still, the
                  * actual (best) value to use will depend on the implementation. */
  vector<mstreal> cost(S, 0.0); // estimated search costs for each possible word length
  vector<int> Niters(S, 0.0);   // corresponding number of cycles needed
  for (int w = 1; w <= S; w++) {
    mstreal p = 1.0; // the probability of a hit with S identities being in the
                     // set of sequences with a w-common with the query sequence.
    for (int k = 0; k < w; k++) {
      p *= (S - k)*1.0/(L - k);
    }
    Niters[w-1] = ceil(log(1 - a)/log(1 - p));
    cost[w-1] = Niters[w-1]*(d + (N - 1)*pow(pe, w));
    if (verb) printf("pe = %f, w = %d, p = %f, n = %d, cost = %f\n", pe, w, p, Niters[w-1], cost[w-1]);
  }
  int minIdx; MstUtils::min(cost, 0, cost.size() - 1, &minIdx);
  int w = minIdx + 1;     // best word length to use
  int n = Niters[minIdx]; // number of repeated lookups needed to reach coverage a
  if (verb) cout << "chose word length " << w << ", and will do " << n << " cycles" << endl;

  // -- prepare C-style arrays to store sorted indices, words and words transposes (for different accept patterns)
  int* indices = new int[N];
  res_t** wordsT = new res_t*[L];
  for (int k = 0; k < w; k++) wordsT[k] = new res_t[N];
  res_t** words = new res_t*[N];
  for (int k = 0; k < N; k++) words[k] = new res_t[w];

  // --- do repeated word-based lookups
  vector<set<int> > resultSets(N);
  for (int c = 0; c < n; c++) {
    // pick random sub-set of positions
    vector<int> pos(L), wordPos(w);
    for (int i = 0; i < L; i++) pos[i] = i;
    MstUtils::shuffle(pos);
    for (int i = 0; i < w; i++) wordPos[i] = pos[i];

    if (verb) tim.start();
    // extract words and their transposes
    for (int k = 0; k < w; k++) {
      for (int i = 0; i < N; i++) {
        wordsT[k][i] = seqs[i][wordPos[k]];
        words[i][k] = wordsT[k][i];
      }
    }

    // sort sequences by word
    SeqTools::sortSequences(wordsT, indices, N, w);
    if (verb) {
      tim.stop();
      cout << "sorting took " << tim.getDuration(MstTimer::msec) << " msec" << endl;
      tim.start();
    }

    // for each sequence, look through other sequences with matching word
    int beg = 0, end = beg;
    for (int i = 1; i < N; i++) {
      if (!areWordsIdentical(words[indices[i]], words[indices[beg]], w)) {
        end = i - 1;
        // mutually compare set between beg and end
        for (int j = beg; j <= end - 1; j++) {
          const Sequence& seqI = seqs[indices[j]];
          for (int k = j + 1; k <= end ; k++) {
            if (areSequencesWithinID(seqI, seqs[indices[k]], S)) {
              resultSets[indices[j]].insert(indices[k]);
              resultSets[indices[k]].insert(indices[j]);
            }
          }
        }
        beg = i;
        end = i;
      }
    }
    if (verb) {
      tim.stop();
      cout << "comparing took " << tim.getDuration(MstTimer::msec) << " msec" << endl;
    }
  }

  // -- clea up
  delete[] indices;
  for (int k = 0; k < w; k++) delete[] wordsT[k];
  delete[] wordsT;
  for (int k = 0; k < N; k++) delete[] words[k];
  delete[] words;

  // --- extract/return results
  for (int i = 0; i < N; i++) {
    int k = 0;
    result[i].resize(resultSets[i].size());
    for (auto it = resultSets[i].begin(); it != resultSets[i].end(); ++it, ++k) result[i][k] = *it;
  }
  return result;
}

bool SeqTools::areWordsIdentical(res_t* wordA, res_t* wordB, int L) {
  for (int i = 0; i < L; i++) {
    if (wordA[i] != wordB[i]) return false;
  }
  return true;
}


vector<int> SeqTools::sortSequences(const vector<Sequence>& _seqs) {
  int N = _seqs.size();
  if (N == 0) return vector<int>();
  int L = _seqs[0].length();
  int i, j;

  // take transpose of sequence list as we will need to access columns
  // (this speeds things up A LOT)
  res_t** seqs = new res_t*[L];
  for (i = 0; i < L; i++) {
    seqs[i] = new res_t[N];
    for (j = 0; j < N; j++) seqs[i][j] = _seqs[j][i];
  }

  // sort
  int* sortedIndicesArray = new int[N];
  sortSequences(seqs, sortedIndicesArray, N, L);
  vector<int> sortedIndices(N);
  for (int i = 0; i < N; i++) sortedIndices[i] = sortedIndicesArray[i];
  delete[] sortedIndicesArray;

  // cleanup
  for (i = 0; i < L; i++) delete[] seqs[i];
  delete[] seqs;

  return sortedIndices;
}

void SeqTools::sortSequences(res_t** seqs, int* sortedIndices, int N, int L) {
  int i, j, k, c, idx, sz;

  // radix sort
  vector<vector<int> > buckets(SeqTools::maxIndex());
  for (i = 0; i < N; i++) sortedIndices[i] = i;
  for (k = L-1; k >= 0; k--) {
    // sort list elements into buckets
    res_t* col = seqs[k];
    for (i = 0; i < N; i++) {
      idx = sortedIndices[i];
      buckets[col[idx]].push_back(idx);
    }

    // rebuild list from buckets
    c = 0;
    for (i = 0; i < buckets.size(); i++) {
      vector<int>& bucket = buckets[i];
      int sz = bucket.size();
      for (j = 0; j < sz; j++, c++) {
        sortedIndices[c] = bucket[j];
      }
    }

    // empty buckets
    for (i = 0; i < buckets.size(); i++) {
      buckets[i].resize(0); // should not change capacity, so futer push_back() will be fast
    }
  }
}

bool SeqTools::sortSequencesTest(vector<Sequence> seqs) {
  vector<int> indsStd = MstUtils::sortIndices(seqs);
  vector<int> inds = sortSequences(seqs);
  if (indsStd.size() != inds.size()) return false;
  for (int i = 0; i < inds.size(); i++) {
    if (seqs[inds[i]] != seqs[indsStd[i]]) return false;
  }
  return true;
}

bool SeqTools::areSequencesWithinID(const Sequence& seqA, const Sequence& seqB, int numID) {
  int nRem = numID, L = seqA.size();
  for (int i = 0; (i < L - nRem + 1) && (nRem > 0); i++) {
    if (seqA[i] == seqB[i]) nRem--;
  }
  return nRem <= 0; // < 0 should never happen, but just in case
}

bool SeqTools::areSequencesWithinID(const Sequence& seqA, const Sequence& seqB, mstreal idCut) {
  return areSequencesWithinID(seqA, seqB, (int) ceil(seqA.size() * idCut));
}

int SeqTools::sequenceIdentity(const Sequence& seqA, const Sequence& seqB) {
  int numID = 0, L = seqA.size();
  for (int i = 0; i < L; i++) {
    if (seqA[i] == seqB[i]) numID++;
  }
  return numID;
}

mstreal SeqTools::complexity(const vector<int>& seq, int mutSite, int mutAA) {
  if (seq.empty()) return 0;

  // get counts of all letters
  int mi = 9999, ma = -1;
  for (int aa : seq) {
    mi = min(mi, aa);
    ma = max(ma, aa);
  }
  vector<int> counts(ma - mi + 1, 0);
  for (int aa : seq) counts[aa - mi]++;
  mstreal e = exp(1.0);
  int L = seq.size();

  // if need mutational difference, fewer calculations
  if (mutSite >= 0) {
    mstreal del = 0;
    int k = seq[mutSite] - mi;
    del += 0.5*log(2*M_PI*counts[k]) + counts[k]*log(counts[k]/e);
    if (counts[k] > 1) del -= 0.5*log(2*M_PI*(counts[k] - 1)) + (counts[k] - 1)*log((counts[k] - 1)/e);
    k = mutAA - mi;
    del -= 0.5*log(2*M_PI*(counts[k] + 1)) + (counts[k] + 1)*log((counts[k] + 1)/e);
    if (counts[k] > 0) del += 0.5*log(2*M_PI*counts[k]) + counts[k]*log(counts[k]/e);
    return del;
  }

  // compute complexity using Stirling's approximation of the factorial
  mstreal C = 0.5*log(2*M_PI*L) + L*log(L/e);
  for (int i = 0; i < counts.size(); i++) {
    if (counts[i]) C -= 0.5*log(2*M_PI*counts[i]) + counts[i]*log(counts[i]/e);
  }

  return C;
}

/* ------------ Sequence ------------ */

Sequence::Sequence(const Chain& C) {
  seq.resize(C.residueSize());
  for (int i = 0; i < C.residueSize(); i++) seq[i] = SeqTools::aaToIdx(C[i].getName());
  name = C.getID();
}

Sequence::Sequence(const Structure& S) {
  vector<Residue*> residues = S.getResidues();
  seq.resize(residues.size());
  for (int i = 0; i < residues.size(); i++) seq[i] = SeqTools::aaToIdx(residues[i]->getName());
  name = S.getName();
}

Sequence::Sequence(const string& _seq, const string& _name, const string& delim) {
  vector<string> seqChar = MstUtils::split(_seq, delim);
  seq.resize(seqChar.size());
  for (int i = 0; i < seqChar.size(); i++) seq[i] = SeqTools::aaToIdx(seqChar[i]);
  name = _name;
}

Sequence::Sequence(const Sequence& S) {
  name = S.name;
  seq = S.seq;
}

Sequence::Sequence(int L, const string& _name) {
  seq.resize(L, SeqTools::unknownIdx());
  name = _name;
}

void Sequence::appendResidue(const string& aa) {
  seq.push_back(SeqTools::aaToIdx(aa));
}

string Sequence::toString(bool triple, const string& delim) const {
	string s;
	for (int i = 0; i < seq.size(); i++) {
		if (triple) {
			s += SeqTools::idxToTriple(seq[i]);
		} else {
			s += SeqTools::idxToSingle(seq[i]);
		}
		if (i < seq.size() - 1) s += delim;
	}
	return s;
}

vector<string> Sequence::toStringVector(bool triple) const {
	vector<string> s(size());
	for (int i = 0; i < seq.size(); i++) {
		if (triple) {
			s[i] = SeqTools::idxToTriple(seq[i]);
		} else {
			s[i] = SeqTools::idxToSingle(seq[i]);
		}
	}
	return s;
}

string Sequence::getResidue(int i, bool triple) const {
  if (triple) return SeqTools::idxToTriple(seq[i]);
  return SeqTools::idxToSingle(seq[i]);
}

Sequence Sequence::subSequence(const vector<int>& inds) const {
  Sequence sub(inds.size());
  for (int i = 0; i < inds.size(); i++) sub[i] = (*this)[inds[i]];
  return sub;
}

Sequence Sequence::extractRange(int min, int max) const {
  if (max < min) MstUtils::error("max ("+MstUtils::toString(max)+") must be > than min ("+MstUtils::toString(min)+")","Sequence::extractRange");
  if ((max >= seq.size())||(min < 0)) MstUtils::error("max ("+MstUtils::toString(max)+") and min ("+MstUtils::toString(min)+") must be in the range [0,L)","Sequence::extractRange");
  Sequence sub(max-min+1);
  for (int i = 0; i+min <= max; i++) sub[i] = (*this)[i+min];
  return sub;
}

void Sequence::resize(int newLen, res_t newIdx) {
  if (newIdx == -1) seq.resize(newLen, SeqTools::gapIdx());
  else seq.resize(newLen, newIdx);
}

void Sequence::write(ostream& _os) const {
  MstUtils::writeBin(_os, getName());
  MstUtils::writeBin(_os, length());
  for (int i = 0; i < length(); i++) MstUtils::writeBin(_os, (*this)[i]);
}

void Sequence::read(istream& _is) {
  MstUtils::readBin(_is, name);
  int len; MstUtils::readBin(_is, len);
  seq.resize(len);
  for (int i = 0; i < seq.size(); i++) MstUtils::readBin(_is, seq[i]);
}

ostream& MST::operator<<(ostream &_os, const Sequence& _seq) {
  _os << "> " << _seq.getName() << endl;
  int k = 0;
  for (int i = 0; i < _seq.length(); i++) {
    if (k >= 40) { _os << endl; k = 0; }
    _os << SeqTools::idxToSingle(_seq[i]); k++;
  }
  return _os;
}

bool Sequence::operator==(const Sequence& other) const {
  if (length() != other.length()) return false;
  for (int i = 0; i < length(); i++) {
    if ((*this)[i] != other[i]) return false;
  }
  return true;
}

bool Sequence::operator!=(const Sequence& other) const {
  if (length() != other.length()) return true;
  for (int i = 0; i < length(); i++) {
    if ((*this)[i] != other[i]) return true;
  }
  return false;
}
