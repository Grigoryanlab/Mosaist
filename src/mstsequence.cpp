#include "mstsequence.h"

// statics must be defined in the cpp file
vector<string> Sequence::aa1, Sequence::aa3;
map<string, int> Sequence::aa3ToIdx, Sequence::aa1ToIdx;
map<int, string> Sequence::idxToAA3, Sequence::idxToAA1;

int Sequence::aaToIdx(const string aa) {
  if (aa.length() == 1) {
    MstUtils::assert(aa1ToIdx.find(aa) != aa1ToIdx.end(), "unrecognized single-letter amino acid '" + aa + "'", "Sequence::aaToIdx");
    return aa1ToIdx[aa];
  } else if (aa.length() == 3) {
    MstUtils::assert(aa3ToIdx.find(aa) != aa3ToIdx.end(), "unrecognized three-letter amino acid '" + aa + "'", "Sequence::aaToIdx");
    return aa3ToIdx[aa];
  }
  MstUtils::error("uknown amino acid '" + aa + "'", "Sequence::aaToIdx");
  return 0; // for the compiler to be happy
}

string Sequence::idxToTriple(int idx) {
  MstUtils::assert(idxToAA3.find(idx) != idxToAA3.end(), "unknown amino-acid index '" + MstUtils::toString(idx) + "'", "Sequence::idxToTriple");
  return idxToAA3[idx];
}

string Sequence::idxToSingle(int idx) {
  MstUtils::assert(idxToAA1.find(idx) != idxToAA1.end(), "unknown amino-acid index '" + MstUtils::toString(idx) + "'", "Sequence::idxToSingle");
  return idxToAA1[idx];
}

string Sequence::toTriple(const string aa) {
  return idxToTriple(aaToIdx(aa));
}

string Sequence::toSingle(const string aa) {
  return Sequence::idxToSingle(Sequence::aaToIdx(aa));
}

string Sequence::tripleToSingle(const string triple, string del) {
  vector<string> seq = MstUtils::split(triple, del);
  string single;
  for (int i = 0; i < seq.size(); i++) {
    single += Sequence::toSingle(seq[i]);
  }
  return single;
}

string singleToTriple(const string single, string del) {
  vector<string> seq = MstUtils::split(single, del);
  string triple;
  for (int i = 0; i < seq.size(); i++) {
    triple += Sequence::toTriple(seq[i]);
    if (i != seq.size() - 1) triple += " ";
  }
  return triple;
}

vector<int> seqToIdx(const string seqstr, string del) {
  vector<string> seq = MstUtils::split(seqstr, del);
  vector<int> indices(seq.size());
  for (int i = 0; i < seq.size(); i++) {
    indices[i] = Sequence::aaToIdx(seq[i]);
  }
  return indices;
}
