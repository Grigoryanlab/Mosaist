#include "mstsequence.h"

// statics must be defined in the cpp file
vector<string> Sequence::aa1, Sequence::aa3;
map<string, int> Sequence::aa3ToIdx, Sequence::aa1ToIdx;
map<int, string> Sequence::idxToAA3, Sequence::idxToAA1;
bool initialized = Sequence::initConstants();

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

bool Sequence::initConstants() {
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

  for (int i = 0; i < aa3.size(); i++) {
    aa3ToIdx[aa3[i]] = i;
    // single single-letter code is ambiguous, take the first amino acid
    // when going from single-letter code to index
    if (aa1ToIdx.find(aa1[i]) == aa1ToIdx.end()) aa1ToIdx[aa1[i]] = i;
    idxToAA3[i] = aa3[i];
    idxToAA1[i] = aa1[i];
  }
  return true;
}
