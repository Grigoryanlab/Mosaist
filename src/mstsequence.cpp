#include "mstsequence.h"
using namespace MST;

/* ------------ SeqTools ------------ */

// statics must be defined in the cpp file
vector<string> SeqTools::aa1, SeqTools::aa3;
map<string, int> SeqTools::aa3ToIdx, SeqTools::aa1ToIdx;
vector<string> SeqTools::idxToAA3, SeqTools::idxToAA1;
bool initialized = SeqTools::initConstants();

int SeqTools::aaToIdx(const string aa) {
  if (aa.length() == 1) {
    MstUtils::assert(aa1ToIdx.find(aa) != aa1ToIdx.end(), "unrecognized single-letter amino acid '" + aa + "'", "SeqTools::aaToIdx");
    return aa1ToIdx[aa];
  } else if (aa.length() == 3) {
    MstUtils::assert(aa3ToIdx.find(aa) != aa3ToIdx.end(), "unrecognized three-letter amino acid '" + aa + "'", "SeqTools::aaToIdx");
    return aa3ToIdx[aa];
  }
  MstUtils::error("uknown amino acid '" + aa + "'", "SeqTools::aaToIdx");
  return 0; // for the compiler to be happy
}

string SeqTools::idxToTriple(int idx) {
  MstUtils::assert((idx >= 0) && (idx < idxToAA3.size()), "unknown amino-acid index '" + MstUtils::toString(idx) + "'", "SeqTools::idxToTriple");
  return idxToAA3[idx];
}

string SeqTools::idxToSingle(int idx) {
  MstUtils::assert((idx >= 0) && (idx < idxToAA1.size()), "unknown amino-acid index '" + MstUtils::toString(idx) + "'", "SeqTools::idxToSingle");
  return idxToAA1[idx];
}

string SeqTools::toTriple(const string aa) {
  return idxToTriple(aaToIdx(aa));
}

string SeqTools::toSingle(const string aa) {
  return SeqTools::idxToSingle(SeqTools::aaToIdx(aa));
}

string SeqTools::tripleToSingle(const string triple, string del) {
  vector<string> seq = MstUtils::split(triple, del);
  string single;
  for (int i = 0; i < seq.size(); i++) {
    single += SeqTools::toSingle(seq[i]);
  }
  return single;
}

string singleToTriple(const string single, string del) {
  vector<string> seq = MstUtils::split(single, del);
  string triple;
  for (int i = 0; i < seq.size(); i++) {
    triple += SeqTools::toTriple(seq[i]);
    if (i != seq.size() - 1) triple += " ";
  }
  return triple;
}

vector<int> seqToIdx(const string seqstr, string del) {
  vector<string> seq = MstUtils::split(seqstr, del);
  vector<int> indices(seq.size());
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

  idxToAA3.resize(aa3.size());
  idxToAA1.resize(aa1.size());
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

void SeqTools::readFasta(const string& fastaFile, vector<Sequence>& seqs) {
	fstream file;
	MstUtils::openFile(file, fastaFile, fstream::in, "SeqTools::readFasta " + fastaFile);

	string id, seq, line;
	while (getline(file, line)) {
    line = MstUtils::trim(line);
    if (line.empty()) continue;
		if (line[0] == '>') { // identifier lines should start with '>'
			if (id.length() > 0) { // add previous (id, sequence) if it is not the first identifier (i.e. id.length > 0)
				MstUtils::assert((seq.length() > 0), "Sequence " + MstUtils::toString(seqs.size() + 1) + " appears to be missing", "SeqTools::readFast");
				seqs.push_back(Sequence(seq, id));
			}
			id = MstUtils::trim(line.substr(1, line.length() - 1)); // skip ">"
			MstUtils::assert((id.length() > 0), "Identifier for sequence " + MstUtils::toString(seqs.size() + 1) + " appears to be missing", "SeqTools::readFast");
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


/* ------------ Sequence ------------ */

Sequence::Sequence(Chain& C) {
  seq.resize(C.residueSize());
  for (int i = 0; i < C.residueSize(); i++) seq[i] = SeqTools::aaToIdx(C[i].getName());
  name = C.getID();
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

string Sequence::toString(bool triple, const string& delim) {
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
