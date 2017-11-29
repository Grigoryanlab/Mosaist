#include "msttypes.h"

namespace MST {

class Sequence {
  public:
    Sequence() {}
    Sequence(Chain& C);
    // allows one to mixes of 3- and 1-letter codes, with a delimiter
    Sequence(const string& _seq, const string& _name = "", const string& delim = "");
    Sequence(const Sequence& S);
    ~Sequence() {}

    string getName() const { return name; }
    void setName(const string& _name) { name = _name; }
    string toString(bool triple = false, const string& delim = "");
    string getElement(int i, bool triple = false);
    int& operator[] (int i) { return seq[i]; }
    int operator[] (int i) const { return seq[i]; }
    int length() const { return seq.size(); }
    void appendResidue(const string& aa);
    void appendResidue(int aai) { seq.push_back(aai); }
    friend ostream & operator<<(ostream &_os, const Sequence& _seq);

  private:
    vector<int> seq;
    string name;
};

class SeqTools {
  public:
    // A kind of a constructor for all the static variables, so they do not need
    // an instance to exist. This will be called automatically.
    static bool initConstants();

    static string tripleToSingle(const string triple, string del = " ");
    static string singleToTriple(const string single, string del = "");
    static vector<int> seqToIdx(const string single, string del = "");
    static int aaToIdx(const string aa);
    static int unknownIdx() { return unkIdx; }
    static string idxToTriple(int idx);
    static string idxToSingle(int idx);
    static string toTriple(const string aa);
    static string toSingle(const string aa);
    static vector<Sequence> readFasta(const string& fastaFile);
    static void readFasta(const string& fastaFile, vector<Sequence>& seqs);

  private:
    static vector<string> aa1, aa3;
    static map<string, int> aa3ToIdx, aa1ToIdx;
    static vector<string> idxToAA3, idxToAA1;
    static bool initialized;
    static int unkIdx;
};


}
