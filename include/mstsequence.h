#include "msttypes.h"

class Sequence {
  public:
    Sequence() {}
//    Sequence(const Chain& C);
//    Sequence(const Residue& R);
//    Sequence(const System& S);
//    Sequence(string seq, int letWidth = 1, char sep = ' ');
//    ~Sequence() {}

    // A kind of a constructor for all the static variables, so they do not need
    // an instance to exist. This will be called automatically.
    static bool initConstants();

    static string tripleToSingle(const string triple, string del = " ");
    static string singleToTriple(const string single, string del = "");
    static vector<int> seqToIdx(const string single, string del = "");
    static int aaToIdx(const string aa);
    static string idxToTriple(int idx);
    static string idxToSingle(int idx);
    static string toTriple(const string aa);
    static string toSingle(const string aa);

  private:
    vector<int> _seq;
    static vector<string> aa1, aa3;
    static map<string, int> aa3ToIdx, aa1ToIdx;
    static map<int, string> idxToAA3, idxToAA1;
    static bool initialized;
};
