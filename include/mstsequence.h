#include "msttypes.h"

class Sequence {
  public:
    Sequence() {}
//    Sequence(const Chain& C);
//    Sequence(const Residue& R);
//    Sequence(const System& S);
//    Sequence(string seq, int letWidth = 1, char sep = ' ');
//    ~Sequence() {}

    static string tripleToSingle(const string triple, string del = " ");
    static string singleToTriple(const string single, string del = "");
    static vector<int> seqToIdx(const string single, string del = "");
    static int aaToIdx(const string aa);
    static string idxToTriple(int idx);
    static string idxToSingle(int idx);
    static string toTriple(const string aa);
    static string toSingle(const string aa);

  protected:
    static void fillArrays();

  private:
    vector<int> _seq;
    static vector<string> aa1, aa3;
    static map<string, int> aa3ToIdx, aa1ToIdx;
    static map<int, string> idxToAA3, idxToAA1;

    static class _initSequence {
      public:
        void _init() {
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
        }
    } SequenceDataInitializer;

};
