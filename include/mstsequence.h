#ifndef _MSTSEQUENCE_H
#define _MSTSEQUENCE_H

#include "msttypes.h"
#include <unordered_map>
#include <unordered_set>

namespace MST {

typedef short res_t;
class Sequence {
  public:
    Sequence() {}
    Sequence(const Structure& S);
    Sequence(const Chain& C);
    // allows one to mix 3- and 1-letter codes, with a delimiter
    Sequence(const string& _seq, const string& _name = "", const string& delim = "");
    Sequence(const Sequence& S);
    Sequence(const vector<res_t> _seq, const string& _name = "") { seq = _seq; name = _name; }
    Sequence(int L, const string& _name = "");

    string getName() const { return name; }
    void setName(const string& _name) { name = _name; }
    string toString(bool triple = false, const string& delim = "") const;
    vector<string> toStringVector(bool triple = false) const;
    string getResidue(int i, bool triple = false) const;
    res_t& operator[] (int i) { return seq[i]; }
    res_t operator[] (int i) const { return seq[i]; }
    Sequence subSequence(const vector<int>& inds) const;
    int length() const { return seq.size(); }
    int size() const { return seq.size(); }
    void appendResidue(const string& aa);
    void appendResidue(res_t aai) { seq.push_back(aai); }
    void resize(int newLen, res_t newIdx = -1);
    void write(ostream& _os) const; // write Sequence to a binary stream
    void read(istream& _is);  // read Sequence from a binary stream

    // to allow Sequences to be used in maps, sets, etc. and to be sorted
    friend bool operator<(const Sequence& s1, const Sequence& s2) {
      if (s1.length() != s2.length()) return s1.length() < s2.length();
      for (int i = 0; i < s1.length(); i++) {
        if (s1[i] != s2[i]) return s1[i] < s2[i];
      }
      return false;
    }

    friend ostream & operator<<(ostream &_os, const Sequence& _seq);
    bool operator==(const Sequence& other) const;
    bool operator!=(const Sequence& other) const;

  private:
    vector<res_t> seq;
    string name;
};
ostream & operator<<(ostream &_os, const Sequence& _seq); // this is just to silence a silly compiler warning

class SeqTools {
  public:
    // A kind of a constructor for all the static variables, so they do not need
    // an instance to exist. This will be called automatically.
    static bool initConstants();

    static string tripleToSingle(const string& triple, const string& del = " ");
    static string singleToTriple(const string& single, const string& del = "");
    static vector<res_t> seqToIdx(const string& single, const string& del = "");
    static res_t aaToIdx(const string& aa);
    static res_t unknownIdx() { return _unkIdx; }
    static bool isUnknown(const string& aa) { return aaToIdx(aa) == _unkIdx; }
    static res_t gapIdx() { return _gapIdx; }
    static string idxToTriple(res_t idx);
    static string idxToSingle(res_t idx);
    static res_t maxIndex() { return idxToAA1.size() - 1; }
    static string toTriple(const string& aa);
    static string toSingle(const string& aa);
    static vector<Sequence> readFasta(const string& fastaFile);
    static void readFasta(const string& fastaFile, vector<Sequence>& seqs);
    static vector<Sequence> readSequences(const string& seqsFile);
    static void readSequences(const string& seqsFile, vector<Sequence>& seqs);
    static vector<string> getAA1() {return aa3;}
    static vector<string> getAA3() {return aa1;}

    /* For these two functions, seqA and seqB must be of the same length. This
     * is not checked, for efficiency. The comparison is smart, in that it stops
     * when it is clear that the identity level cannot be reached with the rest
     * of the sequences. The first function looks for a given number of identities
     * and the second one for a given fraction identity. */
    static bool areSequencesWithinID(const Sequence& seqA, const Sequence& seqB, int numID);
    static bool areSequencesWithinID(const Sequence& seqA, const Sequence& seqB, mstreal idCut);

    /* Performs an all-by-all sequence identity search (ungapped) using a
     * randomized algorithm inspired by the USEARCH method (Robert C. Edgar,
     * Bioinformatics, v. 26 (19), 2010, 2460–2461). The idea is that if we look
     * for all sequences that have a particular word (i.e., a sub-sequence of
     * not necessarily contiguous positions) in common with a query sequence,
     * there is a (calculatable) chance that any match to the query (passing the
     * sequence identity threshold) would be within this set. That chance is,
     * of course, not 100%, but it is some computable quantity (at least, the
     * expectation is computable). So if we repeat this procedure some number of
     * times (i.e., find all sequences with a common word with the query, where
     * the positions for the word are randomly generated each time, and compare
     * the query to each of these sequences), we will eventually reach a pretty
     * high chance of recovering all matches. This algorithm chooses a word size
     * and number of iterations to get to any desired level of accuracy a (i.e.,
     * expected fracton of matches found) while minimizing the running time. */
    static vector<vector<int> > rSearch(const vector<Sequence>& seqs, mstreal idCut, mstreal a = 0.99, bool verb = false);

    /* A fast sort algorithm specialized for sequence data. Uses radix sort,
     * since residue values have an a priori known maximum value. Sorts in
     * ascending order by residue, treating residues earlier in the sequence as
     * more significant "digits". All sequences must be of the same length (not
     * explicitly checked!) */
    static vector<int> sortSequences(const vector<Sequence>& _seqs);

  protected:
    // helper function for sequence sorting
    static bool areWordsIdentical(res_t* wordA, res_t* wordB, int L);

    // Sort helper function. Implements radix sort for sequences. Parameter seqs
    // is expected to be an L x N 2D array, where L is sequence length and N is
    // the number of sequences. Note that this is "transpose" of a typical array
    // of sequences (this speeds things A LOT due to memory access patterns).
    static void sortSequences(res_t** seqs, int* sortedIndices, int N, int L);

    static bool sortSequencesTest(vector<Sequence> seqs);

  private:
    static vector<string> aa1, aa3;
    static map<string, res_t> aa3ToIdx, aa1ToIdx;
    static vector<string> idxToAA3, idxToAA1;
    static bool initialized;
    static res_t _unkIdx, _gapIdx;
};


}

/* specialize std::hash for objects of type Sequence, to be able to use
 * unordered_map and other associative STL containters with them. */
namespace std {
  template <>
  struct hash<MST::Sequence> {
    std::size_t operator()(const MST::Sequence& S) const {
      std::size_t ret = 0;
      for (int i = 0; i < S.length(); i++) {
        // Shifting by 5 bits left is akin to multiplying by 32, which is a
        // decent base for sequence elements (typically expect about 20 or so).
        // Of course, there will be collisions for long sequences due to
        // wrapping of the size_t value.
        ret = ret << 5;
        ret = ret + S[i];
      }
      return ret;
    }
  };
}

#endif
