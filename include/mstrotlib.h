#ifndef _MSTLIB_H
#define _MSTLIB_H

#include "msttypes.h"
#include "msttransforms.h"
#include <algorithm> // needed for sort

namespace MST {

class RotamerLibrary {
  public:
    RotamerLibrary() {}
    RotamerLibrary(string rotLibFile);
    ~RotamerLibrary();

    void readRotamerLibrary(string rotLibFile);

    /* returns the index of the bin into which the given phi/psi value combination
     * maps for the given amino acid. */
    int getBackboneBin(string aa, real phi, real psi, bool assumeDefault = true);

    /* finds the phi/psi values corresponding to the given bin index for the given amino acid */
    pair<real, real> getBinPhiPsi(string aa, int bi);

    /* places the specified rotamer into the given Residue. NOTE: the original residue
     * is modified, with some of its atoms potentially destroyed (as needed), so if you
     * want to be able to go back to the wild type, first make a copy of the residue
     * before placing the rotamer. If a destination residue is specified, rather than
     * modifying the target residue in place, it modifies the destination residue.
     * Expects that the destination residue will either be empty (i.e., no atoms) OR
     * will be filled with precisely the correct atoms for the amino acid. The latter
     * corresponds to the case when the destination residue was already previously built
     * by this function, with perhaps a different rotamer; this case is for efficiency. */
    void placeRotamer(Residue& res, string aa, int rotIndex, Residue* dest_ptr = NULL, bool strict = false);

    /* decides whether the atom is a backbone atom basded on the name */
    static bool isBackboneAtom(string atomName);
    static bool isBackboneAtom(Atom& atom) { return isBackboneAtom(atom.getName()); }
    static bool isBackboneAtom(Atom* atom) { return isBackboneAtom(atom->getName()); }
    static bool isHydrogen(string atomName);
    static bool isHydrogen(Atom& atom) { return isHydrogen(atom.getName()); }
    static bool isHydrogen(Atom* atom) { return isHydrogen(atom->getName()); }

    int numberOfRotamers(string aa, real phi = Residue::badDihedral, real psi = Residue::badDihedral, bool strict = false);
    real rotamerProbability(string aa, int ri, real phi = Residue::badDihedral, real psi = Residue::badDihedral, bool strict = false);
    vector<string> availableAminoAcids() { return keys(rotamers); }

  protected:
    /* given an array of angles, stored in ascending order (i.e., in the counter-clockwise
     * direction), find the array index with the angle closest to the given angle */
    int findClosestAngle(vector<real>& array, real value);

    /* make newAtoms be a vector of atoms corresponding to the given rotamer, upon
     * transformation according to the given Transform. */
    void transformRotamerAtoms(Transform& T, Residue& rots, int rotIndex, vector<Atom*>& newAtoms);

    // computes the difference between two angles, choosing the closest direction
    // (i.e., either clockwise, indicated by a negative difference or counter-clockwise,
    // indicated by a positive difference). The order of subtraction is a - b.
    real angleDiff(real a, real b);

    // map a given angle, in degrees to the "standard" range of [-180, 180)
    real angleToStandardRange(real angle);

    // return a vector of keys given a map
    template<class T1, class T2>
    vector<T1> keys(map<T1, T2>& mymap, bool sorted = false) {
      vector<T1> vec;
      for (typename map<T1, T2>::iterator it = mymap.begin(); it != mymap.end(); ++it) {
        vec.push_back(it->first);
      }
      if (sorted) {
        sort(vec.begin(), vec.end());
      }
      return vec;
    }

  private:
    /* the map is keyed by amino-acid name, the vector goes over phi/psi bins, and
     * the inner Residue object stores atom information and coordinates for all
     * rotamers of the amino acid in the backbone bin. The first rotamer is stored
     * in the main coordinates of the Residue's constituent atoms, and the remaining
     * rotamers, if any, in the alternative coordinates. */
    map<string, vector<Residue*> > rotamers;

    /* rotamer probabilities. As above, the map is keyed by amino acid, the outer vector
     * goes over phi/psi bins and the inner vector over rotamers. */
    map<string, vector<vector<real> > > prob;

    /* for a given amino acid, binFreq[aa] stores the frequencies of each phi/psi bin. These
     * are stored as reals, so they can be either counts (i.e., number of occurrences) as with
     * Dunbrack's rotamer library or true frequencies (i.e., probabilities). */
    map<string, vector<real> > binFreq;

    /* the default phi/psi bin for each amino acid; assumed in the absence of a valid phi/psi pair. */
    map<string, int> defaultBin;

    /* here the structure is similar, except rather than storing a Residue object with
     * all the rotamers, here we store a vector<vector<real> > that defines the chi
     * angles of each rotamer. The outer vector is over rotamers and the inner over chi
     * angles (e.g., chi1 through chi4). So,
     * chi["ARG"][i][13][2].first is the chi3 value for the 14th rotamer of ARG in phi/psi bin i
     * and
     * chi["ARG"][i][13][2].second -- is the standard deviation around this value (sigma). */
    map<string, vector<vector<vector<pair<real, real> > > > > chi;

    /* definitions of chi angles (via atom names). The key is the amino acid name, the
     * outer vector is over chi angles and the inner vector is over atoms comprising each
     * chi angle. E.g., chidef["ASP"][1][0] is the first atom of chi2 for ASP. */
    map<string, vector<vector<string> > > chidef;

    /* phi/psi bins are required to be on a grid, but the grid lines do not need to be
     * unifirm, do not need to be the same between phi and psi, and can vary between
     * different amino acids. The two variables below store where the grid lines in phi
     * and psi lie, such that the total number of bins for amino-acid aa is
     * binPhiCenters[aa].size() * binPsiCenters[aa].size() */
    map<string, vector<real> > binPhiCenters;
    map<string, vector<real> > binPsiCenters;
};

}
#endif
