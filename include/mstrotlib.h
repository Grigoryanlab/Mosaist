#ifndef _MSTLIB_H
#define _MSTLIB_H

#include "msttypes.h"

namespace MST {

class RotamerLibrary {
  public:
    RotamerLibrary();
    RotamerLibrary(string rotLibFile);
    ~RotamerLibrary();
  
    void readRotamerLibrary(string rotLibFile);

    /* returns the index of the bin into which the given phi/psi value combination
     * maps for the given amino acid. */
    int getBackboneBin(string aa, real phi, real psi);

    /* places the specified rotamer into the given Residue. NOTE: the original residue,
     * along with its atoms are destroyed, so if you want to be able to go back to
     * the wild type, first make a copy of the residue before placing the rotamer. */
    bool placeRotamer(Residue& res, string aa, int rotIndex, bool strict = true);

  protected:
    /* given an array of angles, stored in ascending order (i.e., in the counter-clockwise
     * direction), find the array index with the angle closest to the given angle */
    int findClosestAngle(vector<real>& array, real value);

    // computes the difference between two angles, choosing the closest direction
    // (i.e., either clockwise, indicated by a negative difference or counter-clockwise,
    // indicated by a positive difference). The order of subtraction is a - b.
    real angleDiff(real a, real b);

  private:
    /* the map is keyed by amino-acid name, the vector goes over phi/psi bins, and
     * the inner Residue object stores atom information and coordinates for all
     * rotamers of the amino acid in the backbone bin. The first rotamer is stored
     * in the main coordinates of the Residue's constituent atoms, and the remaining
     * rotamers, if any, in the alternative coordinates.
     */
    map<string, vector<Residue*> > rotamers;

    /* here the structure is similar, except rather than storing a Residue object with
     * all the rotamers, here we store a vector<vector<real> > that defines the chi
     * angles of each rotamer. The outer vector is over rotamers and the inner one over
     * chi angles (e.g., chi1 through chi4). So, :
     * chi[i]["ARG"][13][2] -- is the chi3 value for the 14th rotamer of ARG in phi/psi bin i. */
    map<string, vector<vector<vector<real> > > chi;

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
