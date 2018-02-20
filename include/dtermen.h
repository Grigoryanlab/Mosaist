#ifndef _DTERMEN_H
#define _DTERMEN_H

#include "msttypes.h"

class TERMen {
  public:
    string type;
    mstreal value;
};

class dTERMen {
  public:
    // Implementation approach:
    // 1. DONE: get phi, psi, and environment values into the database (fasstDB does this now)
    // 2. TODO: perform all-by-all local sequence window comparison within the database, using the
    //    uclust-like algorithm (not guaranteed, but pretty good for close homology), and get
    //    two things for every position: its multiplicity (i.e., the number of times closely
    //    homologous local sequences are found in the database) and the local cluster index
    //    (i.e., such that all closely homologous local sequence windows belong to the same
    //    cluster).
    // 3. TODO: make a modification in FASST, so that it can optionally use the "local cluster
    //    index" property to remove redundancy (no need to calculate alighments, so very fast),
    //    if the property is available within the database.
    // 4. TODO: implement on-the-fly backbone and environment potential calculation (upon reading
    //    the FASST database), by accounting for multiplicity of each residue.
    // TODO:
    // 1. get locations of phi, psi, omega, environment potentials and read these
    // 2. get locations of phi, psi, omega, env databases and read these
    // 3. get location of FASST database and initialize a built-in FASST object
    dTERMen(const string& configFile);

    vector<TERMen> selfEnergy(Residue* R, vector<Residue*> C);

  private:
    FASTT F;
};

#endif
