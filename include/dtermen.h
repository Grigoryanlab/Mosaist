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
