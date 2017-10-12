#include "dtermen.h"

dTERMen::dTERMen(const string& configFile) {
  ...
  F.readDatabase(...);
  F.setRedundancyRemoval(0.30);
}
