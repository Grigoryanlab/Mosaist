#include "msttypes.h"
#include "dtermen.h"
#include <chrono>

int main(int argc, char *argv[]) {
  dTERMen D("testfiles/dtermen.conf");
  cout << "omega energy for ALA with omega 0 is: " << D.bbOmegaEner(0, SeqTools::aaToIdx("ALA")) << endl;
  cout << "omega energy for ALA with omega 200 is: " << D.bbOmegaEner(200, SeqTools::aaToIdx("ALA")) << endl;
  cout << "omega energy for ALA with omega -5 is: " << D.bbOmegaEner(-5, SeqTools::aaToIdx("ALA")) << endl;
  cout << "omega energy for LEU with omega -180 is: " << D.bbOmegaEner(-180, SeqTools::aaToIdx("ALA")) << endl;
  cout << "env energy for GLU with freedom 0.2 is: " << D.envEner(0.2, SeqTools::aaToIdx("GLU")) << endl;
  cout << "env energy for GLU with freedom 1.0 is: " << D.envEner(1.0, SeqTools::aaToIdx("GLU")) << endl;
}
