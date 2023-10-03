#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <random>
#include <cmath>
#include <iostream>
#include "mstcondeg.h"
#include "mstoptions.h"
#include "mstrotlib.h"
#include "msttypes.h"
#include "msttransforms.h"
#include "mstfasst.h"
#include "mstsequence.h"

int main(int argc, char** argv) {

    MstOptions op;
    op.setTitle("Rehydration test.");
    op.addOption("a", "CA path", true);
    op.addOption("b", "CB path", true);
    op.addOption("t", "translation 3 values", true);
    op.addOption("r", "rotation 9 values", true);
    op.addOption("o", "sructOutfile", true);
    op.setOptions(argc, argv);

// load structures CA & CB, make structure C

    Structure CA(op.getString("a"));
    Structure CB(op.getString("b"));

    Structure C;
    for (int i = 0; i < CA.chainSize(); i++) {
        C.appendChain(&CA.getChain(i));
    }
    for (int i = 0; i < CB.chainSize(); i++) {
        C.appendChain(&CB.getChain(i));
    }

    // make translation & rotation matrixes

    //*** fine

    vector <mstreal> translationA = MstUtils::splitToReal(op.getString("t"), ",",false,false);

    vector <mstreal> rotationsBase = MstUtils::splitToReal(op.getString("r"), ",",false,false);

    vector<vector<mstreal>> rotationsA(3);
    for (int i = 0 ; i < 3 ; i++) {
        rotationsA[i].resize(3);
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j ++) {
            rotationsA[i][j] = rotationsBase[i*3 + j];
        }
    }

    //*** segFault by here

    // make transformation matrix, and print to check it

    Transform T(rotationsA, translationA);

    cout << T(0,0) << endl;
    cout << T(0,1) << endl;
    cout << T(0,2) << endl;
    cout << T(0,3) << endl;
    cout << T(1,0) << endl;
    cout << T(1,1) << endl;
    cout << T(1,2) << endl;
    cout << T(1,3) << endl;
    cout << T(2,0) << endl;
    cout << T(2,1) << endl;
    cout << T(2,2) << endl;
    cout << T(2,3) << endl;
    cout << T(3,0) << endl;
    cout << T(3,1) << endl;
    cout << T(3,2) << endl;
    cout << T(3,3) << endl;

    // apply to CA

    T.apply(CA);

    // save C --> open & compare to true struct!

    C.writePDB(op.getString("o"));

}