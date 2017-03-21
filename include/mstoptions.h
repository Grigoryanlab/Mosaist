#ifndef _MSTOPTIONS_H
#define _MSTOPTIONS_H

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <map>
#include <set>
#include "msttypes.h"

using namespace std;
using namespace MST;

/* Options parsing/usage statement class. */
class MstOptions {
  public:
    MstOptions() { w = 80; p1 = 3; p2 = p1+8; }
    MstOptions(int argc, char** argv) : MstOptions() { setOptions(argc, argv); }

    // formatting of usage information
    void addOption(string opt, string info, bool req = false);
    void setTitle(string _title) { title = _title; }
    string usage();
    void setUsageWidth(int _w) { w = _w; }
    void setUsageParOffset(int _p1) { p1 = _p1; }
    void setUsageSecondParOffset(int _p2) { p2 = _p2; }
    string formatOptInfo(string opt, string mes, int _p1 = -1, int _p2 = -1);

    // parsing command-line options
    void setOptions(int argc, char** argv);
    bool isInt(const string& opt);
    int getInt(const string& opt, int defVal = 0);
    bool isReal(const string& opt);
    real getReal(const string& opt, real defVal = 0.0);
    string getString(const string& opt, string defVal = "");
    bool getBool(const string& opt) { return isGiven(opt); }
    bool isGiven(const string& opt);

  private:
    int w, p1, p2;
    string execName, title;
    map<string, string> givenOptions;
    vector<string> options;
    vector<string> optionsInfo;
    set<string> required;
};

#endif
