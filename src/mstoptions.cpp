#include "mstoptions.h"

void MstOptions::addOption(string opt, string info, bool req) {
  options.push_back(opt);
  optionsInfo.push_back(info);
  if (req) required.insert(opt);
  // if long option names are specified, move the usage info tab over
  if (p1 + opt.length() + 3 > p2) p2 = p1 + opt.length() + 3;
}

string MstOptions::usage() {
  string usageString;
  usageString = "\n" + formatOptInfo("", title, 0, 0) + "\n";
  for (int i = 0; i < options.size(); i++) {
    usageString += formatOptInfo(options[i], optionsInfo[i]) + "\n";
  }
  return usageString;
}

string MstOptions::formatOptInfo(string opt, string mes, int _p1, int _p2) {
  if (_p1 < 0) _p1 = p1;
  if (_p2 < 0) _p2 = p2;

  // first print the name of the option
  string text(_p1, ' ');
  if (!opt.empty()) text += "--" + opt;
  if (_p2 > text.size()) text += string(_p2 - text.size(), ' ');

  // next print the description text
  int i = 0, k, L = text.size(), n;
  while (i < mes.size()) {
    k = mes.find_first_of(" ", i);
    if (k == string::npos) k = mes.size();
    n = k - i;
    if ((L + n >= w) && (L > 0)) { text += "\n" + string(_p2, ' '); L = _p2; }
    text += mes.substr(i, n) + " ";
    L += n + 1;
    i = k+1;
  }
  return text;
}

void MstOptions::setOptions(int argc, char** argv) {
  execName.clear();
  givenOptions.clear();
  set<string> missingRequired = required;
  if (argc <= 0) return;
  execName = argv[0];
  for (int i = 1; i < argc; ) {
    string token = argv[i];
    string nextToken = ((i < argc - 1) ? argv[i+1] : "");
    if (token.find("--") != 0) MstUtils::error("could not parse options, could not understand the token '" + token + "'", "MstOptions::setOptions");
    token = token.substr(2);
    if (nextToken.find("--") != 0) {
      givenOptions[token].push_back(nextToken);
      i += 2;
    } else {
      givenOptions[token].push_back("");
      i++;
    }
    missingRequired.erase(token);
  }
  if (missingRequired.size() != 0) {
    cerr << usage() << endl;
    for (auto it = missingRequired.begin(); it != missingRequired.end(); ++it) cerr << "missing required option: " << *it << endl;
    cerr << endl;
    MstUtils::error("not all required options were specified!", "MstOptions::setOptions");
  }
}

int MstOptions::getInt(const string& opt, int defVal, int idx) {
  if (givenOptions.find(opt) == givenOptions.end()) return defVal;
  return MstUtils::toInt(givenOptions[opt][idx]);
}

mstreal MstOptions::getReal(const string& opt, mstreal defVal, int idx) {
  if (givenOptions.find(opt) == givenOptions.end()) return defVal;
  return MstUtils::toReal(givenOptions[opt][idx]);
}

bool MstOptions::isInt(const string& opt, int idx) {
  if (givenOptions.find(opt) == givenOptions.end()) return false;
  return MstUtils::isInt(givenOptions[opt][idx]);
}

bool MstOptions::isReal(const string& opt, int idx) {
  if (givenOptions.find(opt) == givenOptions.end()) return false;
  return MstUtils::isReal(givenOptions[opt][idx]);
}

string MstOptions::getString(const string& opt, string defVal, int idx) {
  if (givenOptions.find(opt) == givenOptions.end()) return defVal;
  return givenOptions[opt][idx];
}

bool MstOptions::getBool(const string& opt, int idx) {
  if (givenOptions.find(opt) == givenOptions.end()) return false;
  if (!MstUtils::isInt(givenOptions[opt][idx])) return false;
  return (MstUtils::toInt(givenOptions[opt][idx]) != 0);
}

bool MstOptions::isGiven(const string& opt) const {
  if (givenOptions.find(opt) == givenOptions.end()) return false;
  return true;
}

int MstOptions::timesGiven(const string& opt) const {
  if (givenOptions.find(opt) == givenOptions.end()) return 0;
  return givenOptions.at(opt).size();
}
