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
  usageString = "\n" + title + "\n";
  for (int i = 0; i < options.size(); i++) {
    usageString += formatOptInfo(options[i], optionsInfo[i]) + "\n";
  }
  return usageString;
}

string MstOptions::formatOptInfo(string opt, string mes) {
  // first print the name of the option
  string text(p1, ' ');
  text += "--" + opt;
  if (p2 > text.size()) text += string(p2 - text.size(), ' ');

  // next print the description text
  int i = 0, k, L = text.size(), n;
  while (i < mes.size()) {
    k = mes.find_first_of(" ", i);
    if (k == string::npos) k = mes.size();
    n = k - i;
    if ((L + n >= w) && (L > 0)) { text += "\n" + string(p2, ' '); L = p2; }
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
    if (token.find("--") != 0) MstUtils::error("could not parse options, could not understand the token '" + token + "'", "MstOptions::setOptions");
    token = token.substr(2);
    if (i < argc - 1) {
      givenOptions[token] = (string) argv[i+1];
      i += 2;
    } else {
      givenOptions[token] = "";
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

int MstOptions::getInt(const string& opt, int defVal) {
  if (givenOptions.find(opt) == givenOptions.end()) return defVal;
  return MstUtils::toInt(givenOptions[opt]);
}

real MstOptions::getReal(const string& opt, real defVal) {
  if (givenOptions.find(opt) == givenOptions.end()) return defVal;
  return MstUtils::toReal(givenOptions[opt]);
}

bool MstOptions::isInt(const string& opt) {
  if (givenOptions.find(opt) == givenOptions.end()) return false;
  return MstUtils::isInt(givenOptions[opt]);
}

bool MstOptions::isReal(const string& opt) {
  if (givenOptions.find(opt) == givenOptions.end()) return false;
  return MstUtils::isReal(givenOptions[opt]);
}

string MstOptions::getString(const string& opt, string defVal) {
  if (givenOptions.find(opt) == givenOptions.end()) return defVal;
  return givenOptions[opt];
}

bool MstOptions::isGiven(const string& opt) {
  if (givenOptions.find(opt) == givenOptions.end()) return false;
  return true;
}
