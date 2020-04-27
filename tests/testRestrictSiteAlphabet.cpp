#include "msttypes.h"
#include "dtermen.h"
#include "mstoptions.h"
#include "mstsystem.h"
#include "mstrotlib.h"
//#include <chrono>

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Tests EnergyTable::restrictSiteAlphabet(). Options:");
  op.addOption("etab", "the path to the energy table that will be modified",true);
  op.addOption("rsa", "a string with the site alphabet specified for each position. e.g., for an energy table with 3 positions, the input could be 'ALA,LEU ALA,GLY ALA,LEU,GLY'");
  op.addOption("p", "template PDB file. can be passed to the function to constrain allowed residues types at each position");
  op.addOption("s", "selection of residues in the template to be set to UNK");
  op.addOption("c", "if provided, this energy table will be formatted to act as a constraint on the non-restricted version");
  op.addOption("o", "output base.", true);
  op.setOptions(argc, argv);
  
  if ((op.isGiven("rsa")) == (op.isGiven("p"))) MstUtils::error("Either a structure or string must be provided for constraints, but not both");
  
  EnergyTable etab(op.getString("etab"));
  
  vector<vector<string>> restrictedSiteAlphabets;
  
  if (op.isGiven("p")) {
    cout << "--p provided, using the protein sequence to constrain the etab" << endl;
    Structure S(op.getString("p"));
    vector<Residue*> all_residues = S.getResidues();
    
    if (op.isGiven("s")) {
      selector sel(S);
      vector<Residue*> unk_res = sel.selectRes(op.getString("s"));
      cout << "selecting " << unk_res.size() << " residues to be set to 'UNK'" << endl;
      for (Residue* R: unk_res) R->setName("UNK");
    }
    
    vector<string> etab_sites = etab.getSites();
    
    cout << "number of sites in the energy table: " << MstUtils::toString(etab_sites.size()) << endl;
    for (string site_name : etab_sites) {
      vector<string> split = MstUtils::split(site_name,",");
      if (split.size() != 2) MstUtils::error("Site name should be a CHAINID,RESNUM");
      string chain_ID = split[0];
      int res_num = MstUtils::toInt(split[1]);
      cout << "searching for position: " << chain_ID << " " << MstUtils::toString(res_num) << " in protein residues..."<< endl;
      for (Residue* R : all_residues) {
        if ((R->getChainID() == chain_ID) && (R->getNum() == res_num)) {
          cout << "found" << endl;
          restrictedSiteAlphabets.push_back({R->getName()});
          break;
        }
      }
    }
  } else {
    //use the provided string to constrain the etab
    cout << "--rsa provided, using the string to constrain the etab" << endl;
    vector<string> site_split = MstUtils::split(op.getString("rsa")," ");
    restrictedSiteAlphabets.resize(site_split.size());
    cout << "there are " << site_split.size() << " positions specified in the selection string" << endl;
    if (site_split.size() != etab.numSites()) MstUtils::error("The number of positions in the restrictedSiteAlphabet and energy table do not match: ("+MstUtils::toString(site_split.size())+") and ("+MstUtils::toString(etab.numSites())+")");
    for (string site : site_split) {
      vector<string> aa_split = MstUtils::split(site,",");
      restrictedSiteAlphabets.push_back(aa_split);
    }
  }
  
  bool constraint_etab = op.isGiven("c");
  EnergyTable new_etab = etab.restrictSiteAlphabet(restrictedSiteAlphabets,constraint_etab);
  new_etab.writeToFile("./"+op.getString("o")+".etab");
}

