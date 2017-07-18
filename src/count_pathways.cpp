#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include "utilities.h"



int main(int argc, char* argv[]) {
  if (argc < 2) {
		cerr << "Error: Too few arguments. example ./filter_data \"Bladder\" 1 " << endl;
		return 1;
	}

  /**
   *
   *  Map all the existing patways in a tissue
   *
   */

  map <string, set<string>> pathwaysTissue;
  set <string> tissueNames;

  cout << "Starting" << endl;
  string tissueName = argv[1];


  /**
   * Filtered files
   */

  ifstream nodesFile;

  nodesFile.open(tissueName);

  //printf("\r%s", "Loading network 2...");
  //fflush(stdout);
  string line = "";
  string tissue1 = "";

  while (getline(nodesFile, line)){
    vector <string> values = splitTSV(line);

    string& tissueName = values[0];
    tissueNames.insert(tissueName);
    tissue1 = values[0];
    //string& mitG = values[1];
    //cout << "Tissue " << tissueName << endl;
    string& pathwayName = values[2];

    if(pathwayName != "NONE"){
      pathwaysTissue[tissueName].insert(pathwayName);
    }
  }

  /**
   *
   *  Write header with tissueNames
   *
   */

  /**
  *       Tecido1 | Tecido2  | ... | Tecido N |
  * Tecido1   -       2         ...     0
  * Tecido2   2   |   -         ...     4
  */

  ofstream csvOutput;
  ofstream occurrences;
  ofstream analysis;
  csvOutput.open("../data/output/camacho-ficheiros/heatmap.csv");
  occurrences.open("../data/output/camacho-ficheiros/occurrences.csv");
  analysis.open("../data/output/camacho-ficheiros/analysis.csv");

  string header =",";
  for (const auto& tecido : tissueNames) {
    header+= tecido + ",";
  }

  if(header.length() != 1)
    header.pop_back();

  analysis << header << endl;
  occurrences << header << endl;

  /**
   *
   *  Get all pathways in this file
   *
   */

  set <string> allPathways;
  header = "";
  for (const auto& tecido : pathwaysTissue) {
    const auto& pathways = tecido.second;
    for (const auto& pathway : pathways) {
        allPathways.insert(pathway);
    }
  }

  for (auto& path : allPathways) {
    header += path + ",";
  }

  if(header.length() != 1)
    header.pop_back();

  csvOutput << header << endl;


  for (const auto& tecido : tissueNames) {
    string linha = tecido+",";
    string linha2 = tecido+",";
    string linha3 = tecido+",";

    /**
     *
     *  Check if the pathway is on a tissue
     */
    for (auto& name : allPathways) {
      if (pathwaysTissue[tecido].find(name) != pathwaysTissue[tecido].end()){
        linha +="1,";
      }
      else
      {
        linha += "0,";
      }

    }


    /**
     *
     *  Match the common pathways and count how many in common
     */
    for (const auto& outroTecido : tissueNames) {
      if(tecido == outroTecido){
        linha3 += "NA,";
        linha2 += "0,";
        continue;
      }

      set <string> networkIntersection;
      std::set_intersection(pathwaysTissue[tecido].begin(),pathwaysTissue[tecido].end(),pathwaysTissue[outroTecido].begin(),pathwaysTissue[outroTecido].end(),std::inserter(networkIntersection, networkIntersection.begin()));
      string pathwaysLista = "";
      for (auto& pathway : networkIntersection) {
        pathwaysLista += pathway + " ";
      }
      if(pathwaysLista.length() != 1)
        pathwaysLista.pop_back();
      linha3 += /*to_string(networkIntersection.size())\*/ pathwaysLista + ",";
      linha2 += to_string(networkIntersection.size()) + ",";
    }

    linha.pop_back();
    linha2.pop_back();
    linha3.pop_back();
    csvOutput << linha << endl;
    occurrences << linha2 << endl;
    analysis << linha3 << endl;
  }


  cout << "\nAll done" << endl;

}
