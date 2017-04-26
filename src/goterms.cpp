#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <algorithm>
#include "utilities.h"

map<string, vector<string> > parsePathway(const string& directory){
  /*
  * Map genesSet to pathwayList
  */
  map<string, vector<string>> pathwayToGene;
  ifstream pathwayList;
  pathwayList.open(directory);

  if(pathwayList.is_open()){
    string line;
    while (getline(pathwayList, line)){
      vector<string> lineTokens = splitTSV(line);
      string pathwayName = lineTokens[0];
      //Discard Line Token number 2 because it has an URL
      lineTokens.erase(lineTokens.begin(),lineTokens.begin()+1);
      pathwayToGene[pathwayName] = lineTokens;
    }
  }
  pathwayList.close();

  return pathwayToGene;
}

int getPathwayListForGene(const map<string, vector <string>>& pathwayToGene, const string& geneName, const string& geneName2) {
  vector<string> pathwaysGene1;
  vector<string> pathwaysGene2;

  vector<string> commonPathways;

  for(auto& pathways : pathwayToGene){
    auto& pathwayName = pathways.first;
    auto& geneSet = pathways.second;

    if(isPresent(geneSet, geneName)){
      pathwaysGene1.push_back(pathwayName);
    }

    if(isPresent(geneSet, geneName2)){
      pathwaysGene2.push_back(pathwayName);
    }
  }

  auto it = set_intersection(pathwaysGene1.begin(), pathwaysGene1.end(), pathwaysGene2.begin(), pathwaysGene2.end(), back_inserter(commonPathways));

  return commonPathways.size();
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cerr << "Error: Too few arguments" << endl;
		return 1;
	}

  /*
  * argv[0] - pathway file directory | argv[1] - network file directory | argv[2] - enrich network output directory
  */

  //Parse pathway file
  const map<string, vector <string>>& pathwayToGene = parsePathway(argv[1]);

  cout << "Starting" << endl;
  cout << argv[1] << endl;
  cout << argv[2] << endl;
  cout << argv[3] << endl;

  //Enrich network Table
  ifstream networkFile;
  networkFile.open(argv[2]);

  ofstream networkEnrich;
  networkEnrich.open(argv[3]);

  if (!networkEnrich.is_open())
    cerr << endl << "ERROR: Could not open output file." << endl;

  string tempGene;

  networkEnrich << "P1" << "\t" << "P2" << "\t" << "PairScore" << "\t" << "GEDScore" << "\t" << "GrletSigScore" << "\t" << "Info" << "\t" <<  "GOSum" << "\t" << "gedevoGOPercScore" << endl;
  if(networkFile.is_open()){
    string line;

    //Waste first 9 lines
    for (int i = 0; i < 9; i++) {
      //cout << "hit" << endl;
      getline(networkFile, line);
    }
    //P1	P2	PairScore	GEDScore	GrletSigScore	Info
    while (getline(networkFile, line)){
      vector<string> lineTk = splitTSV(line);
      const string& P1 = lineTk[0];
      const string& P2 = lineTk[1];
      const string& PairScore = lineTk[2];
      const string& GEDScore = lineTk[3];
      const string& GrletSigScore = lineTk[4];
      const string& Info = lineTk[5];

      //cout << line << endl;
      const int& similarity = getPathwayListForGene(pathwayToGene, P1, P2);

      ostringstream outputStream;
      outputStream << P1 << "\t" << P2 << "\t" << PairScore << "\t" << GEDScore << "\t" << GrletSigScore << "\t" << Info << "\t" <<  similarity << "\t" << similarity << endl;
      networkEnrich << outputStream.str();
    }
  }

  networkEnrich.close();
  networkFile.close();

  cout << "All done" << endl;

}
