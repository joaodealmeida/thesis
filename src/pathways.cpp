#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
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

string getPathwayListForGene(const map<string, vector <string>>& pathwayToGene, const string& geneName) {
  string pathwayList = "";
  for(auto& pathways : pathwayToGene){
    auto& pathwayName = pathways.first;
    auto& geneSet = pathways.second;

    if(isPresent(geneSet, geneName))
      pathwayList+=pathwayName+",";
  }

  if(pathwayList.length() != 1)
    pathwayList.pop_back();

  //pathwayList +="]";

  return pathwayList;
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

  networkEnrich << "MitoGene\tProteinGene\tCorrelation\tPathwayList\tPathwayList\tisMitocondrial\tisMitocondrial\n";


  /*
  * Read mitocondrial genes list
  */
  ifstream mitocondrialGenesFileIn;
  mitocondrialGenesFileIn.open("../data/input/mitocondrialGenes.list");

  string tempGene;
  vector<string> mitocondrialGenes;

  vector<string> visitedGenes;

  while (getline(mitocondrialGenesFileIn, tempGene))
    mitocondrialGenes.push_back(tempGene);

  mitocondrialGenesFileIn.close();

  if(networkFile.is_open()){
    string line;
    while (getline(networkFile, line)){
      vector<string> lineTk = splitTSV(line);
      const string& mitoGene = lineTk[0];
      const string& proteinGene = lineTk[1];
      const string& correlation = lineTk[2];
      string pathway1 ="";
      if( !isPresent(visitedGenes, mitoGene)){
        pathway1 = getPathwayListForGene(pathwayToGene, mitoGene);
        visitedGenes.push_back(mitoGene);
      }
      const string& pathway2 = getPathwayListForGene(pathwayToGene, proteinGene);

      string isMito = "";

      if(isPresent(mitocondrialGenes, proteinGene))
        isMito = "true";
      else
        isMito = "false";
      /*
      * GENE1 GENE2 CORRELATION [PATHWAYS1] [PATHWAY2]
      */

      ostringstream outputStream;
      outputStream << mitoGene << "\t" << proteinGene << "\t" << correlation << "\t" << pathway1 << "\t" << pathway2  << "\t" << "true" << "\t" << isMito << endl;
      networkEnrich << outputStream.str();
    }
  }

  networkEnrich.close();
  networkFile.close();

  cout << "All done" << endl;

}
