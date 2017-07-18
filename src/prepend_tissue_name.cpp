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
  string newName = argv[2];


  /**
   * Filtered files
   */

  ifstream nodesFile;
  ofstream csvOutput;
  csvOutput.open("../data/output/camacho-ficheiros/"+newName+".txt");
  nodesFile.open("../data/output/camacho-ficheiros/"+tissueName+".txt");

  string line = "";

  while (getline(nodesFile, line)){
    vector <string> values = splitTSV(line);

    string& mitGene = values[0];
    string& pathwayName = values[1];
    string& conGenes = values[2];

    csvOutput << newName << "\t" << mitGene << "\t" << pathwayName << "\t" << conGenes << endl;
  }


  cout << "\nAll done" << endl;

}
