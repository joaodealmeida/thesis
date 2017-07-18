#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "utilities.h"

using namespace std;

struct tissueFile_t {
	vector<string> genes;
	vector<vector<double>> values;
};

vector<vector<double>> eliminateOutlier(vector<double> x, vector<double> y){
	assert(x.size() == y.size());

	double n = x.size();
	double sumX = 0, sumY = 0;
	double meanX = 0;
	double meanY = 0;

	double sdX = 0, sdY = 0;


	for (unsigned int i = 0; i < n; i++) {
		sumX += x[i];
		sumY += y[i];
	}

	meanX = sumX / x.size();
	meanY = sumY / y.size();

	for (unsigned int i = 0; i < n; i++) {
		//Standard Deviation
		sdX += pow(x[i] - meanX, 2);
		sdY += pow(y[i] - meanY, 2);
	}

	sdX = sqrt(sdX / x.size());
	sdY = sqrt(sdY / y.size());


	//Check if any values is bigger that 4 comparing to his sd, if so eliminate it.
	for (unsigned int i = 0; i < n; i++) {


		//[ux-4SD, ux+4sd ]
		if( (meanX - 4*sdX > x[i] || meanX + 4*sdX < x[i]) || (meanY - 4*sdY > y[i] || meanY + 4*sdY < y[i])){
			x.erase(x.begin() + i);
			y.erase(y.begin() + i);
			n = x.size();
			i--;
		}
	}

	vector<vector<double>> ret;
	ret.push_back(x);
	ret.push_back(y);

	return ret;
}

double calcPearson(const vector<double>& x, const vector<double>& y) {
	assert(x.size() == y.size());


	double n = x.size();

	double sumXY = 0;
	double sumX = 0, sumY = 0;
	double sumXX = 0, sumYY = 0;



	for (unsigned int i = 0; i < n; i++) {
		//numerator
		sumXY += x[i] * y[i];

		sumX += x[i];
		sumY += y[i];

		// denominator
		sumXX += x[i] * x[i];
		sumYY += y[i] * y[i];
	}

	/*
	* Taken from:
	* 	https://upload.wikimedia.org/math/b/c/7/bc7fa889f31ddbec0a05df656de340a5.png
	*/

	double numerator = n * sumXY - sumX * sumY;
	double sdX = sqrt(n * sumXX - sumX * sumX);
	double sdY = sqrt(n * sumYY - sumY * sumY);
	double denominator = sdX * sdY;

	double correlation = denominator == 0 ? 0 : numerator / denominator;

	return correlation;
}

/**
*
* Compile:
*  g++ correlations.cpp -o correlations.out -Wall -std=c++11
*
* Usage:
*  ./correlations.out <tissue name> <mitocondrial genes list>
*
* Example:
*  ./correlations.out Bladder ../data/input/mitocondrialGenes.list
*
*/
int main(int argc, char* argv[]) {
	if (argc < 3) {
		cerr << "Error: Too few arguments" << endl;
		return 1;
	}


	/*
	* Read mitocondrial genes list
	*/
	ifstream mitocondrialGenesFileIn;
	mitocondrialGenesFileIn.open(argv[2]);

	string tempGene;
	vector<string> mitocondrialGenes;

	while (getline(mitocondrialGenesFileIn, tempGene))
		mitocondrialGenes.push_back(tempGene);

	mitocondrialGenesFileIn.close();


	/*
	* Read tissue file
	*/
	string tissueName = argv[1];

	ifstream tissueFileIn;
	tissueFileIn.open("../data/output/tissues-output/" + tissueName + ".txt");

	cout << tissueName << " ... " << "Loading" << flush;

	if (tissueFileIn.is_open()) {
		string line;

		// discarding useless info
		getline(tissueFileIn, line);
		getline(tissueFileIn, line);
		getline(tissueFileIn, line);

		tissueFile_t tissueFile;

		while (getline(tissueFileIn, line)) {
			vector<string> lineTokens = splitTSV(line);

			// save gene name
			tissueFile.genes.push_back(lineTokens[0]);

			// save gene data
			vector<double> valuesRow;

			for (unsigned int i = 1; i < lineTokens.size(); i++)
				valuesRow.push_back(stod(lineTokens[i], nullptr));

			tissueFile.values.push_back(valuesRow);
		}


		/*
		* Map tissue gene name to data row
		*/
		map<string, int> geneToRow;
		for (unsigned int i = 0; i < tissueFile.genes.size(); i++)
			geneToRow[tissueFile.genes[i]] = i;

		/*
		* Map mitocondrialGene to overlayGenes
		*/
		map<string, vector<string>> mitocondrialGeneOverlay;

		/*
		* Read overlay genes file
		*/

		ifstream overlayGenesFileIn;
		overlayGenesFileIn.open("../data/output/overlays/overlayGenesList.txt");

		cout << "Overlay Genes file" << " ... " << "Loading" << flush;

		if (overlayGenesFileIn.is_open()) {
			string line;

			while (getline(overlayGenesFileIn, line)){
				vector<string> overlayPair = splitOverlay(line);

				mitocondrialGeneOverlay[overlayPair[0]].push_back(overlayPair[1]);
			}

			cout << "\nOverlay Genes - All done! -" << endl;
		}

		/*
		* Read genome info file
		*/

		vector<string> proteinGenes;
		ifstream proteinGenesFileIn;
		proteinGenesFileIn.open("../data/input/proteinEncodingGenes.txt");

		cout << "Protein Genes file" << " ... " << "Loading" << flush;

		if (proteinGenesFileIn.is_open()) {
			string line;
			while (getline(proteinGenesFileIn, line)){
				proteinGenes.push_back(line);
			}

			cout << "\nProtein Genes - All done! - " << proteinGenes.size() << endl;
		}



		/*
		* Calculate correlations
		*/
		ofstream allCorrelationsOut, posCorrelationsOut, negCorrelationsOut, overlayCorrelationsOut;
		allCorrelationsOut.open("../data/output/correlations-output/" + tissueName + "-all.txt");
		posCorrelationsOut.open("../data/output/correlations-output/" + tissueName + "-pos.txt");
		negCorrelationsOut.open("../data/output/correlations-output/" + tissueName + "-neg.txt");
		overlayCorrelationsOut.open("../data/output/correlations-output/" + tissueName + "-overlay.txt");
		if (!allCorrelationsOut.is_open()) {
			cerr << endl << "ERROR: Could not open 'all' output file for: " << tissueName << endl;
		} else if (!posCorrelationsOut.is_open()) {
			cerr << endl << "ERROR: Could not open 'pos' output file for: " << tissueName << endl;
		} else if (!negCorrelationsOut.is_open()) {
			cerr << endl << "ERROR: Could not open 'neg' output file for: " << tissueName << endl;
		} else if (!overlayCorrelationsOut.is_open()) {
			cerr << endl << "ERROR: Could not open 'overlay' output file for: " << tissueName << endl;
		}
			else {
				int progress = 0;

			for (const auto& mitocondrialGene : mitocondrialGenes) {
				for (const auto& tissueGene : tissueFile.genes) {
					// skip self-correlation
					if (mitocondrialGene == tissueGene)
						continue;

					vector<vector<double>> out = eliminateOutlier(tissueFile.values[geneToRow[mitocondrialGene]], tissueFile.values[geneToRow[tissueGene]]);
					//double pearson = calcPearson(tissueFile.values[geneToRow[mitocondrialGene]], tissueFile.values[geneToRow[tissueGene]]);
					double pearson = calcPearson(out[0], out[1]);

					ostringstream outputStream;
					outputStream << mitocondrialGene << "\t" << tissueGene << "\t" << pearson << endl;

					allCorrelationsOut << outputStream.str();

					if (pearson > 0.7){
						if(!isPresent(mitocondrialGeneOverlay[mitocondrialGene], tissueGene) && isPresent(proteinGenes, tissueGene)){
							overlayCorrelationsOut << outputStream.str();
						}
						else
							posCorrelationsOut << outputStream.str();

					}
					else if (pearson < -0.7)
						negCorrelationsOut << outputStream.str();




				}
				printf("\r%s ... %5.1f %%", tissueName.c_str(), (++progress) * 100.0 / mitocondrialGenes.size());
				fflush(stdout);
			}

			cout << "\r" << tissueName << " ... OK!    " << endl;

			allCorrelationsOut.close();
			posCorrelationsOut.close();
			negCorrelationsOut.close();
			overlayCorrelationsOut.close();
		}

		tissueFileIn.close();
	} else {
		cerr << endl << "ERROR: Could not open file of tissue: " << tissueName << endl;
	}

	cout << "- All done! -" << endl;

	return 0;
}
