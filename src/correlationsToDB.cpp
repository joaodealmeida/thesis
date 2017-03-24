#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include <sqlite3.h>

#include "utilities.h"

using namespace std;

const unsigned int BUFFER_SIZE = 256;

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
	double denominator = sqrt(n * sumXX - sumX * sumX) * sqrt(n * sumYY - sumY * sumY);

	return denominator == 0 ? 0 : numerator / denominator;
}

void resetTableDB(const string& tableName, sqlite3* db) {
	cout << "Resetting DB table ... " << flush;

	// drop table
	string statement = "DROP TABLE IF EXISTS '" + tableName + "'";

	sqlite3_exec(db, statement.c_str(), nullptr, nullptr, nullptr);

	// create table
	statement = "CREATE TABLE '" + tableName + "' ("
	"\"id\" integer not null primary key autoincrement,"
	"\"gene1\" varchar not null,"
	"\"gene2\" varchar not null,"
	"\"correlation\" float not null)";

	sqlite3_exec(db, statement.c_str(), nullptr, nullptr, nullptr);

	cout << "OK!" << endl;
}

/**
*
* Compile:
*  g++ correlationsToDB.cpp -o correlationsToDB.out -Wall -lsqlite3 -std=c++11
*
* Usage:
*  ./correlationsToDB.out <tissue name> <mitocondrial genes list>
*
* Example:
*  ./correlationsToDB.out Bladder ../data/input/mitocondrialGenes.list
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

	cout << "Reading " << tissueName << " file ... " << flush;

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

		cout << "OK!" << endl;


		/*
		* Map tissue gene name to data row
		*/
		map<string, int> geneToRow;
		for (unsigned int i = 0; i < tissueFile.genes.size(); i++)
			geneToRow[tissueFile.genes[i]] = i;


		/*
		* Open DB
		*/
		tissueName = str_slug(tissueName);

		ostringstream ss;
		ss << "../../coexpr/database/" << tissueName << ".sqlite";

		sqlite3* db;

		cout << "Opening DB ... ";
		if (sqlite3_open(ss.str().c_str(), &db) == SQLITE_OK) {
			cout << "OK!" << endl;

			resetTableDB(tissueName, db);

			sqlite3_exec(db, "PRAGMA synchronous = OFF", nullptr, nullptr, nullptr);
			sqlite3_exec(db, "PRAGMA journal_mode = MEMORY", nullptr, nullptr, nullptr);

			ss.str(std::string());
			ss << "INSERT INTO 'correlations' VALUES (NULL, @gene1, @gene2, @correlation)";

			sqlite3_stmt* stmt;
			sqlite3_prepare_v2(db, ss.str().c_str(), BUFFER_SIZE, &stmt, nullptr);

			int progress = 0;
			sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);

			cout << tissueName << " ... " << "Loading" << flush;

			/*
			* Calculate correlations
			*/
			for (const auto& mitocondrialGene : mitocondrialGenes) {
				for (const auto& tissueGene : tissueFile.genes) {
					// skip self-correlation
					if (mitocondrialGene == tissueGene)
						continue;

					vector<vector<double>> out = eliminateOutlier(tissueFile.values[geneToRow[mitocondrialGene]], tissueFile.values[geneToRow[tissueGene]]);
					//double pearson = calcPearson(tissueFile.values[geneToRow[mitocondrialGene]], tissueFile.values[geneToRow[tissueGene]]);
					double pearson = calcPearson(out[0], out[1]);

					ss.str(std::string());
					ss << pearson;

					sqlite3_bind_text(stmt, 1, mitocondrialGene.c_str(), -1, SQLITE_TRANSIENT);
					sqlite3_bind_text(stmt, 2, tissueGene.c_str(), -1, SQLITE_TRANSIENT);
					sqlite3_bind_text(stmt, 3, ss.str().c_str(), -1, SQLITE_TRANSIENT);

					sqlite3_step(stmt);

					sqlite3_clear_bindings(stmt);
					sqlite3_reset(stmt);
				}

				printf("\r%s ... %5.1f %%", tissueName.c_str(), (++progress) * 100.0 / mitocondrialGenes.size());
				fflush(stdout);
			}

			sqlite3_exec(db, "END TRANSACTION", nullptr, nullptr, nullptr);

			cout << "\r" << tissueName << " ... OK!    " << endl;

			cout << "Creating DB index (hang on just a bit more!) ... " << flush;

			ss.str(std::string());
			ss << "CREATE INDEX 'correlations_index' ON 'correlations' ('correlation')";

			sqlite3_exec(db, ss.str().c_str(), nullptr, nullptr, nullptr);

			cout << "OK!" << endl;

			sqlite3_finalize(stmt);
			sqlite3_close(db);
		} else {
			cerr << endl << "ERROR: Failed to open DB." << endl;
		}

		tissueFileIn.close();
	} else {
		cerr << endl << "ERROR: Could not open file of tissue: " << tissueName << endl;
	}

	cout << "- All done! -" << endl;

	return 0;
}
