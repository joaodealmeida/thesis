#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include "utilities.h"




struct Node {
  string name;
  std::set <string> pathways;
  bool isMitocondrial;
  double betweenessCentrality;
  double closenessCentrality;
  double degree;
  double numberOfDirectedEdges;
  double radiality;
  double topologicalCoefficient;

  bool operator==(const Node& n) const
  {
      return (name == n.name);
  }

  bool operator<(const Node& n) const
  {
    return (name < n.name);
  }

};

struct Edge
{
  double correlation;
  string mitocondrialGene;
  string proteinGene;

};

struct equal_node
{
  inline bool operator() (const Node& node1, const Node& node2)
  {
    return (node1.name != node2.name);
  }
};

struct order_by_name
{
  inline bool operator() (const Node& node1, const Node& node2)
  {
    return (node1.name < node2.name);
  }
};


struct order_by_betweeness
{
  inline bool operator() (const Node& node1, const Node& node2)
  {
    return (node1.betweenessCentrality > node2.betweenessCentrality);
  }
};

struct order_by_closeness
{
  inline bool operator() (const Node& node1, const Node& node2)
  {
    return (node1.closenessCentrality > node2.closenessCentrality);
  }
};

struct order_by_degree
{
  inline bool operator() (const Node& node1, const Node& node2)
  {
    return (node1.degree > node2.degree);
  }
};

struct order_by_numberofedges
{
  inline bool operator() (const Node& node1, const Node& node2)
  {
    return (node1.numberOfDirectedEdges > node2.numberOfDirectedEdges);
  }
};

struct order_by_radiality
{
  inline bool operator() (const Node& node1, const Node& node2)
  {
    return (node1.radiality > node2.radiality);
  }
};

struct order_by_topological
{
  inline bool operator() (const Node& node1, const Node& node2)
  {
    return (node1.topologicalCoefficient > node2.topologicalCoefficient);
  }
};

//Example how to sort by operator
//std::sort(vec.begin(), vec.end(), less_than_key());


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

/**
 * EXAMPLE: ./filter_data "Bladder" 1
 * 1  - betweenessCentrality
 * 2  - closenessCentrality
 * 3  - degree
 * 4  - numberOfDirectedEdges
 * 5  - radiality
 * 6  - topologicalCoefficient
 */

int main(int argc, char* argv[]) {
  if (argc < 2) {
		cerr << "Error: Too few arguments. example ./filter_data \"Bladder\" 1 " << endl;
		return 1;
	}

  /*
  * argv[1] - tissueName
  */

  cout << "Starting" << endl;
  cout << argv[1] << endl;
  string tissueName = argv[1];




  /**
   * Map pathway name - nodes
   */

  map <string, vector<string>> pathwayNodes;

  /**
   * Parse nodes file
   */
  vector <Node> networkNodes;

  ifstream nodesFile;
  string line;
  nodesFile.open("../data/output/cytoscape-networks/"+ tissueName + "-nodes.csv");

  /**
   * Discard first row
   */

   printf("\r%s", "Loading nodes...");
   fflush(stdout);
  getline(nodesFile, line);

  while (getline(nodesFile, line)){
    vector <string> values = splitCyto(line);

    double betweenessCentrality = stod(ReplaceAll(values[2], "\"" ,""));
    double closenessCentrality = stod(ReplaceAll(values[3], "\"" ,""));
    double degree = stod(ReplaceAll(values[5], "\"" ,""));
    bool isMitocondrial = to_bool(ReplaceAll(values[7], "\"" ,""));
    string name = ReplaceAll(values[9], "\"" ,"");
    double numberOfDirectedEdges = stod(ReplaceAll(values[11], "\"" ,""));

    set <string> pathwayList;
    for (auto& pathway : splitPathway(ReplaceAll(values[14],"\"","")) ){
      pathwayList.insert(pathway);
      pathwayNodes[pathway].push_back(name);
    }

    double radiality = stod(ReplaceAll(values[15], "\"" ,""));
    double topologicalCoefficient = stod(ReplaceAll(values[20], "\"" ,""));

    Node n;
    n.betweenessCentrality = betweenessCentrality;
    n.closenessCentrality = closenessCentrality;
    n.degree = degree;
    n.isMitocondrial = isMitocondrial;
    n.name = name;
    n.numberOfDirectedEdges = numberOfDirectedEdges;
    n.pathways = pathwayList;
    n.radiality = radiality;
    n.topologicalCoefficient = topologicalCoefficient;
    networkNodes.push_back(n);
  }
  printf("\r%s", "Loading edges...");
  fflush(stdout);
  /**
   * Parse edges file
   */
  vector <Edge> networkEdges;
  map <string, vector <Node>> nodeInteractions;


  ifstream edgesFile;

  edgesFile.open("../data/output/cytoscape-networks/"+ tissueName+"-edges.csv");

  /**
   * Discard first row
   */
  getline(edgesFile, line);

  while (getline(edgesFile, line)){
    vector <string> values = splitCyto(line);

    double correlation = stod(ReplaceAll(values[1], "\"" ,""));
    vector <string> interaction = splitOverlay(ReplaceAll(values[4], "\"", ""));


    string node1 = interaction[0];
    string node2 = interaction[3];

    Edge edg;

    auto it = find_if(networkNodes.begin(), networkNodes.end(), [&node2](const Node& obj) {return obj.name == node2;});
    auto it2 = find_if(networkNodes.begin(), networkNodes.end(), [&node1](const Node& obj) {return obj.name == node1;});

    //if(it != networkNodes.end()){
      nodeInteractions[node1].push_back(*it);
      nodeInteractions[node2].push_back(*it2);
    //}

    edg.correlation = correlation;
    edg.mitocondrialGene = node1;
    edg.proteinGene = node2;
    networkEdges.push_back(edg);

  }


  /**
   * Tecido | Medida  | Genes | Genes path x-1 | Genes pathway x
   * GENE1    degree    GENE2,GENE3 GENE2       GENE2, GENE3
   */

  ofstream networkEnrich;
  networkEnrich.open("../data/output/camacho-ficheiros/"+ tissueName + ".txt");

  sort(networkNodes.begin(), networkNodes.end(), order_by_degree());

   for (auto& node : networkNodes){
     if(node.isMitocondrial){

       string neighboorGeneOnPathway = "";
       vector<Node> genesWithPathway;
       vector<Node> genesWithoutPathway;
       vector<string> neighboorGenes;

       for (auto& pathway : pathwayNodes){
         string pathwayName = pathway.first;
         auto search = node.pathways.find(pathwayName);
         genesWithPathway.push_back(node);

         if(search != node.pathways.end() ){
           /**
            *
            * Gene is on a pathway, lets find the neighboor genes that share the same pathway
            *
            */
            neighboorGeneOnPathway = "";

            std::set<Node> nodeInteractionsSet(nodeInteractions[node.name].begin(), nodeInteractions[node.name].end());

            for (auto& neighboorGene : nodeInteractionsSet){

              if(neighboorGene.pathways.find(pathwayName) != neighboorGene.pathways.end()){
                neighboorGeneOnPathway += neighboorGene.name + ",";
                genesWithPathway.push_back(neighboorGene);
              }
            }

            if(neighboorGeneOnPathway.size() > 0) {
              neighboorGeneOnPathway.pop_back();
              networkEnrich << node.name << "\t" << pathwayName << "\t" << neighboorGeneOnPathway << endl;
            }

         }

       }
       /**
        *
        *  Find genes that are neighbor but don't belong to any pathway
        *
        */
        std::transform(std::begin(nodeInteractions[node.name]),
                std::end(nodeInteractions[node.name]),
                std::back_inserter(neighboorGenes),
                [](Node d) { return d.name; }
               );

        std::set<Node> neighboorGenesNodeSet(nodeInteractions[node.name].begin(), nodeInteractions[node.name].end());
        std::set<Node> genesWithPathwayNodeSet(genesWithPathway.begin(), genesWithPathway.end());

        set_difference(neighboorGenesNodeSet.begin(), neighboorGenesNodeSet.end(), genesWithPathwayNodeSet.begin(), genesWithPathwayNodeSet.end(), back_inserter(genesWithoutPathway));

       //networkEnrich << "SIZE: " << neighboorGenesNodeSet.size() << "PathG: " << genesWithPathwayNodeSet.size() << "No PathG: " << genesWithoutPathway.size() << "\t";
       neighboorGeneOnPathway = "";
       for (auto& gene : genesWithoutPathway) {
         neighboorGeneOnPathway += gene.name + ",";
       }

       if(neighboorGeneOnPathway.size() > 0 ){
         neighboorGeneOnPathway.pop_back();
         networkEnrich << node.name << "\t" << "NONE" << "\t" << neighboorGeneOnPathway << endl;
       }
    }
  }

  if (!networkEnrich.is_open())
    cerr << endl << "ERROR: Could not open output file." << endl;

  networkEnrich.close();

  cout << "\nAll done" << endl;

}
