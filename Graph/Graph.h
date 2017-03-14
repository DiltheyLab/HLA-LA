/*
 * Graph.h
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef GRAPH_H_
#define GRAPH_H_

class Graph;
class Node;
class Edge;

#include <vector>
#include <string>
#include <set>
#include <utility>

#include "Node.h"
#include "Edge.h"
#include "LocusCodeAllocation.h"
#include "HaplotypePanel.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/map.hpp>

using namespace std;

struct diploidEdgePointerPath
{
	vector<Edge*> h1;
	vector<Edge*> h2;
};

struct levelInfo
{
	int nodes;
	int edges;
	int symbols;
	int symbols_CODE;
};

extern std::string separatorForSerialization;


// Forward declaration of class boost::serialization::access
namespace boost {
namespace serialization {
class access;
}
}

class Graph {
private:
    friend class boost::serialization::access;

    template<typename Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & Nodes;
		ar & Edges;
		ar & NodesPerLevel;
		ar & filename_last_read;
		ar & haveComputedGapEdgePaths;
		ar & completedGapEdgePaths;
		ar & pseudoNodes;
		ar & pseudoEdges_correspondingToGapEdgePaths;
		ar & gapEdgePaths_connectedNodes_forwards;
		ar & gapEdgePaths_connectedNodes_backwards;
	}

public:
	Graph();
	virtual ~Graph();

	set<Node*> Nodes;
	set<Edge*> Edges;
	vector< set<Node*> > NodesPerLevel;

	void registerNode(Node* n, unsigned int level);
	void registerEdge(Edge* e);
	void unRegisterNode(Node* n);
	void unRegisterEdge(Edge* e);

	void checkConsistency(bool terminalCheck);
	void checkLocusOrderConsistency(vector<string> loci);
	vector<string> getAssignedLoci();
	
	void freeMemory();

	void writeToFile(string filename);
	void readFromFile(string filename);
	void printComplexity (string filename);


	int trimGraph(bool remove2DHLA = false);

	std::vector<std::string> simulateHaplotypes(int number, bool includeGaps = false);
	void makeEdgesGaps(double proportion);

	void removeStarPaths();


	void graphViz(int level_start, int level_stop, std::string output_filename);
	void graphViz2(std::string locus_string, std::string output_filename);

	vector<levelInfo> getLevelInfo();
	
	std::string getOneLocusIDforLevel(unsigned int level);

	diploidEdgePointerPath simulateRandomDiploidPath();

	std::string filename_last_read;

	void buildFromHaplotypes(HaplotypePanel hp, bool wantPGFprotection, int want_suffix_length);

	std::set<Edge*> getEdgesEmanatingFromLevel(unsigned int lI);
	
	void checkStructure();
	void regenerateNodeIncomingOutgoingEdges();

	std::set<Node*> pseudoNodes;
	std::vector<std::vector<Edge*>> completedGapEdgePaths;
	std::map<Edge*, unsigned int> pseudoEdges_correspondingToGapEdgePaths;
	std::map<Node*, std::map<Node*, Edge*>> gapEdgePaths_connectedNodes_forwards;
	std::map<Node*, std::map<Node*, Edge*>> gapEdgePaths_connectedNodes_backwards;

	void computeGapEdgePaths();

	bool get_haveComputedGapEdgePaths() {return haveComputedGapEdgePaths;};

	static std::vector<std::string> readGraphLoci(std::string graphDir);

	void checkAlignmentBackbonePresences(std::vector<std::pair<std::string, std::vector<int>>> backbones);
	std::set<unsigned int> checkPartialAlignmentBackbonePresence(std::pair<std::string, std::vector<int>> alignmentBackbone, bool greedy, bool verbose = false);
	std::set<unsigned int> checkPartialSequencePresence(std::string S, bool greedy, bool verbose = false);
	bool checkSequencePresence(std::string S, int startingAt, bool verbose = false);
	bool checkSequencePresence(std::string S, bool verbose = false);

	void checkAlignmentBackbonePresences_ignoreGraphGaps(std::vector<std::pair<std::string, std::vector<int>>> backbones);
	std::set<unsigned int> checkPartialAlignmentBackbonePresence_ignoreGraphGaps(std::pair<std::string, std::vector<int>> alignmentBackbone, bool greedy, bool verbose = false);
	std::set<unsigned int> checkPartialSequencePresence_ignoreGraphGaps(std::string S, bool greedy, bool verbose = false);
	bool checkSequencePresence_ignoreGraphGaps(std::string S, int startingAt, bool verbose = false);

protected:
	LocusCodeAllocation constructCODEfromGraph();
	bool haveComputedGapEdgePaths;

};

#endif /* GRAPH_H_ */
