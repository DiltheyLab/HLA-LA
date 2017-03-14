/*
 * GraphAndEdgeIndex.h
 *
 *  Created on: 29.07.2013
 *      Author: AlexanderDilthey
 */

#ifndef GRAPHANDEDGEINDEX_H_
#define GRAPHANDEDGEINDEX_H_

#include <vector>
#include <map>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "../Graph/Graph.h"

namespace boost {
namespace serialization {
class access;
}
}

class kMerInfo {
public:
	vector<unsigned char> kMer_coded;
	vector<string> kMer_deCoded;
	double p;
	vector<Edge*> traverseEdges;
	string traverseEdges_string;
	bool gapEdge;
	bool allPGF;

private:
    friend class boost::serialization::access;

    template<typename Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & kMer_coded;
		ar & kMer_deCoded;
		ar & p;
		ar & traverseEdges;
		ar & traverseEdges_string;
		ar & gapEdge;
		ar & allPGF;
	}
};

class kMerAtNode {
public:
	vector<unsigned char> km1Mer_coded;
	vector<string> km1Mer_deCoded;
	vector<Edge*> traverseEdges;
	string traverseEdges_string;
	Node* lastNewNode;
	bool allPGF;

private:
    friend class boost::serialization::access;

    template<typename Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & km1Mer_coded;
		ar & km1Mer_deCoded;
		ar & traverseEdges;
		ar & traverseEdges_string;
		ar & lastNewNode;
		ar & allPGF;
	}
};

class kMerInGraphSpec {
public:
	std::vector<Edge*> traversedEdges;
	std::string getID()
	{
		std::string forReturn;
		std::ostringstream idString;
		idString << setw(15) << (void const *) traversedEdges.front();
		idString << " -- ";
		idString << setw(15) << (void const *) traversedEdges.back();
		return idString.str();
	}

private:
    friend class boost::serialization::access;

    template<typename Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & traversedEdges;
	}
};

class kMerEdgeChain {
public:
	int sequence_begin;
	int sequence_end;

	std::vector<Edge*> traversedEdges;
	kMerEdgeChain() : sequence_begin(-1), sequence_end(-1)
	{

	}

private:
    friend class boost::serialization::access;

    template<typename Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & sequence_begin;
		ar & sequence_end;
		ar & traversedEdges;
	}
};

class GraphAndEdgeIndex {
protected:
	int kMerSize;
	Graph* g;

	std::map<std::string, std::vector<kMerInGraphSpec> > kMers;
	std::map<Node*, std::vector<Edge*>> nodes_jumpOverGaps;

	std::set<Node*> generated_nodes;
	std::set<Edge*> generated_edges;


public:
	GraphAndEdgeIndex(Graph* graph, int k);
	~GraphAndEdgeIndex();
	void Index();
	void fillEdgeJumper();
	void printIndex();
	std::vector<kMerEdgeChain*> findChains(std::string sequence);
	std::vector<kMerInGraphSpec> queryIndex(std::string kMer);
	std::vector<std::string> getIndexedkMers();

private:
    friend class boost::serialization::access;

    template<typename Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
//		ar & kMerSize;
//		ar & g;
//		ar & kMers;
//		ar & nodes_jumpOverGaps;
//		ar & generated_nodes;
//		ar & generated_edges;
	}
};

#endif /* GRAPHANDEDGEINDEX_H_ */
