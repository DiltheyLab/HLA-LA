/*
 * VirtualNW.h
 *
 *  Created on: 30.06.2013
 *      Author: AlexanderDilthey
 */

#ifndef VIRTUALNWUNIQUE_H_
#define VIRTUALNWUNIQUE_H_

#include <set>
#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "../../Graph/Graph.h"
#include "alignmentContext.h"
#include "../reads/verboseSeedChain.h"

#include <assert.h>

namespace mapper {
namespace aligner {

class localExtension_pathDescription {
public:
	std::vector<std::vector<int> > coordinates;
	std::vector<Edge*> usedEdges;
	double Score;
	std::string alignedSequence;
	std::string alignedGraph;
	std::vector<int> alignedGraph_levels;

	void _printExtension()
	{
		std::cout << "localExtension_pathDescription object:\n" << std::flush;
		for(unsigned int cI = 0; cI < coordinates.size(); cI++)
		{
			for(unsigned int cII = 0; cII < coordinates.at(cI).size(); cII++)
			{
				std::cout << coordinates.at(cI).at(cII) << " " << std::flush;
			}
			std::cout << "\n" << std::flush;
			if(cI < usedEdges.size())
			{
				std::cout << "\t " << usedEdges.at(cI) << "\t" << std::flush;
				if(usedEdges.at(cI) != 0)
				{
					Edge* e = usedEdges.at(cI);
					Graph* g = usedEdges.at(cI)->From->g;
					assert(g != 0);
					std::cout << e->getEmission() << " " << std::flush;
				}
				std::cout << "\n" << std::flush;
			}
		}
	}

	mapper::reads::verboseSeedChain toVerboseSeedChain();


};

class NWPath;
class NWEdge;
class NWPath;

class VirtualNWTable_Unique {
public:
	std::set<NWPath*> paths;

	std::map<int, std::set<NWEdge*>> index_edges_start_x;
	std::map<int, std::set<NWEdge*>> index_edges_stop_x;
	std::map<int, std::set<NWEdge*>> index_edges_start_y;
	std::map<int, std::set<NWEdge*>> index_edges_stop_y;

	std::map<Edge*, std::set<NWEdge*> > index_graphEdges_2_NWEdges;

	std::map<int, std::map<int, std::map<int, std::set<NWEdge*> > > > index_edges_from;
	std::map<int, std::map<int, std::map<int, std::set<NWEdge*> > > > index_edges_to;

	alignmentContext* Context;

	// for debugging
	std::map<NWPath*, int> _ordered_paths;
	std::map<NWPath*, std::map<NWEdge*, int> > _ordered_paths_edges;

	VirtualNWTable_Unique(alignmentContext* alignmentContext_) : Context(alignmentContext_)
	{

	}

	void addPath(NWPath* p);
	void removePath(NWPath* p);

	void checkConsistency();
	void freeMemory();

	void _testTracePath(NWPath* p, std::string& reconstructed_sequence, std::string& reconstructed_graph, std::vector<int>& reconstructed_graph_levels);

	bool hasEdgeEmanatingFrom(int x, int y, int z);
	bool hasEdgeGoingInto(int x, int y, int z);

	bool edgePresentInPath(int from_x, int from_y, int from_z, int to_x, int to_y, int to_z, Edge* usedEdge);

	std::set<NWEdge*> getEdgesEmanatingFrom(int x, int y, int z);
	std::set<NWEdge*> getEdgesGoingInto(int x, int y, int z);

	unsigned int getNumPaths()
	{
		return paths.size();
	}

	unsigned int getNumEdges();

	unsigned int getNumEntryEgdes();
	unsigned int getNumExitEdges();

	std::vector<NWEdge*> getEntryEdges();
	std::vector<NWEdge*> getExitEdges();

	std::set<NWEdge*> getAllEdges();

	void print();
	void _enumeratePathsAndEdges();

	void printSequenceCoverageStats();
};

class NWPath {
public:
	NWPath();

	std::map<int, std::map<int, std::map<int, std::set<NWEdge*> > > > edges_from;
	std::map<int, std::map<int, std::map<int, std::set<NWEdge*> > > > edges_to;

	std::set<NWEdge*> first_edges;
	std::set<NWEdge*> last_edges;
	std::set<NWEdge*> all_edges;
	std::set<NWEdge*> entry_edges;
	std::set<NWEdge*> exit_edges;

	VirtualNWTable_Unique* table;

	void checkConsistency();
	void freeMemory();

	NWPath* clonePathWithoutTable();

	std::set<NWEdge*> forwardEdges(int x, int y, int z);
	std::set<NWEdge*> forwardEdges(NWEdge* e);

	std::set<NWEdge*> backwardEdges(int x, int y, int z);
	std::set<NWEdge*> backwardEdges(NWEdge* e);

	void verify_edges();
	NWEdge* createAndAddEdge(int from_x, int from_y, int from_z, int to_x, int to_y, int to_z, Edge* graphEdge, int firstOrLast = 0, int entryExitStatus = 0);
	void addStatusEdge(int from_x, int from_y, int from_z, int to_x, int to_y, int to_z, Edge* graphEdge, int firstOrLast = 0, int entryExitStatus = 0);
	void eraseStatusEdge(int from_x, int from_y, int from_z, int to_x, int to_y, int to_z, Edge* graphEdge);
	void eraseAllEdgeStatus();

	bool edgeExists(int from_x, int from_y, int from_z, int to_x, int to_y, int to_z, Edge* graphEdge);
	NWEdge* retrieveEdge(int from_x, int from_y, int from_z, int to_x, int to_y, int to_z, Edge* graphEdge);

	void addEdge(NWEdge* e);

	std::vector<Edge*> graphEdgesPath();

	void takeInExtensionPath(localExtension_pathDescription* pD, int entryExit, int completeToEnd);

	unsigned int getNumEdges()
	{
		return all_edges.size();
	}

	void recalculateFirstLast();

	void _printPath();

	void extendToZero(int whichEnd, int lastYCoordinate = 0);

};

class NWEdge {
public:
	int from_x;
	int from_y;
	int from_z;

	int to_x;
	int to_y;
	int to_z;

	Edge* usedGraphEdge;
	NWPath* path;

	bool scoreComputed;
	double scoreAfterEdge;
	NWEdge* scoreBacktrack;
	bool endsFree_previousEdgeAffineSequenceGap;

	void checkConsistency();
	NWEdge();

	void _printEdge();

	void makeExitEdge();

	double getScoreAfterEdge()
	{
		if(! scoreComputed)
		{
			_printEdge();
		}
		assert(scoreComputed);
		return scoreAfterEdge;
	}

	NWEdge* getScoreBacktrack()
	{
		assert(scoreComputed);
		return scoreBacktrack;
	}

	double calculateScore(bool comingFromAffineSequenceGap, bool comingFromAffineGraphGap, double S_match, double S_mismatch, double S_gapOpen, double S_gapExtend);
	double calculateScore_endsFree(bool comingFromAffineSequenceGap, bool comingFromAffineGraphGap, double S_match, double S_mismatch, double S_gapOpen, double S_gapExtend);

	void takeScore(double Score, NWEdge* backtrackEdge, double S_match, double S_mismatch, double S_gap);
	void takeScore_affine(double Score, NWEdge* backtrackEdge, double S_match, double S_mismatch, double S_gapOpen, double S_gapExtend);
	void takeScore_endsFree(double Score, bool previousEdgeAffineSequenceGap, NWEdge* backtrackEdge, double S_match, double S_mismatch, double S_gapOpen, double S_gapExtend);
	void takeScore_endsFree_nonCritical(double Score, bool previousEdgeAffineSequenceGap, NWEdge* backtrackEdge, double S_match, double S_mismatch, double S_gapOpen, double S_gapExtend, unsigned int* thisSeed);

	bool isSequenceGap();
	bool isSequenceGap_affine();
	bool isSequenceGap_endsFree();
	bool isGraphGap();
};

}
}

#endif /* VIRTUALNWUNIQUE_H_ */
