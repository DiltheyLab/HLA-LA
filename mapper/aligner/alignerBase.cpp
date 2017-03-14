/*
 * alignerBase.cpp
 *
 *  Created on: 26.09.2015
 *      Author: AlexanderDilthey
 */

#include "alignerBase.h"

#include <assert.h>
#include <iostream>

namespace mapper {
namespace aligner {

alignerBase::alignerBase(Graph* g_) {
	g = g_;

	S_match = 2;
	S_mismatch = -5;
	S_gap = -2;
	S_graphGap = 0;

	S_openGap = -4;
	S_extendGap = -2;

	unsigned int levels = g->NodesPerLevel.size();
	nodesPerLevel_ordered.resize(levels);
	nodesPerLevel_ordered_rev.resize(levels);
	for(unsigned int levelI = 0; levelI < levels; levelI++)
	{
		nodesPerLevel_ordered.at(levelI) = std::vector<Node*>(g->NodesPerLevel.at(levelI).begin(), g->NodesPerLevel.at(levelI).end());
		for(unsigned int nodeI = 0; nodeI < nodesPerLevel_ordered.at(levelI).size(); nodeI++)
		{
			nodesPerLevel_ordered_rev.at(levelI)[nodesPerLevel_ordered.at(levelI).at(nodeI)] = nodeI;
		}
	}
}

alignerBase::~alignerBase() {
	// TODO Auto-generated destructor stub
}


std::set<int> alignerBase::alignedReadPair_properGraphDistance(const mapper::reads::verboseSeedChain& alignment_read1, const mapper::reads::verboseSeedChain& alignment_read2) const
{

	std::set<int> forReturn;
	int scanPositions = 2;
	if(alignment_read1.alignment_firstLevel() < alignment_read2.alignment_firstLevel())
	{
		// read1 << read2
		Node* read1_endNode = alignment_read1.alignment_lastNode();
		Node* read2_beginNode = alignment_read2.alignment_firstNode();

		return properGraphDistance(read1_endNode, read2_beginNode);
	}
	else
	{
		assert(alignment_read1.alignment_firstLevel() >= alignment_read2.alignment_firstLevel());
		// read1 << read2
		Node* read2_endNode = alignment_read2.alignment_lastNode();
		Node* read1_beginNode = alignment_read1.alignment_firstNode();

		return properGraphDistance(read2_endNode, read1_beginNode);
	}
}

std::set<int> alignerBase::properGraphDistance(Node* n1, Node* n2) const
{
	if(n2->level < n1->level)
	{
		std::set<int> distance_rightOrientiation = properGraphDistance(n2, n1);
		std::set<int> forReturn;
		for(std::set<int>::iterator distanceIt = distance_rightOrientiation.begin(); distanceIt != distance_rightOrientiation.end(); distanceIt++)
		{
			int distance = *distanceIt;
			forReturn.insert(-1 * distance);
		}
		return forReturn;
	}

	if(n1->level == n2->level)
	{
		std::set<int> forReturn;
		if(n1 == n2)
		{
			forReturn.insert(0);
		}
		return forReturn;
	}

	assert(n1->level < n2->level);

	std::map<Node*, std::pair<int, int>> runningNodes;
	runningNodes[n1] = make_pair(0, 0);
	for(int levelI = (int)n1->level; levelI < (int)n2->level; levelI++)
	{
		std::map<Node*, std::pair<int, int>> runningNodes_nextLevel;
		for(std::map<Node*, std::pair<int, int>>::iterator nodeIt = runningNodes.begin(); nodeIt != runningNodes.end(); nodeIt++)
		{
			Node* originNode = nodeIt->first;
			const std::set<Edge*>& edges = nodeIt->first->Outgoing_Edges;
			for(std::set<Edge*>::const_iterator edgeIt = edges.begin(); edgeIt != edges.end(); edgeIt++)
			{
				Node* targetNode = (*edgeIt)->To;
				bool isGapEdge = ((*edgeIt)->getEmission() == "_");

				if(runningNodes_nextLevel.count(targetNode) == 0)
				{
					if(isGapEdge)
					{
						runningNodes_nextLevel[targetNode] = runningNodes.at(originNode);
					}
					else
					{
						runningNodes_nextLevel[targetNode] = make_pair(runningNodes.at(originNode).first + 1, runningNodes.at(originNode).second + 1);
					}
				}
				else
				{
					int nextValue_min = runningNodes.at(originNode).first + (isGapEdge ? 0 : 1);
					int nextValue_max = runningNodes.at(originNode).second + (isGapEdge ? 0 : 1);

					if(nextValue_min < runningNodes_nextLevel[targetNode].first)
					{
						runningNodes_nextLevel[targetNode].first = nextValue_min;
					}

					if(nextValue_max > runningNodes_nextLevel[targetNode].second)
					{
						runningNodes_nextLevel[targetNode].second = nextValue_max;
					}
				}
			}
		}
		runningNodes = runningNodes_nextLevel;
	}

	std::set<int> forReturn;
	if(runningNodes.count(n2))
	{
		forReturn.insert(runningNodes.at(n2).first);
		forReturn.insert(runningNodes.at(n2).second);
	}
	return forReturn;
}

std::vector<std::pair<int, Edge*> > alignerBase::_graph_get_previous_z_values_and_edges(int x, int z) const
{
	assert(g != 0);
	std::vector<std::pair<int, Edge*> > forReturn;
	unsigned int uZ = (unsigned int)z;
	assert(x > 0);
	assert(x < (int)nodesPerLevel_ordered.size());
	assert(z >= 0);
	assert(z < (int)nodesPerLevel_ordered.at(x).size());
	Node* thisZ = nodesPerLevel_ordered.at(x).at(uZ);
	set<Edge*> incomingEdges = thisZ->Incoming_Edges;
	set<Node*> nodesPreviousLevel;
	assert(incomingEdges.size() > 0);
	for(std::set<Edge*>::iterator eIt = incomingEdges.begin(); eIt != incomingEdges.end(); eIt++)
	{
		Edge* e = *eIt;
		Node* fromNode = e->From;
		int z_for_fromNode = nodesPerLevel_ordered_rev.at(x-1).at(fromNode);
		forReturn.push_back(std::pair<int, Edge*>(z_for_fromNode, e));
	}

	return forReturn;
}

std::vector<std::pair<int, Edge*> > alignerBase::_graph_get_next_z_values_and_edges(int x, int z) const
{
	assert(g != 0);
	std::vector<std::pair<int, Edge*> > forReturn;
	unsigned int uZ = (unsigned int)z;
	assert(x >= 0);
	assert(x < (int)nodesPerLevel_ordered.size());
	assert(z >= 0);
	assert(z < (int)nodesPerLevel_ordered.at(x).size());
	Node* thisZ = nodesPerLevel_ordered.at(x).at(uZ);
	set<Edge*> outgoingEdges = thisZ->Outgoing_Edges;
	set<Node*> nodesNextLevel;
	assert(outgoingEdges.size() > 0);
	for(std::set<Edge*>::iterator eIt = outgoingEdges.begin(); eIt != outgoingEdges.end(); eIt++)
	{
		Edge* e = *eIt;
		Node* toNode = e->To;
		int z_for_toNode = nodesPerLevel_ordered_rev.at(x+1).at(toNode);
		forReturn.push_back(std::pair<int, Edge*>(z_for_toNode, e));
	}

	return forReturn;
}

std::vector<int> alignerBase::_graph_get_previous_z_values(int x, int z) const
{
	assert(g != 0);
	std::vector<std::pair<int, Edge*> > previousZs = _graph_get_previous_z_values_and_edges(x, z);

	std::set<int> forReturn;
	for(unsigned int i = 0; i < previousZs.size(); i++)
	{
		int previousZ = previousZs.at(i).first;
		forReturn.insert(previousZ);
	}

	return std::vector<int>(forReturn.begin(), forReturn.end());
}


bool alignerBase::alignedReadPair_strandsValid(const mapper::reads::verboseSeedChain& alignment_read1, const mapper::reads::verboseSeedChain& alignment_read2)
{
	if ((alignment_read1.alignment_firstLevel() != -1) && (alignment_read2.alignment_firstLevel() != -1) && (alignment_read1.reverse != alignment_read2.reverse))
	{
		if(! alignment_read1.reverse)
		{
			if(alignment_read1.alignment_firstLevel() < alignment_read2.alignment_firstLevel())
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			if(alignment_read1.alignment_lastLevel() > alignment_read2.alignment_lastLevel())
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}
	else
	{
		return false;
	}
}

int alignerBase::alignedReadPair_pairsDistanceInGraphLevels(const mapper::reads::verboseSeedChain& alignment_read1, const mapper::reads::verboseSeedChain& alignment_read2)
{
	// std::cerr << "alignedReadPair_pairsDistanceInGraphLevels(..)\n";
	// std::cerr << "\t" << "first reverse: " << p.first.reverse << "\n";
	// std::cerr << "\t" << "first coordinates: " << p.first.alignment_firstLevel() << " - " << p.first.alignment_lastLevel() << "\n";
	// std::cerr << "\t" << "second reverse: " << p.second.reverse << "\n";
	// std::cerr << "\t" << "second coordinates: " << p.second.alignment_firstLevel() << " - " << p.second.alignment_lastLevel() << "\n";

	if(alignment_read1.alignment_firstLevel() < alignment_read2.alignment_firstLevel())
	{
		int D = (alignment_read2.alignment_firstLevel() - alignment_read1.alignment_lastLevel() - 1);
		// std::cerr << "\t" << "first in front -- distance " << D << "\n";

		// if(alignedReadPair_strandsValid(p) && (! p.first.reverse))
		// {
			// std::cerr << "\t" << "OK" << "\n";
		// }
		// else
		// {
			// std::cerr << "\t" << "WARNING!" << "\n";
		// }
		// std::cerr << std::flush;
		return D;

	}
	else
	{
		assert(alignment_read1.alignment_firstLevel() >= alignment_read2.alignment_firstLevel());
		int D = (alignment_read1.alignment_firstLevel() - alignment_read2.alignment_lastLevel() - 1);
		// std::cerr << "\t" << "second in front -- distance " << D << "\n";

		// if(alignedRead Pair_strandsValid(p) && (p.first.reverse))
		// {
			// std::cerr << "\t" << "OK" << "\n";
		// }
		// else
		// {
			// std::cerr << "\t" << "WARNING!" << "\n";
		// }
		// std::cerr << std::flush;
		return D;
	}
}

std::set<int> alignerBase::alignedReadPair_pairsDistancesUnderlyingSequences(const mapper::reads::verboseSeedChain& alignment_read1, const mapper::reads::verboseSeedChain& alignment_read2, const std::vector<std::map<int, int>>& graphLevel_2_underlyingSequencePositions) const
{
	std::set<int> forReturn;
	int scanPositions = 2;
	if(alignment_read1.alignment_firstLevel() < alignment_read2.alignment_firstLevel())
	{
		// read1 << read2
		std::map<int, int> read1_endAnchors = alignment_read1.alignment_end_originalSequenceAnchors(scanPositions, graphLevel_2_underlyingSequencePositions);
		std::map<int, int> read2_beginAnchors = alignment_read2.alignment_begin_originalSequenceAnchors(scanPositions, graphLevel_2_underlyingSequencePositions);

		for(std::map<int, int>::iterator sequenceIDit = read1_endAnchors.begin(); sequenceIDit != read1_endAnchors.end(); sequenceIDit++)
		{
			int sequenceID = sequenceIDit->first;
			if(read2_beginAnchors.count(sequenceIDit->first))
			{
				int D = (read2_beginAnchors.at(sequenceID) - read1_endAnchors.at(sequenceID) - 1);
				forReturn.insert(D);
			}
		}
	}
	else
	{
		assert(alignment_read1.alignment_firstLevel() >= alignment_read2.alignment_firstLevel());
		// read2 << read1
		std::map<int, int> read2_endAnchors = alignment_read2.alignment_end_originalSequenceAnchors(scanPositions, graphLevel_2_underlyingSequencePositions);
		std::map<int, int> read1_beginAnchors = alignment_read1.alignment_begin_originalSequenceAnchors(scanPositions, graphLevel_2_underlyingSequencePositions);

		for(std::map<int, int>::iterator sequenceIDit = read2_endAnchors.begin(); sequenceIDit != read2_endAnchors.end(); sequenceIDit++)
		{
			int sequenceID = sequenceIDit->first;
			if(read1_beginAnchors.count(sequenceIDit->first))
			{
				int D = (read1_beginAnchors.at(sequenceID) - read2_endAnchors.at(sequenceID) - 1);
				forReturn.insert(D);
			}
		}
	}

	return forReturn;
}
} /* namespace aligner */
} /* namespace mapper */
