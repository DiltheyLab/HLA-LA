/*
 * alignerBase.h
 *
 *  Created on: 26.09.2015
 *      Author: AlexanderDilthey
 */

#ifndef MAPPER_ALIGNER_ALIGNERBASE_H_
#define MAPPER_ALIGNER_ALIGNERBASE_H_

#include <vector>
#include <map>
#include <utility>
#include <string>

#include "../../Graph/Graph.h"
#include "../reads/verboseSeedChain.h"
#include "api/BamAlignment.h"

namespace mapper {
namespace aligner {

class backtraceStep {
public:
	int x;
	int y;
	int z;
	Edge* usedEdge;

	backtraceStep() : x(-1), y(-1), z(-1), usedEdge(0)
	{

	}
};

class backtraceStep_affine {
public:
	int x;
	int y;
	int z;
	int sourceMatrix;
	Edge* usedEdge;
	bool _border_lastStep_affine;

	backtraceStep_affine() : x(-1), y(-1), z(-1), sourceMatrix(-1), usedEdge(0), _border_lastStep_affine(false)
	{

	}
};

class alignerBase {
public:
	alignerBase(Graph* g_);
	virtual ~alignerBase();

	Graph* g;
	double S_match;
	double S_mismatch;
	double S_gap;
	double S_openGap;
	double S_extendGap;
	double S_graphGap;

	std::vector<std::vector<Node*> > nodesPerLevel_ordered;
	std::vector<std::map<Node*, unsigned int> > nodesPerLevel_ordered_rev;

	std::vector<std::pair<int, Edge*> > _graph_get_previous_z_values_and_edges(int x, int z) const;
	std::vector<std::pair<int, Edge*> > _graph_get_next_z_values_and_edges(int x, int z) const;
	std::vector<int> _graph_get_previous_z_values(int x, int z) const;

	static bool alignedReadPair_strandsValid(const mapper::reads::verboseSeedChain& alignment_read1, const mapper::reads::verboseSeedChain& alignment_read2);
	static int alignedReadPair_pairsDistanceInGraphLevels(const mapper::reads::verboseSeedChain& alignment_read1, const mapper::reads::verboseSeedChain& alignment_read2);

	static bool alignedReadPair_strandsValid(const mapper::reads::verboseSeedChainPair& alignment)
	{
		return alignedReadPair_strandsValid(alignment.chains.first, alignment.chains.second);
	};

	static int alignedReadPair_pairsDistanceInGraphLevels(const mapper::reads::verboseSeedChainPair& alignment)
	{
		return alignedReadPair_pairsDistanceInGraphLevels(alignment.chains.first, alignment.chains.second);
	}


	std::set<int> alignedReadPair_pairsDistancesUnderlyingSequences(const mapper::reads::verboseSeedChain& alignment_read1, const mapper::reads::verboseSeedChain& alignment_read2, const std::vector<std::map<int, int>>& graphLevel_2_underlyingSequencePositions) const;

	std::set<int> alignedReadPair_properGraphDistance(const mapper::reads::verboseSeedChain& alignment_read1, const mapper::reads::verboseSeedChain& alignment_read2) const;
	std::set<int> properGraphDistance(Node* n1, Node* n2) const;
};

} /* namespace aligner */
} /* namespace mapper */

#endif /* MAPPER_ALIGNER_ALIGNERBASE_H_ */
