/*
 * extensionAligner.h
 *
 *  Created on: 26.09.2015
 *      Author: AlexanderDilthey
 */

#ifndef MAPPER_ALIGNER_EXTENSIONALIGNER_H_
#define MAPPER_ALIGNER_EXTENSIONALIGNER_H_

#include "alignerBase.h"

#include "../reads/oneRead.h"
#include "../reads/verboseSeedChain.h"
#include "../../Graph/Graph.h"
#include "VirtualNWUnique.h"
#include "../reads/protoSeeds.h"

#include <string>

namespace mapper {
namespace aligner {

class extensionAligner : public alignerBase {
protected:
public:
	extensionAligner(Graph* g_);
	virtual ~extensionAligner();

	mapper::reads::verboseSeedChain extendSeedChain(const std::string& sequence, const mapper::reads::verboseSeedChain& seedChain) const;
	std::vector<mapper::reads::verboseSeedChain> fullNeedleman_diagonal_extension(const std::string& sequence, int start_sequence, int startLevel_graph, int startZ_graph, int maxLevel_graph, int maxPosition_sequence, int diagonal_stop_threshold, VirtualNWTable_Unique* blockedPathsTable, bool directionPositive, bool returnGlobalScore, bool preferSequenceCompleAlignments) const;
	std::vector<mapper::reads::verboseSeedChain> fullNeedleman_diagonal_extension_gapJumper(const std::string& sequence, int start_sequence, int startLevel_graph, int startZ_graph, int maxLevel_graph, int maxPosition_sequence, int diagonal_stop_threshold, VirtualNWTable_Unique* blockedPathsTable, bool directionPositive, bool returnGlobalScore, bool preferSequenceCompleAlignments) const;

	double scoreOneAlignment(const mapper::reads::verboseSeedChain& alignment, const mapper::reads::oneRead& underlyingRead, std::string longReadMode = "") const;

	bool paranoid;

	mutable std::vector<unsigned int> rng_seeds;

	unsigned int threads_initialized;

	std::vector<std::pair<int, int>> _graph_get_jumpPrevious_x_and_z_values(int x, int z) const;
	std::vector<std::pair<int, std::pair<int, Edge*>>> _graph_get_jumpNext_x_and_z_values_and_edges(int x, int z) const;
	std::vector<std::pair<int, std::pair<int, Edge*>>> _graph_get_jumpPrevious_x_and_z_values_and_edges(int x, int z) const;

	void init_for_threads(unsigned int threads);
};

} /* namespace aligner */
} /* namespace mapper */

#endif /* MAPPER_ALIGNER_EXTENSIONALIGNER_H_ */
