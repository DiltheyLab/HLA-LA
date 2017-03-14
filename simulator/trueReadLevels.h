/*
 * trueReadLevels.h
 *
 *  Created on: 01.10.2015
 *      Author: AlexanderDilthey
 */

#ifndef SIMULATOR_TRUEREADLEVELS_H_
#define SIMULATOR_TRUEREADLEVELS_H_

#include <map>
#include <vector>
#include <utility>
#include <string>

#include "../mapper/reads/verboseSeedChain.h"
#include "../mapper/reads/oneReadPair.h"
#include "../mapper/aligner/extensionAligner.h"

namespace simulator {

class trueReadLevels {
protected:
	std::map<std::string, std::pair<std::vector<int>, std::vector<int>>> trueLevels;
	std::map<std::string, std::pair<std::string, std::string>> true_underlyingEdgeLabels;

	std::map<std::string, std::pair<std::vector<int>, std::vector<int>>> true_fullAlignment_coordinates_edgePath;
	std::map<std::string, std::pair<std::string, std::string>> true_fullAlignment_underlyingEdgeLabels;
	std::map<std::string, std::pair<std::string, std::string>> true_fullAlignment_sequence;

	std::pair<size_t, size_t> total_and_correct;

public:
	trueReadLevels(std::string R1_levels, std::string R2_levels);
	virtual ~trueReadLevels();

	std::pair<size_t, size_t> evaluateAlignment(const mapper::reads::verboseSeedChainPair& alignment, const mapper::reads::oneReadPair& rP, mapper::aligner::extensionAligner* eA);

	std::pair<size_t, size_t> get_total_and_correct();

};


} /* namespace simulator */

#endif /* SIMULATOR_TRUEREADLEVELS_H_ */
