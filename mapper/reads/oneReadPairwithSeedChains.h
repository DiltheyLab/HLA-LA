/*
 * oneReadPairwithSeedChains.h
 *
 *  Created on: 22.09.2015
 *      Author: AlexanderDilthey
 */

#ifndef MAPPER_READS_ONEREADPAIRWITHSEEDCHAINS_H_
#define MAPPER_READS_ONEREADPAIRWITHSEEDCHAINS_H_

#include <vector>
#include <utility>

#include "oneReadPair.h"
#include "seedChain.h"

namespace mapper {
namespace reads {

class oneReadPair_withSeedChains {
public:
	oneReadPair rP;
	std::pair<std::vector<seedChain>, std::vector<seedChain>> seeds;

	oneReadPair_withSeedChains();
	virtual ~oneReadPair_withSeedChains();
};

} /* namespace reads */
} /* namespace mapper */

#endif /* MAPPER_READS_ONEREADPAIRWITHSEEDCHAINS_H_ */
