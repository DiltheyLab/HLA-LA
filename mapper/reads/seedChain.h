/*
 * seedChain.h
 *
 *  Created on: 21.09.2015
 *      Author: AlexanderDilthey
 */

#ifndef MAPPER_SEEDCHAIN_H_
#define MAPPER_SEEDCHAIN_H_

#include <vector>
#include "../../Graph/Edge.h"

namespace mapper {
namespace reads {
class seedChain {
public:
	int sequence_begin;
	int sequence_end;

	bool reverse;

	std::vector<Edge*> traversedEdges;
	seedChain() : sequence_begin(-1), sequence_end(-1), reverse(false)
	{

	}

	virtual ~seedChain();
};
}
} /* namespace mapper */

#endif /* MAPPER_SEEDCHAIN_H_ */

