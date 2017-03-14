/*
 * alignmentContext.h
 *
 *  Created on: 28.09.2015
 *      Author: AlexanderDilthey
 */

#ifndef MAPPER_ALIGNER_ALIGNMENTCONTEXT_H_
#define MAPPER_ALIGNER_ALIGNMENTCONTEXT_H_

#include "alignerBase.h"
#include <string>

namespace mapper {
namespace aligner {

class alignmentContext {
public:
	alignmentContext();
	virtual ~alignmentContext();

	const alignerBase* aligner;
	std::string const * sequence;
};

} /* namespace aligner */
} /* namespace mapper */

#endif /* MAPPER_ALIGNER_ALIGNMENTCONTEXT_H_ */
