/*
 * BAMalignment.h
 *
 *  Created on: 22.09.2015
 *      Author: AlexanderDilthey
 */

#ifndef MAPPER_READS_PRGCONTIGBAMALIGNMENT_H_
#define MAPPER_READS_PRGCONTIGBAMALIGNMENT_H_

#include <string>
#include <vector>

namespace mapper {
namespace reads {

class PRGContigBAMAlignment {
public:
	std::vector<int> graph_aligned_levels;
	std::string graph_aligned;
	std::string sequence_aligned;
	int sequence_aligned_startInRaw;
	int sequence_aligned_stopInRaw;
	bool reverse;

	PRGContigBAMAlignment();
	virtual ~PRGContigBAMAlignment();

	void removeSequenceCharacters(bool left, int n);

	void print();

	void checkAlignmentConcordanceWithSequence(const std::string& sequence) const;
};

} /* namespace reads */
} /* namespace mapper */

#endif /* MAPPER_READS_PRGCONTIGBAMALIGNMENT_H_ */
