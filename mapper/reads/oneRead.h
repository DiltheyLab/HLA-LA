/*
 * oneRead.h
 *
 *  Created on: 21.09.2015
 *      Author: AlexanderDilthey
 */

#ifndef MAPPER_READS_ONEREAD_H_
#define MAPPER_READS_ONEREAD_H_

#include <assert.h>
#include <string>
#include <algorithm>
#include "../../Utilities.h"

namespace mapper {
namespace reads {

class oneRead {
public:
	std::string name;
	std::string sequence;
	std::string quality;

	oneRead()
	{

	}

	oneRead(std::string read_name, std::string read_sequence, std::string read_qualities) : name(read_name), sequence(read_sequence), quality(read_qualities)
	{
		assert(read_sequence.length() == read_qualities.length());
	}

	void invert()
	{
		sequence = Utilities::seq_reverse_complement(sequence);
		std::reverse(quality.begin(), quality.end());
	}

	virtual ~oneRead();
};

} /* namespace reads */
} /* namespace mapper */

#endif /* MAPPER_READS_ONEREAD_H_ */
