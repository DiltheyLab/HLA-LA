/*
 * oneReadPair.h
 *
 *  Created on: 21.09.2015
 *      Author: AlexanderDilthey
 */

#ifndef MAPPER_READS_ONEREADPAIR_H_
#define MAPPER_READS_ONEREADPAIR_H_

#include "oneRead.h"
#include <utility>

namespace mapper {
namespace reads {

class oneReadPair {
public:
	oneReadPair();
	virtual ~oneReadPair();

	std::pair<oneRead, oneRead> reads;

	oneReadPair(oneRead r1, oneRead r2, unsigned int difference_starting_coordinates) : reads(std::pair<oneRead, oneRead>(r1, r2))
	{

	}

	void invert()
	{
		oneRead t("", "", "");
		t = reads.first;
		reads.first = reads.second;
		reads.second = t;
		reads.first.invert();
		reads.second.invert();
	}
};

} /* namespace reads */
} /* namespace mapper */

#endif /* MAPPER_READS_ONEREADPAIR_H_ */
