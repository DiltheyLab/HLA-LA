/*
 * protoSeeds.h
 *
 *  Created on: 21.09.2015
 *      Author: AlexanderDilthey
 */

#ifndef MAPPER_PROTOSEEDS_H_
#define MAPPER_PROTOSEEDS_H_

#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include <utility>
#include <tuple>

namespace mapper {

namespace reads {
class protoSeeds {
friend class seed;

public:
	protoSeeds();
	virtual ~protoSeeds();

	void takeAlignment(const BamTools::BamAlignment& a, int which, std::string referenceID, int reference2level_offset_0based, int whichReader);
	bool isComplete() const;

	void printDebug(const BamTools::BamReader& R) const;

	size_t read1_getPrimaryAlignmentI() const;
	size_t read2_getPrimaryAlignmentI() const;

	std::vector<std::tuple<std::string, int, BamTools::BamAlignment, int>> read1_alignments;
	std::vector<std::tuple<std::string, int, BamTools::BamAlignment, int>> read2_alignments;

	static void printAlignmentInfo(size_t i, const BamTools::BamAlignment& al);

	void refreshPrimaryStatus();

protected:


	bool read1_havePrimary;
	bool read2_havePrimary;



};

}
} /* namespace mapper */

#endif /* MAPPER_PROTOSEEDS_H_ */
