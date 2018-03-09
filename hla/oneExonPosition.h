/*
 * oneExonPosition.h
 *
 *  Created on: 28.10.2015
 *      Author: AlexanderDilthey
 */

#ifndef HLA_ONEEXONPOSITION_H_
#define HLA_ONEEXONPOSITION_H_

#include <string>

namespace hla {

class oneExonPosition {
public:
	unsigned int positionInExon;
	int graphLevel;
	std::string genotype;
	std::string alignment_edgelabels;
	std::string qualities;

	std::string thisRead_ID;
	double thisRead_fractionOK;
	double thisRead_WeightedCharactersOK;

	std::string pairedRead_ID;
	double pairedRead_fractionOK;
	double pairedRead_WeightedCharactersOK;

	std::string read1_ID;

	bool pairs_strands_OK;
	double pairs_strands_distance;

	double mapQ;
	double mapQ_genomic;
	double mapQ_position;
	
	int alignmentColumnsWithAtLeastOneNonGap;
	
	int runningNovelGapEitherDirection;

	bool reverse;
};

} /* namespace hla */

#endif /* HLA_ONEEXONPOSITION_H_ */
