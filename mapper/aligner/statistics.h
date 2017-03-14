/*
 * statistics.h
 *
 *  Created on: Mar 12, 2017
 *      Author: diltheyat
 */

#ifndef MAPPER_ALIGNER_STATISTICS_H_
#define MAPPER_ALIGNER_STATISTICS_H_

#include <string>

namespace mapper {
namespace aligner {

class statistics {
public:
	statistics();
	virtual ~statistics();
	size_t complete_seeds;
	size_t incomplete_seeds;
	size_t complete_seeds_n_chains;
	size_t complete_seeds_chains_combined_length;
	size_t complete_seeds_n_primaryChains;
	size_t complete_seeds_primaryChain_combined_length;


	size_t n_calls_alignOneReadPair;
	size_t alignOneReadPair_considered_chains = 0;
	size_t alignOneReadPair_n_calledChainExtension = 0;
	size_t alignOneReadPair_n_calledChainExtension_primary = 0;

	size_t alignOneReadPair_return_selectedChains_totalColumns = 0;
	size_t alignOneReadPair_return_selectedChains_totalColumns_fromBWASeed = 0;

	size_t alignOneReadPair_selectedPair_removed_columns_noGap_restriction = 0;

	double alignOneReadPair_selectedPair_improvement_through_bt = 0;

	size_t readPairs_used_for_HLAtyping;

	size_t HLA_n_calls_alignOneReadPair;
	size_t HLA_alignOneReadPair_considered_chains;
	size_t HLA_alignOneReadPair_n_calledChainExtension;
	size_t HLA_alignOneReadPair_n_calledChainExtension_primary;

	size_t HLA_alignOneReadPair_return_selectedChains_totalColumns;
	size_t HLA_alignOneReadPair_return_selectedChains_totalColumns_fromBWASeed;

	size_t HLA_alignOneReadPair_selectedPair_removed_columns_noGap_restriction;

	double HLA_alignOneReadPair_selectedPair_improvement_through_bt;

	void printStatistics();

	void takeInHLARelatedDiff(statistics* Sbefore);

};


} /* namespace aligner */
} /* namespace mapper */

#endif /* MAPPER_ALIGNER_STATISTICS_H_ */
