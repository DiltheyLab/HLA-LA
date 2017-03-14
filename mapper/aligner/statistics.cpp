/*
 * statistics.cpp
 *
 *  Created on: Mar 12, 2017
 *      Author: diltheyat
 */

#include "statistics.h"
#include <iostream>

namespace mapper {
namespace aligner {

statistics::statistics() {
	// TODO Auto-generated constructor stub

	complete_seeds = 0;
	incomplete_seeds = 0;
	complete_seeds_n_chains = 0;
	complete_seeds_chains_combined_length = 0;
	complete_seeds_n_primaryChains = 0;
	complete_seeds_primaryChain_combined_length = 0;

	n_calls_alignOneReadPair = 0;
	alignOneReadPair_considered_chains = 0;
	alignOneReadPair_n_calledChainExtension = 0;
	alignOneReadPair_n_calledChainExtension_primary = 0;

	alignOneReadPair_return_selectedChains_totalColumns = 0;
	alignOneReadPair_return_selectedChains_totalColumns_fromBWASeed = 0;

	alignOneReadPair_selectedPair_removed_columns_noGap_restriction = 0;

	alignOneReadPair_selectedPair_improvement_through_bt = 0;

	readPairs_used_for_HLAtyping = 0;

	HLA_n_calls_alignOneReadPair = 0;
	HLA_alignOneReadPair_considered_chains = 0;
	HLA_alignOneReadPair_n_calledChainExtension = 0;
	HLA_alignOneReadPair_n_calledChainExtension_primary = 0;

	HLA_alignOneReadPair_return_selectedChains_totalColumns = 0;
	HLA_alignOneReadPair_return_selectedChains_totalColumns_fromBWASeed = 0;
	HLA_alignOneReadPair_selectedPair_removed_columns_noGap_restriction = 0;
	HLA_alignOneReadPair_selectedPair_improvement_through_bt = 0;
}

statistics::~statistics() {
	// TODO Auto-generated destructor stub
}

void statistics::printStatistics()
{
	size_t n_reads = complete_seeds + incomplete_seeds;

	std::cout << "Output from mapper::aligner::statistics" << "\n";
	std::cout << "=======================================" << "\n\n";

	std::cout << "Read / chain statistics:\n";
	std::cout << "\t" << "Seen reads" << ": " << n_reads << "\n";
	std::cout << "\t" << "Complete chains" << ": " << complete_seeds << " (" << (float)complete_seeds/n_reads << ")\n";
	std::cout << "\t" << "Incomplete chains" << ": " << incomplete_seeds << " (" << (float)complete_seeds/n_reads << ")\n";

	std::cout << "For complete seeds:\n";
	std::cout << "\t" << "Total chains" << ": " << complete_seeds_n_chains << ", i.e. " << (float)complete_seeds_n_chains/complete_seeds << " per complete seed." << "\n";
	std::cout << "\t" << "Average chain length: " << (float)complete_seeds_chains_combined_length/complete_seeds_n_chains << "\n";
	std::cout << "\t" << "Total primary-derived chains" << ": " << complete_seeds_n_primaryChains << ", i.e. " << (float)complete_seeds_n_primaryChains/complete_seeds << " per complete seed." << "\n";
	std::cout << "\t" << "Average primary-derived chain length: " << (float)complete_seeds_primaryChain_combined_length/complete_seeds_n_primaryChains << "\n";

	std::cout << "(All reads) alignOneReadPair stats:\n";
	std::cout << "\t" << "Calls (i.e. read pairs)" << ": " << n_calls_alignOneReadPair << "\n";
	std::cout << "\t" << "Total considered chains" << ": " << alignOneReadPair_considered_chains << ", i.e. " << (float)alignOneReadPair_considered_chains/n_calls_alignOneReadPair << " per call." << "\n";
	std::cout << "\t" << "Called chain extension" << ": " << alignOneReadPair_n_calledChainExtension << ", i.e. " << (float)alignOneReadPair_n_calledChainExtension/n_calls_alignOneReadPair << " per call." << "\n";
	std::cout << "\t" << "Called chain extension on primary" << ": " << alignOneReadPair_n_calledChainExtension_primary << ", i.e. " << (float)alignOneReadPair_n_calledChainExtension_primary/alignOneReadPair_n_calledChainExtension << " of calls are on primary." << "\n";

	std::cout << "(All reads) alignOneReadPair output stats:\n";
	std::cout << "\t" << "Total removed columns due to gap-related checks" << ": " << alignOneReadPair_selectedPair_removed_columns_noGap_restriction << ", i.e. " << (float)alignOneReadPair_selectedPair_removed_columns_noGap_restriction/n_calls_alignOneReadPair << " per read pair/call." << "\n";
	std::cout << "\t" << "Total improvement due to backtrace" << ": " << alignOneReadPair_selectedPair_improvement_through_bt << ", i.e. " << (float)alignOneReadPair_selectedPair_improvement_through_bt/n_calls_alignOneReadPair << " per read pair/call." << "\n";
	std::cout << "\t" << "Total returned columns pre-extension" << ": " << alignOneReadPair_return_selectedChains_totalColumns_fromBWASeed << ", i.e. " << (float)alignOneReadPair_return_selectedChains_totalColumns_fromBWASeed/n_calls_alignOneReadPair << " per read pair/call." << "\n";
	std::cout << "\t" << "Total returned columns" << ": " << alignOneReadPair_return_selectedChains_totalColumns << ", i.e. " << (float)alignOneReadPair_return_selectedChains_totalColumns/n_calls_alignOneReadPair << " per read pair/call." << "\n";

	std::cout << "For HLA typing:\n";
	std::cout << "\t" << "Reads used" << ": " << readPairs_used_for_HLAtyping << "\n";

	std::cout << "(HLA) alignOneReadPair stats:\n";
	std::cout << "\t" << "Calls (i.e. read pairs)" << ": " << HLA_n_calls_alignOneReadPair << "\n";
	std::cout << "\t" << "Total considered chains" << ": " << HLA_alignOneReadPair_considered_chains << ", i.e. " << (float)HLA_alignOneReadPair_considered_chains/HLA_n_calls_alignOneReadPair << " per call." << "\n";
	std::cout << "\t" << "Called chain extension" << ": " << HLA_alignOneReadPair_n_calledChainExtension << ", i.e. " << (float)HLA_alignOneReadPair_n_calledChainExtension/HLA_n_calls_alignOneReadPair << " per call." << "\n";
	std::cout << "\t" << "Called chain extension on primary" << ": " << HLA_alignOneReadPair_n_calledChainExtension_primary << ", i.e. " << (float)HLA_alignOneReadPair_n_calledChainExtension_primary/HLA_alignOneReadPair_n_calledChainExtension << " of calls are on primary." << "\n";

	std::cout << "(HLA) alignOneReadPair output stats:\n";
	std::cout << "\t" << "Total removed columns due to gap-related checks" << ": " << HLA_alignOneReadPair_selectedPair_removed_columns_noGap_restriction << ", i.e. " << (float)HLA_alignOneReadPair_selectedPair_removed_columns_noGap_restriction/HLA_n_calls_alignOneReadPair << " per read pair/call." << "\n";
	std::cout << "\t" << "Total improvement due to backtrace" << ": " << HLA_alignOneReadPair_selectedPair_improvement_through_bt << ", i.e. " << (float)HLA_alignOneReadPair_selectedPair_improvement_through_bt/HLA_n_calls_alignOneReadPair << " per read pair/call." << "\n";
	std::cout << "\t" << "Total returned columns pre-extension" << ": " << HLA_alignOneReadPair_return_selectedChains_totalColumns_fromBWASeed << ", i.e. " << (float)HLA_alignOneReadPair_return_selectedChains_totalColumns_fromBWASeed/HLA_n_calls_alignOneReadPair << " per read pair/call." << "\n";
	std::cout << "\t" << "Total returned columns" << ": " << HLA_alignOneReadPair_return_selectedChains_totalColumns << ", i.e. " << (float)HLA_alignOneReadPair_return_selectedChains_totalColumns/HLA_n_calls_alignOneReadPair << " per read pair/call." << "\n";


	std::cout << "\n==\n\n" << std::flush;

}

void statistics::takeInHLARelatedDiff(statistics* Sbefore)
{

	size_t d_n_calls_alignOneReadPair = n_calls_alignOneReadPair - Sbefore->n_calls_alignOneReadPair;
	size_t d_alignOneReadPair_considered_chains = alignOneReadPair_considered_chains - Sbefore->alignOneReadPair_considered_chains;
	size_t d_alignOneReadPair_n_calledChainExtension = alignOneReadPair_n_calledChainExtension - Sbefore->alignOneReadPair_n_calledChainExtension;
	size_t d_alignOneReadPair_n_calledChainExtension_primary = alignOneReadPair_n_calledChainExtension_primary - Sbefore->alignOneReadPair_n_calledChainExtension_primary;
	size_t d_alignOneReadPair_return_selectedChains_totalColumns = alignOneReadPair_return_selectedChains_totalColumns - Sbefore->alignOneReadPair_return_selectedChains_totalColumns;
	size_t d_alignOneReadPair_return_selectedChains_totalColumns_fromBWASeed = alignOneReadPair_return_selectedChains_totalColumns_fromBWASeed - Sbefore->alignOneReadPair_return_selectedChains_totalColumns_fromBWASeed;
	size_t d_alignOneReadPair_selectedPair_removed_columns_noGap_restriction = alignOneReadPair_selectedPair_removed_columns_noGap_restriction - Sbefore->alignOneReadPair_selectedPair_removed_columns_noGap_restriction;
	double d_alignOneReadPair_selectedPair_improvement_through_bt = alignOneReadPair_selectedPair_improvement_through_bt - Sbefore->alignOneReadPair_selectedPair_improvement_through_bt;

	HLA_n_calls_alignOneReadPair += d_n_calls_alignOneReadPair;
	HLA_alignOneReadPair_considered_chains += d_alignOneReadPair_considered_chains;
	HLA_alignOneReadPair_n_calledChainExtension += d_alignOneReadPair_n_calledChainExtension;
	HLA_alignOneReadPair_n_calledChainExtension_primary += d_alignOneReadPair_n_calledChainExtension_primary;
	HLA_alignOneReadPair_return_selectedChains_totalColumns += d_alignOneReadPair_return_selectedChains_totalColumns;
	HLA_alignOneReadPair_return_selectedChains_totalColumns_fromBWASeed += d_alignOneReadPair_return_selectedChains_totalColumns_fromBWASeed;
	HLA_alignOneReadPair_selectedPair_removed_columns_noGap_restriction += d_alignOneReadPair_selectedPair_removed_columns_noGap_restriction;
	HLA_alignOneReadPair_selectedPair_improvement_through_bt += d_alignOneReadPair_selectedPair_improvement_through_bt;
}

} /* namespace aligner */
} /* namespace mapper */
