/*
 * fullLengthHMM.h
 *
 *  Created on: 24.01.2022
 *      Author: Alexa
 */

#ifndef FULLLENGTHHMM_FULLLENGTHHMM_H_
#define FULLLENGTHHMM_FULLLENGTHHMM_H_

#include <vector>
#include <string>
#include <map>
#include <set>
#include <utility>


struct HMMstate
{
	unsigned int level;
	std::pair<std::string, std::string> haplotypes_alleles;
	std::pair<unsigned int, unsigned int> copyingFrom;
	size_t readAssignmentState;
};

class fullLengthHMM {
protected:
	std::map<std::string, unsigned int> gene_length;
	std::map<std::string, std::string> reads_2_genes;
	std::map<std::string, std::map<std::string, std::pair<unsigned int, unsigned int>>> read_start_stop_positions;
	std::map<std::string, std::map<unsigned int, std::set<std::string>>> read_start_per_position;
	std::map<std::string, std::map<unsigned int, std::set<std::string>>> read_stop_per_position;
	std::map<std::string, std::map<std::string, std::map<unsigned int, std::string>>> read_genotypes_per_position;
	std::map<std::string, std::map<unsigned int, std::set<std::string>>> activeAlleles_per_position;
	std::map<std::string, std::map<std::string, std::string>> MSA_reference_sequences;
	std::map<std::string, std::map<std::string, std::string>> MSA_reference_sequences_whichHap;


	std::vector<std::string> readAssignmentStates;
	std::map<std::string, size_t> readAssignmentState_2_index;

	std::vector<std::map<size_t, std::set<size_t>>> level_readAssignmentState_2_states;

	std::map<std::string, std::size_t> readID_2_index;
	std::vector<std::string> computeReadAssignmentSets(const std::set<std::string>& runningReadIDs);
	std::set<std::string> nextLevel_compatibleReadAssignments(const std::string& thisReadAssignment, const std::set<size_t>& disappearingReadIDs_nextLevel, const std::set<size_t>& newReadIDs_nextLevel);
	std::string readAssignmentTemplate;

public:
	fullLengthHMM(
		std::map<std::string, unsigned int> _gene_length,
		std::map<std::string, std::string> _reads_2_genes,
		std::map<std::string, std::map<std::string, std::pair<unsigned int, unsigned int>>> _read_start_stop_positions,
		std::map<std::string, std::map<unsigned int, std::set<std::string>>> _read_start_per_position,
		std::map<std::string, std::map<unsigned int, std::set<std::string>>> _read_stop_per_position,
		std::map<std::string, std::map<std::string, std::map<unsigned int, std::string>>> _read_genotypes_per_position,
		std::map<std::string, std::map<unsigned int, std::set<std::string>>> _activeAlleles_per_position,
		std::map<std::string, std::map<std::string, std::string>> _MSA_reference_sequences,
		std::map<std::string, std::map<std::string, std::string>> _MSA_reference_sequences_whichHap
	);

	void makeInference(std::string geneID);

	std::vector<std::vector<HMMstate>> statesByLevel;

};

#endif /* FULLLENGTHHMM_FULLLENGTHHMM_H_ */
