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
#include <fstream>

struct HMMstate
{
	unsigned int level;
	std::pair<std::string, std::string> haplotypes_alleles;
	std::pair<unsigned int, unsigned int> copyingFrom;
	size_t readAssignmentState;

	double fw_p;
	double bw_p;
	double viterbi_p;
	size_t viterbi_p_whereFrom;
};

struct HMMtransition
{
	size_t from_level;
	size_t from_state;
	size_t to_state;
	double P;
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
	std::map<std::string, std::set<std::pair<std::string, unsigned char>>> readAssignment_2_activeReads;

	std::vector<std::map<size_t, std::set<size_t>>> level_readAssignmentState_2_states;
	std::map<std::string, std::size_t> readID_2_index;
	std::vector<std::string> readIndex_2_ID;


	double tr_change_h1h2;
	double tr_noChange_h1h2;

	double tr_within_h1h2_changeTemplate;
	double tr_within_h1h2_NoChangeTemplate;

	double tr_change_oneh1;
	double tr_change_oneh2;
	double tr_change_oneh;

	// no A1 / A2 allele group information
	double tr_noHaplotypeResolution_remain_oneh;
	double tr_noHaplotypeResolution_change_oneh;

	// with A1 / A2 allele group information
	double tr_haplotypeResolution_remain_h1;
	double tr_haplotypeResolution_change_h1_withinGroup;
	double tr_haplotypeResolution_change_h1_outOfGroup;

	double tr_haplotypeResolution_remain_h2;
	double tr_haplotypeResolution_change_h2_withinGroup;
	double tr_haplotypeResolution_change_h2_outOfGroup;


	std::vector<std::string> computeReadAssignmentSets(const std::set<std::string>& runningReadIDs) const;
	std::set<std::string> nextLevel_compatibleReadAssignments(const std::string& thisReadAssignment, const std::set<size_t>& disappearingReadIDs_nextLevel, const std::set<size_t>& newReadIDs_nextLevel) const;
	std::set<std::string> previousLevel_compatibleReadAssignments(const std::string& thisReadAssignment, const std::set<size_t>& newReadIDs_previousLevel, const std::set<size_t>& disappearingReadIDs_thisLevel) const;

	std::string readAssignmentTemplate;

	std::vector<double> computeInitialProbabilities() const;
	std::vector<HMMtransition> computeLevelTransitions(size_t first_level, long long limitToFromState = -1) const;
	std::vector<HMMtransition> computeLevelTransitions_backward(size_t next_level, long long limitToToState = -1) const;

	std::vector<double> computeEmissionProbabilities(size_t level) const;

	std::vector<double> initialProbabilities;
	unsigned int currentGene_geneLength;
	std::string currentGene;
	std::set<std::string> currentGene_MSA_ids;
	std::map<std::string, unsigned int> currentGene_MSA_id_2_int;
	std::vector<std::string> currentGene_MSA_int_2_id;
	std::vector<std::string> currentGene_MSA_int_2_HaplotypeID;
	std::set<std::pair<unsigned int, unsigned int>> currentGene_MSA_ids_same_haploGroup;
	std::set<unsigned int> currentGene_MSA_ids_h1;
	std::set<unsigned int> currentGene_MSA_ids_h2;
	std::vector<std::vector<std::string>> currentGene_activeAlleles_perPosition;
	std::vector<std::set<char>> currentGene_activeAlleles_perPosition_allelesLength1_chars;
	bool currentGene_haplotypeResolved;

	double currentGene_novel_allele_p;


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

	void makeInference(std::string geneID, std::ofstream& output_fasta, std::string outputPrefix_furtherOutput);
	void makeInference(std::string geneID, std::ofstream& output_fasta, std::ofstream& output_graphLevels, std::string outputPrefix_furtherOutput);

	std::vector<std::vector<HMMstate>> statesByLevel;

};

#endif /* FULLLENGTHHMM_FULLLENGTHHMM_H_ */
