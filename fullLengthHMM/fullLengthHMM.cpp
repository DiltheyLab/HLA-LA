/*
 * fullLengthHMM.cpp
 *
 *  Created on: 24.01.2022
 *      Author: Alexa
 */

#include "fullLengthHMM.h"

#include <assert.h>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

#include <cmath>

fullLengthHMM::fullLengthHMM(
	std::map<std::string, unsigned int> _gene_length,
	std::map<std::string, std::string> _reads_2_genes,
	std::map<std::string, std::map<std::string, std::pair<unsigned int, unsigned int>>> _read_start_stop_positions,
	std::map<std::string, std::map<unsigned int, std::set<std::string>>> _read_start_per_position,
	std::map<std::string, std::map<unsigned int, std::set<std::string>>> _read_stop_per_position,
	std::map<std::string, std::map<std::string, std::map<unsigned int, std::string>>> _read_genotypes_per_position,
	std::map<std::string, std::map<unsigned int, std::set<std::string>>> _activeAlleles_per_position,
	std::map<std::string, std::map<std::string, std::string>> _MSA_reference_sequences,
	std::map<std::string, std::map<std::string, std::string>> _MSA_reference_sequences_whichHap
) :
	gene_length(_gene_length),
	reads_2_genes(_reads_2_genes),
	read_start_stop_positions(_read_start_stop_positions),
	read_start_per_position(_read_start_per_position),
	read_stop_per_position(_read_stop_per_position),
	read_genotypes_per_position(_read_genotypes_per_position),
	activeAlleles_per_position(_activeAlleles_per_position),
	MSA_reference_sequences(_MSA_reference_sequences),
	MSA_reference_sequences_whichHap(_MSA_reference_sequences_whichHap)
{

}

std::set<std::string> fullLengthHMM::nextLevel_compatibleReadAssignments(const std::string& thisReadAssignment, const std::set<size_t>& disappearingReadIDs_nextLevel, const std::set<size_t>& newReadIDs_nextLevel)
{
	std::vector<std::string> readAssignmentStates = {thisReadAssignment};
	for(auto disappearingReadIDIdx : disappearingReadIDs_nextLevel)
	{
		readAssignmentStates.at(0).at(disappearingReadIDIdx) = 'N';
	}

	readAssignmentStates.resize(std::pow(2, newReadIDs_nextLevel.size()), readAssignmentTemplate);

	size_t vectorElements_filled = 1;
	for(auto newReadIdIdx : newReadIDs_nextLevel)
	{
		for(unsigned int existingElementI = 0; existingElementI < vectorElements_filled; existingElementI++)
		{
			unsigned int copyIntoElementI = vectorElements_filled + existingElementI;
			readAssignmentStates.at(copyIntoElementI) = readAssignmentStates.at(existingElementI);
			readAssignmentStates.at(existingElementI).at(newReadIdIdx) = '1';
			readAssignmentStates.at(copyIntoElementI).at(newReadIdIdx) = '2';
		}
		vectorElements_filled = 2 * vectorElements_filled;
	}

	return std::set<std::string>(readAssignmentStates.begin(), readAssignmentStates.end());
}

std::vector<std::string> fullLengthHMM::computeReadAssignmentSets(const std::set<std::string>& runningReadIDs)
{
	assert(readAssignmentTemplate.size());
	std::vector<std::string> readAssignmentStates;

	if(runningReadIDs.size() == 0)
	{
		readAssignmentStates.push_back(readAssignmentTemplate);
	}
	else
	{
		readAssignmentStates.resize(std::pow(2, runningReadIDs.size()), readAssignmentTemplate);

		size_t vectorElements_filled = 1;
		for(auto readID : runningReadIDs)
		{
			unsigned int readID_index = readID_2_index.at(readID);
			for(unsigned int existingElementI = 0; existingElementI < vectorElements_filled; existingElementI++)
			{
				unsigned int copyIntoElementI = vectorElements_filled + existingElementI;
				readAssignmentStates.at(copyIntoElementI) = readAssignmentStates.at(existingElementI);
				readAssignmentStates.at(existingElementI).at(readID_index) = '1';
				readAssignmentStates.at(copyIntoElementI).at(readID_index) = '2';
			}
			vectorElements_filled = 2 * vectorElements_filled;

		}
		assert(readAssignmentStates.size() == vectorElements_filled);

		if(0 && (runningReadIDs.size() == 3))
		{
			std::cerr << "Read assignment vectors:\n";
			for(unsigned int i = 0; i < readAssignmentStates.size(); i++)
			{
				std::cerr << " " << i << " " << readAssignmentStates.at(i) << "\n";
			}
			std::cerr << std::flush;
			assert(1 == 0);
		}
	}

	return readAssignmentStates;
}

void fullLengthHMM::makeInference(std::string geneID)
{
	assert(gene_length.count(geneID));
	std::cout << "Now making HMM-based inference for gene " << geneID << "\n" << std::flush;


	bool haplotypeResolved = true;
	unsigned int MSA_h1_n = 0;
	unsigned int MSA_h2_n = 0;
	unsigned int geneLength = gene_length.at(geneID);
	std::set<std::string> MSA_ids;
	std::map<std::string, unsigned int> MSA_id_2_int;
	std::vector<std::string> MSA_int_2_id;
	for(auto MSAentry : MSA_reference_sequences.at(geneID))
	{
		MSA_ids.insert(MSAentry.first);
		MSA_id_2_int[MSAentry.first] = MSA_ids.size() - 1;
		MSA_int_2_id.push_back(MSAentry.first);
		assert((MSA_reference_sequences_whichHap.at(geneID).at(MSAentry.first) == "?") || (MSA_reference_sequences_whichHap.at(geneID).at(MSAentry.first) == "1") || (MSA_reference_sequences_whichHap.at(geneID).at(MSAentry.first) == "2"));
		if(MSA_reference_sequences_whichHap.at(geneID).at(MSAentry.first) == "?")
			haplotypeResolved = false;

		if(MSA_reference_sequences_whichHap.at(geneID).at(MSAentry.first) == "1")
			MSA_h1_n++;
		else
			MSA_h2_n++;
	}

	std::set<std::pair<unsigned int, unsigned int>> MSA_ids_same_haploGroup;
	std::set<unsigned int> MSA_ids_h1;

	for(auto MSAentry_1 : MSA_reference_sequences.at(geneID))
	{
		unsigned int MSAentry_1_id = MSA_id_2_int.at(MSAentry_1.first);
		if(MSA_reference_sequences_whichHap.at(geneID).at(MSAentry_1.first) == "1")
		{
			MSA_ids_h1.insert(MSAentry_1_id);
		}
		for(auto MSAentry_2 : MSA_reference_sequences.at(geneID))
		{
			unsigned int MSAentry_2_id = MSA_id_2_int.at(MSAentry_2.first);
			if(haplotypeResolved)
			{
				if(MSA_reference_sequences_whichHap.at(geneID).at(MSAentry_1.first) == MSA_reference_sequences_whichHap.at(geneID).at(MSAentry_2.first))
				{
					MSA_ids_same_haploGroup.insert(std::make_pair(MSAentry_1_id, MSAentry_2_id));
					MSA_ids_same_haploGroup.insert(std::make_pair(MSAentry_2_id, MSAentry_1_id));
				}
			}
		}
	}

	std::map<std::vector<bool>, std::string> myTest;

	std::vector<std::vector<std::string>> activeAlleles_perPosition;
	activeAlleles_perPosition.resize(geneLength);
	for(unsigned int i = 0; i < geneLength; i++)
	{
		std::vector<std::string> activeAlleles_vec(activeAlleles_per_position.at(geneID).at(i).begin(), activeAlleles_per_position.at(geneID).at(i).end());
		activeAlleles_perPosition.at(i) = activeAlleles_vec;
	}

	statesByLevel.clear();
	statesByLevel.resize(geneLength);

	level_readAssignmentState_2_states.clear();
	level_readAssignmentState_2_states.resize(geneLength);

	std::vector<long long> states_per_position;
	states_per_position.resize(geneLength, 0);

	std::set<std::string> readIDs;
	for(auto read2geneEntry : reads_2_genes)
	{
		if(read2geneEntry.second == geneID)
		{
			readIDs.insert(read2geneEntry.first);
		}
	}

	readID_2_index.clear();
	std::vector<std::string> readIDs_vector(readIDs.begin(), readIDs.end());
	for(unsigned int i = 0; i < readIDs_vector.size(); i++)
	{
		readID_2_index[readIDs_vector.at(i)] = i;
	}

	readAssignmentTemplate.clear();
	readAssignmentTemplate.resize(readIDs_vector.size(), 'N');
	assert(readAssignmentTemplate.size() == readIDs_vector.size());

	readAssignmentStates.clear();
	readAssignmentState_2_index.clear();

	std::cout << "Gene " << geneID << ", starting inference with " << readIDs_vector.size() << " reads.\n" << std::flush;

	std::vector<unsigned int> readAssingmentStates;
	std::set<std::string> runningReadIDs;
	size_t n_states = 0;
	int recompute_readAssingmentStates = 0;
	for(unsigned int levelI = 0; levelI < geneLength; levelI++)
	{
		std::cerr << "Round I: Level " << levelI << " / " << geneLength << "\n" << std::flush;

		if(read_start_per_position.at(geneID).count(levelI))
		{
			for(auto readID : read_start_per_position.at(geneID).at(levelI))
			{
				runningReadIDs.insert(readID);
				recompute_readAssingmentStates++;
			}
		}

		std::vector<std::string> possibleReadAssignmentStates_thisLevel = computeReadAssignmentSets(runningReadIDs);
		for(size_t readAssignmentStateI = 0; readAssignmentStateI < possibleReadAssignmentStates_thisLevel.size(); readAssignmentStateI++)
		{
			const std::string& thisReadAssignmentState = possibleReadAssignmentStates_thisLevel.at(readAssignmentStateI);
			if(readAssignmentState_2_index.count(thisReadAssignmentState) == 0)
			{
				readAssignmentStates.push_back(thisReadAssignmentState);
				readAssignmentState_2_index[thisReadAssignmentState] = readAssignmentStates.size() - 1;

			}
		}

		assert(activeAlleles_per_position.at(geneID).count(levelI));
		assert(activeAlleles_per_position.at(geneID).at(levelI).size() >= 1);

		size_t expectedStates_thisLevel =
					((MSA_int_2_id.size() * MSA_int_2_id.size() - MSA_int_2_id.size()) / 2.0 + MSA_int_2_id.size()) *
					possibleReadAssignmentStates_thisLevel.size() *
					(activeAlleles_perPosition.at(levelI).size() * activeAlleles_perPosition.at(levelI).size());
		statesByLevel.at(levelI).reserve(expectedStates_thisLevel);

		for(unsigned int copyFromI_1 = 0; copyFromI_1 < ((long long)MSA_int_2_id.size()-1); copyFromI_1++)
		{
			for(unsigned int copyFromI_2 = copyFromI_1; copyFromI_2 < MSA_int_2_id.size(); copyFromI_2++)
			{
				for(unsigned int readAssignmentStateI = 0; readAssignmentStateI < possibleReadAssignmentStates_thisLevel.size(); readAssignmentStateI++)
				{
					for(unsigned int h1_alleleIndex = 0; h1_alleleIndex < activeAlleles_perPosition.at(levelI).size(); h1_alleleIndex++)
					{
						for(unsigned int h2_alleleIndex = 0; h2_alleleIndex < activeAlleles_perPosition.at(levelI).size(); h2_alleleIndex++)
						{
							const std::string& thisReadAssignmentState = possibleReadAssignmentStates_thisLevel.at(readAssignmentStateI);
							size_t thisReadAssignmentState_index = readAssignmentState_2_index.at(thisReadAssignmentState);

							HMMstate s;
							s.copyingFrom = std::make_pair(copyFromI_1, copyFromI_2);
							s.haplotypes_alleles = std::make_pair(activeAlleles_perPosition.at(levelI).at(h1_alleleIndex), activeAlleles_perPosition.at(levelI).at(h2_alleleIndex));
							s.level = levelI;
							s.readAssignmentState = thisReadAssignmentState_index;
							statesByLevel.at(levelI).push_back(s);
							n_states++;

							size_t s_index = statesByLevel.at(levelI).size() - 1;
							level_readAssignmentState_2_states.at(levelI)[thisReadAssignmentState_index].insert(s_index);
						}
					}
				}
			}
		}

		if(read_stop_per_position.at(geneID).count(levelI))
		{
			for(auto readID : read_stop_per_position.at(geneID).at(levelI))
			{
				assert(runningReadIDs.count(readID));
				runningReadIDs.erase(readID);
				recompute_readAssingmentStates++;
			}
		}
	}


	/*
	 * HMM state transition structure
	 * - The basic scaffold is formed by {reads} x {copied-from haplotypes}
	 * - Within this scaffold, we have a higher probability to transition into states for which the underlying allele is
	 *   identical to the copied-from allele
	 * - Within the copied-from haplotypes, we transition with a particular probability that is higher for haplotypes within the same group
	 * - Within the set of read sets, we transition into compatible read sets
	 * - Should we have a forward lookup function?
	 *   Conditional on being in state s at level i, give all states in state s + 1 with their corresponding transition probabilites
	 *   --> We can freely transition within the set of copied-from haplotypes, and we can say we leave the current haplotype set with
	 *       probability 0.1 * 1/length(gene) (i.e. only 0.1 expected haplotype recombination events). Furthermore, within the same set of
	 *       haplotypes, leave the current haplotype with probability 1/length(gene). If we leave the haplotype, we have uniform transition probabilities.
	 *   --> For each read config / copied-from haplotype configuration, we prefer states in which the copied-from allele and the true allele
	 *       agree with probability 1 - 3/length(gene) - i.e. on average 3 novel positions per gene. This could be made adaptable.
	 *   --> What happens when we transition between read configurations?
	 *       First, we determine the set of reads that disappear after the current position. These are ignored.
	 *       Second, we could determine the set of reads that appear at the next position. We can build arbitrary configurations for these.
	 *       Apart from that, we just copy the current read configuration.
	 *       Perhaps just use a fixed-length string to store the haplotype assignments of reads?
	 */


	double change_h1h2 = 0.1 * (1.0 / geneLength);
	double noChange_h1h2 = 1 - change_h1h2;

	double within_h1h2_changeTemplate = 1.0 / geneLength;
	double within_h1h2_NoChangeTemplate = 1 - within_h1h2_changeTemplate;

	double change_oneh1 = (haplotypeResolved) ? (1.0 / MSA_h1_n) : -1;
	double change_oneh2 = (haplotypeResolved) ? (1.0 / MSA_h2_n) : -1;
	double change_oneh = 1.0/MSA_reference_sequences.size();

	double noHaplotypeResolution_remain_oneh = within_h1h2_NoChangeTemplate + within_h1h2_changeTemplate * change_oneh;
	double noHaplotypeResolution_change_oneh = within_h1h2_changeTemplate * change_oneh;

	double haplotypeResolution_remain_h1 = noChange_h1h2 * within_h1h2_NoChangeTemplate + noChange_h1h2 * within_h1h2_changeTemplate * change_oneh1;
	double haplotypeResolution_change_h1_withinGroup = noChange_h1h2 * within_h1h2_changeTemplate * change_oneh1;
	double haplotypeResolution_change_h1_outOfGroup = change_h1h2 * change_oneh2;

	double haplotypeResolution_remain_h2 = noChange_h1h2 * within_h1h2_NoChangeTemplate + noChange_h1h2 * within_h1h2_changeTemplate * change_oneh2;
	double haplotypeResolution_change_h2_withinGroup = noChange_h1h2 * within_h1h2_changeTemplate * change_oneh2;
	double haplotypeResolution_change_h2_outOfGroup = change_h1h2 * change_oneh1;

	double novel_allele_p = 3.0/geneLength;

	size_t n_jumps = 0;
	for(unsigned int levelI = 0; levelI < (geneLength - 1); levelI++)
	{
		std::cerr << "Round II: Level " << levelI << " / " << geneLength << "\n" << std::flush;

		std::set<size_t> diseappearing_read_IDs_indices;
		if(read_stop_per_position.at(geneID).count(levelI))
		{
			for(auto readID : read_stop_per_position.at(geneID).at(levelI))
			{
				size_t readID_index = readID_2_index.at(readID);
				diseappearing_read_IDs_indices.insert(readID_index);
			}
		}

		std::set<size_t> new_read_IDs_indices_nextLevel;
		if(read_start_per_position.at(geneID).count(levelI+1))
		{
			for(auto readID : read_start_per_position.at(geneID).at(levelI+1))
			{
				size_t readID_index = readID_2_index.at(readID);
				new_read_IDs_indices_nextLevel.insert(readID_index);
			}
		}


		for(size_t stateI = 0; stateI < statesByLevel.at(levelI).size(); stateI++)
		{
			const HMMstate& s = statesByLevel.at(levelI).at(stateI);

			std::set<std::string> forward_read_assignment_states = nextLevel_compatibleReadAssignments(readAssignmentStates.at(s.readAssignmentState), diseappearing_read_IDs_indices, new_read_IDs_indices_nextLevel);

			double _sum_jump_Ps = 0;

			double read_assignment_state_p = 1.0 / forward_read_assignment_states.size();
			unsigned int n_activeAlleles_nextLevel = activeAlleles_perPosition.at(levelI+1).size();

			for(std::string readAssignmentState_nextLevel : forward_read_assignment_states)
			{
				size_t readAssignmentState_nextLevel_index = readAssignmentState_2_index.at(readAssignmentState_nextLevel);
				if(level_readAssignmentState_2_states.at(levelI+1).count(readAssignmentState_nextLevel_index) == 0)
				{
					std::cerr << "Missing a read assingment state at the next level.\n";
					std::cerr << "\t" << "levelI" << ": " << levelI << "\n";
					std::cerr << "\t" << "levelI+1" << ": " << levelI+1 << "\n";
					std::cerr << "\t" << "readAssignmentState_nextLevel_index" << ": " << readAssignmentState_nextLevel_index << "\n";
					std::cerr << "\t" << "readAssignmentState_nextLevel_index" << ": " << readAssignmentState_nextLevel_index << "\n";
					std::cerr << std::flush;
				}
				std::set<size_t> nextLevel_jumpIntoStates = level_readAssignmentState_2_states.at(levelI+1).at(readAssignmentState_nextLevel_index);

				for(auto nextStateIndex : nextLevel_jumpIntoStates)
				{
					const HMMstate& next_s = statesByLevel.at(levelI+1).at(nextStateIndex);

					double haplotype_copy_p = 1;

					if(haplotypeResolved)
					{
						if(s.copyingFrom.first != next_s.copyingFrom.first)
						{
							if(MSA_ids_same_haploGroup.count(std::make_pair(s.copyingFrom.first, next_s.copyingFrom.first)))
							{
								haplotype_copy_p *= (MSA_ids_h1.count(s.copyingFrom.first) ? haplotypeResolution_change_h1_withinGroup : haplotypeResolution_change_h2_withinGroup);
							}
							else
							{
								haplotype_copy_p *= (MSA_ids_h1.count(s.copyingFrom.first) ?  haplotypeResolution_change_h1_outOfGroup : haplotypeResolution_change_h2_outOfGroup);
							}
						}
						else
						{
							haplotype_copy_p *=  (MSA_ids_h1.count(s.copyingFrom.first) ? haplotypeResolution_remain_h1 : haplotypeResolution_remain_h2);
						}


						if(s.copyingFrom.second != next_s.copyingFrom.second)
						{
							if(MSA_ids_same_haploGroup.count(std::make_pair(s.copyingFrom.second, next_s.copyingFrom.second)))
							{
								haplotype_copy_p *= (MSA_ids_h1.count(s.copyingFrom.second) ? haplotypeResolution_change_h1_withinGroup : haplotypeResolution_change_h2_withinGroup);
							}
							else
							{
								haplotype_copy_p *= (MSA_ids_h1.count(s.copyingFrom.second) ?  haplotypeResolution_change_h1_outOfGroup : haplotypeResolution_change_h2_outOfGroup);
							}
						}
						else
						{
							haplotype_copy_p *=  (MSA_ids_h1.count(s.copyingFrom.second) ? haplotypeResolution_remain_h1 : haplotypeResolution_remain_h2);
						}
					}
					else
					{
						bool isSwitch_h1 = (s.copyingFrom.first != next_s.copyingFrom.first);
						bool isSwitch_h2 = (s.copyingFrom.second != next_s.copyingFrom.second);
						haplotype_copy_p =
									(isSwitch_h1 ? noHaplotypeResolution_change_oneh : noHaplotypeResolution_remain_oneh) *
									(isSwitch_h2 ? noHaplotypeResolution_change_oneh : noHaplotypeResolution_remain_oneh);
					}

					bool alleles_match_h1 = (next_s.haplotypes_alleles.first.length() == 1) && (next_s.haplotypes_alleles.first.at(0) == MSA_reference_sequences.at(geneID).at(MSA_int_2_id.at(next_s.copyingFrom.first)).at(levelI));
					bool alleles_match_h2 = (next_s.haplotypes_alleles.second.length() == 1) && (next_s.haplotypes_alleles.second.at(0) == MSA_reference_sequences.at(geneID).at(MSA_int_2_id.at(next_s.copyingFrom.second)).at(levelI));

					double allele_p =
							(alleles_match_h1 ? ((1 - novel_allele_p) + (novel_allele_p * (1.0 / n_activeAlleles_nextLevel))) : (novel_allele_p * (1.0 / n_activeAlleles_nextLevel))) *
							(alleles_match_h2 ? ((1 - novel_allele_p) + (novel_allele_p * (1.0 / n_activeAlleles_nextLevel))) : (novel_allele_p * (1.0 / n_activeAlleles_nextLevel)));

					double p_jump = read_assignment_state_p * haplotype_copy_p * allele_p;
					_sum_jump_Ps += p_jump;
					n_jumps++;
				}
			}

			std::cerr << "Level " << levelI << " state " << stateI << " sum of transition Ps: " << _sum_jump_Ps << "\n";
		}
	}

		/*
		unsigned int activeAlleles = activeAlleles_per_position.at(geneID).at(i).size();
		long long n_states =
				std::pow(MSA_ids.size(), 2) *
				std::pow(2, runningReadIDs.size()) *
				std::pow(activeAlleles, 2);

		states_per_position.at(i) = n_states;

		if(n_states > 1e6)
			std::cout << "\t" << i << "\t" << n_states << " [" << MSA_ids.size() << " " << runningReadIDs.size() << " " << activeAlleles << "]" <<  "\n";

		*/



	std::cout << geneID << " done - " << n_states << " states -- " << n_jumps << " jumgs.\n" << std::flush;

	// auto states_min_max = std::minmax_element(states_per_position.begin(), states_per_position.end());
	// std::cout << "=====" << "\n" << "\tMin: " << *(states_min_max.first) << " - max: " << *(states_min_max.second) << "\n" << std::flush;

}
