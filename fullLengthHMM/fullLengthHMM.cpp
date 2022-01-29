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
#include <omp.h>

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

std::set<std::string> fullLengthHMM::nextLevel_compatibleReadAssignments(const std::string& thisReadAssignment, const std::set<size_t>& disappearingReadIDs_nextLevel, const std::set<size_t>& newReadIDs_nextLevel) const
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

std::set<std::string> fullLengthHMM::previousLevel_compatibleReadAssignments(const std::string& thisReadAssignment, const std::set<size_t>& newReadIDs_previousLevel, const std::set<size_t>& disappearingReadIDs_thisLevel) const
{
	std::vector<std::string> readAssignmentStates = {thisReadAssignment};
	for(auto disappearingReadIDIdx : disappearingReadIDs_thisLevel)
	{
		readAssignmentStates.at(0).at(disappearingReadIDIdx) = 'N';
	}

	readAssignmentStates.resize(std::pow(2, newReadIDs_previousLevel.size()), readAssignmentTemplate);

	size_t vectorElements_filled = 1;
	for(auto newReadIdIdx : newReadIDs_previousLevel)
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

std::vector<std::string> fullLengthHMM::computeReadAssignmentSets(const std::set<std::string>& runningReadIDs) const
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
	omp_set_num_threads(32);
	
	assert(gene_length.count(geneID));
	std::cout << "Now making HMM-based inference for gene " << currentGene << "\n" << std::flush;

	bool currentGene_haplotypeResolved = true;
	unsigned int MSA_h1_n = 0;
	unsigned int MSA_h2_n = 0;

	currentGene_geneLength = gene_length.at(geneID);
	currentGene = geneID;
	currentGene_novel_allele_p = 3.0/currentGene_geneLength;

	currentGene_MSA_ids.clear();
	currentGene_MSA_id_2_int.clear();
	currentGene_MSA_int_2_id.clear();

	for(auto MSAentry : MSA_reference_sequences.at(currentGene))
	{
		currentGene_MSA_ids.insert(MSAentry.first);
		currentGene_MSA_id_2_int[MSAentry.first] = currentGene_MSA_ids.size() - 1;
		currentGene_MSA_int_2_id.push_back(MSAentry.first);
		assert((MSA_reference_sequences_whichHap.at(currentGene).at(MSAentry.first) == "?") || (MSA_reference_sequences_whichHap.at(currentGene).at(MSAentry.first) == "1") || (MSA_reference_sequences_whichHap.at(currentGene).at(MSAentry.first) == "2"));
		if(MSA_reference_sequences_whichHap.at(currentGene).at(MSAentry.first) == "?")
			currentGene_haplotypeResolved = false;

		if(MSA_reference_sequences_whichHap.at(currentGene).at(MSAentry.first) == "1")
			MSA_h1_n++;
		else
			MSA_h2_n++;
	}

	currentGene_MSA_ids_same_haploGroup.clear();
	currentGene_MSA_ids_h1.clear();
	currentGene_MSA_ids_h2.clear();

	for(auto MSAentry_1 : MSA_reference_sequences.at(currentGene))
	{
		unsigned int MSAentry_1_id = currentGene_MSA_id_2_int.at(MSAentry_1.first);
		if(MSA_reference_sequences_whichHap.at(currentGene).at(MSAentry_1.first) == "1")
		{
			currentGene_MSA_ids_h1.insert(MSAentry_1_id);
		}
		if(MSA_reference_sequences_whichHap.at(currentGene).at(MSAentry_1.first) == "2")
		{
			currentGene_MSA_ids_h2.insert(MSAentry_1_id);
		}
		for(auto MSAentry_2 : MSA_reference_sequences.at(currentGene))
		{
			unsigned int MSAentry_2_id = currentGene_MSA_id_2_int.at(MSAentry_2.first);
			if(currentGene_haplotypeResolved)
			{
				if(MSA_reference_sequences_whichHap.at(currentGene).at(MSAentry_1.first) == MSA_reference_sequences_whichHap.at(currentGene).at(MSAentry_2.first))
				{
					currentGene_MSA_ids_same_haploGroup.insert(std::make_pair(MSAentry_1_id, MSAentry_2_id));
					currentGene_MSA_ids_same_haploGroup.insert(std::make_pair(MSAentry_2_id, MSAentry_1_id));
				}
			}
		}
	}

	std::map<std::vector<bool>, std::string> myTest;

	currentGene_activeAlleles_perPosition.clear();
	currentGene_activeAlleles_perPosition.resize(currentGene_geneLength);
	for(unsigned int i = 0; i < currentGene_geneLength; i++)
	{
		std::vector<std::string> activeAlleles_vec(activeAlleles_per_position.at(currentGene).at(i).begin(), activeAlleles_per_position.at(currentGene).at(i).end());
		currentGene_activeAlleles_perPosition.at(i) = activeAlleles_vec;
	}

	statesByLevel.clear();
	statesByLevel.resize(currentGene_geneLength);

	level_readAssignmentState_2_states.clear();
	level_readAssignmentState_2_states.resize(currentGene_geneLength);


	std::vector<long long> states_per_position;
	states_per_position.resize(currentGene_geneLength, 0);

	std::set<std::string> readIDs;
	for(auto read2geneEntry : reads_2_genes)
	{
		if(read2geneEntry.second == currentGene)
		{
			readIDs.insert(read2geneEntry.first);
		}
	}

	readID_2_index.clear();
	readIndex_2_ID.clear();
	std::vector<std::string> readIDs_vector(readIDs.begin(), readIDs.end());
	readIndex_2_ID = readIDs_vector;
	for(unsigned int i = 0; i < readIDs_vector.size(); i++)
	{
		readID_2_index[readIDs_vector.at(i)] = i;
	}

	readAssignmentTemplate.clear();
	readAssignmentTemplate.resize(readIDs_vector.size(), 'N');
	assert(readAssignmentTemplate.size() == readIDs_vector.size());

	readAssignmentStates.clear();
	readAssignmentState_2_index.clear();

	std::cout << "Gene " << currentGene << ", starting inference with " << readIDs_vector.size() << " reads.\n" << std::flush;

	/*
	{
		for(unsigned int copyFromI_1 = 0; copyFromI_1 < (long long)currentGene_MSA_int_2_id.size(); copyFromI_1++)
		{
			HMMstate s;
			s.copyingFrom = std::make_pair(copyFromI_1, 0);

			std::cout << "State " << copyFromI_1 << "\n";
			double outgoing_p = 0;
			for(unsigned int copyFromI_3 = 0; copyFromI_3 < ((long long)currentGene_MSA_int_2_id.size()); copyFromI_3++)
			{
				HMMstate next_s;
				next_s.copyingFrom = std::make_pair(copyFromI_3, 0);

				double haplotype_copy_p = 1;

				double change_h1h2 = 0.1 * (1.0 / currentGene_geneLength);
				double noChange_h1h2 = 1 - change_h1h2;

				double within_h1h2_changeTemplate = 1.0 / currentGene_geneLength;
				double within_h1h2_NoChangeTemplate = 1 - within_h1h2_changeTemplate;

				double change_oneh1 = (currentGene_haplotypeResolved) ? (1.0 / MSA_h1_n) : -1;
				double change_oneh2 = (currentGene_haplotypeResolved) ? (1.0 / MSA_h2_n) : -1;
				double change_oneh = 1.0/MSA_reference_sequences.size();

				double noHaplotypeResolution_remain_oneh = within_h1h2_NoChangeTemplate + within_h1h2_changeTemplate * change_oneh;
				double noHaplotypeResolution_change_oneh = within_h1h2_changeTemplate * change_oneh;

				double haplotypeResolution_remain_h1 = noChange_h1h2 * within_h1h2_NoChangeTemplate + noChange_h1h2 * within_h1h2_changeTemplate * change_oneh1;
				double haplotypeResolution_change_h1_withinGroup = noChange_h1h2 * within_h1h2_changeTemplate * change_oneh1;
				double haplotypeResolution_change_h1_outOfGroup = change_h1h2 * change_oneh2;

				double haplotypeResolution_remain_h2 = noChange_h1h2 * within_h1h2_NoChangeTemplate + noChange_h1h2 * within_h1h2_changeTemplate * change_oneh2;
				double haplotypeResolution_change_h2_withinGroup = noChange_h1h2 * within_h1h2_changeTemplate * change_oneh2;
				double haplotypeResolution_change_h2_outOfGroup = change_h1h2 * change_oneh1;

				if(s.copyingFrom.first != next_s.copyingFrom.first)
				{
					if(currentGene_MSA_ids_same_haploGroup.count(std::make_pair(s.copyingFrom.first, next_s.copyingFrom.first)))
					{
						haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ? haplotypeResolution_change_h1_withinGroup : haplotypeResolution_change_h2_withinGroup);
					}
					else
					{
						haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ?  haplotypeResolution_change_h1_outOfGroup : haplotypeResolution_change_h2_outOfGroup);
					}
				}
				else
				{
					haplotype_copy_p *=  (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ? haplotypeResolution_remain_h1 : haplotypeResolution_remain_h2);
				}

				std::cout << "\t\t ==> state " << copyFromI_3 << ": " << haplotype_copy_p << "\n";
				outgoing_p += haplotype_copy_p;
			}
			std::cout << "\t" << outgoing_p << "\n" << std::flush;

		}

		assert(2 == 3);
	}
	*/
	/*
	{

		for(unsigned int copyFromI_1 = 0; copyFromI_1 < (long long)currentGene_MSA_int_2_id.size(); copyFromI_1++)
		{
			for(unsigned int copyFromI_2 = 0; copyFromI_2 < currentGene_MSA_int_2_id.size(); copyFromI_2++)
			{
				HMMstate s;
				s.copyingFrom = std::make_pair(copyFromI_1, copyFromI_2);

				std::cout << "State " << copyFromI_1 << " / " << copyFromI_2 << "\n";

				double outgoing_p = 0;

				for(unsigned int copyFromI_3 = 0; copyFromI_3 < ((long long)currentGene_MSA_int_2_id.size()); copyFromI_3++)
				{
					for(unsigned int copyFromI_4 = 0; copyFromI_4 < currentGene_MSA_int_2_id.size(); copyFromI_4++)
					{
						HMMstate next_s;
						next_s.copyingFrom = std::make_pair(copyFromI_3, copyFromI_4);

						double haplotype_copy_p = 1;

						double change_h1h2 = 0.1 * (1.0 / currentGene_geneLength);
						double noChange_h1h2 = 1 - change_h1h2;

						double within_h1h2_changeTemplate = 1.0 / currentGene_geneLength;
						double within_h1h2_NoChangeTemplate = 1 - within_h1h2_changeTemplate;

						double change_oneh1 = (currentGene_haplotypeResolved) ? (1.0 / MSA_h1_n) : -1;
						double change_oneh2 = (currentGene_haplotypeResolved) ? (1.0 / MSA_h2_n) : -1;
						double change_oneh = 1.0/MSA_reference_sequences.size();

						double noHaplotypeResolution_remain_oneh = within_h1h2_NoChangeTemplate + within_h1h2_changeTemplate * change_oneh;
						double noHaplotypeResolution_change_oneh = within_h1h2_changeTemplate * change_oneh;

						double haplotypeResolution_remain_h1 = noChange_h1h2 * within_h1h2_NoChangeTemplate + noChange_h1h2 * within_h1h2_changeTemplate * change_oneh1;
						double haplotypeResolution_change_h1_withinGroup = noChange_h1h2 * within_h1h2_changeTemplate * change_oneh1;
						double haplotypeResolution_change_h1_outOfGroup = change_h1h2 * change_oneh2;

						double haplotypeResolution_remain_h2 = noChange_h1h2 * within_h1h2_NoChangeTemplate + noChange_h1h2 * within_h1h2_changeTemplate * change_oneh2;
						double haplotypeResolution_change_h2_withinGroup = noChange_h1h2 * within_h1h2_changeTemplate * change_oneh2;
						double haplotypeResolution_change_h2_outOfGroup = change_h1h2 * change_oneh1;

						if(s.copyingFrom.first != next_s.copyingFrom.first)
						{
							if(currentGene_MSA_ids_same_haploGroup.count(std::make_pair(s.copyingFrom.first, next_s.copyingFrom.first)))
							{
								haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ? haplotypeResolution_change_h1_withinGroup : haplotypeResolution_change_h2_withinGroup);
							}
							else
							{
								haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ?  haplotypeResolution_change_h1_outOfGroup : haplotypeResolution_change_h2_outOfGroup);
							}
						}
						else
						{
							haplotype_copy_p *=  (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ? haplotypeResolution_remain_h1 : haplotypeResolution_remain_h2);
						}


						if(s.copyingFrom.second != next_s.copyingFrom.second)
						{
							if(currentGene_MSA_ids_same_haploGroup.count(std::make_pair(s.copyingFrom.second, next_s.copyingFrom.second)))
							{
								haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.second) ? haplotypeResolution_change_h1_withinGroup : haplotypeResolution_change_h2_withinGroup);
							}
							else
							{
								haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.second) ?  haplotypeResolution_change_h1_outOfGroup : haplotypeResolution_change_h2_outOfGroup);
							}
						}
						else
						{
							haplotype_copy_p *=  (currentGene_MSA_ids_h1.count(s.copyingFrom.second) ? haplotypeResolution_remain_h1 : haplotypeResolution_remain_h2);
						}

						std::cout << "\t\t ==> state " << copyFromI_3 << " / " << copyFromI_4 << ": " << haplotype_copy_p << "\n";
						outgoing_p += haplotype_copy_p;
					}
				}

				std::cout << "\t" << outgoing_p << "\n" << std::flush;
			}
		}

		assert(1 == 0);
	}

	*/

	std::vector<unsigned int> readAssingmentStates;
	std::set<std::string> runningReadIDs;
	size_t n_states = 0;
	int recompute_readAssingmentStates = 0;
	for(unsigned int first_level = 0; first_level < currentGene_geneLength; first_level++)
	{
		if((first_level % 10) == 0)
			std::cerr << "Round I: Level " << first_level << " / " << currentGene_geneLength << "\n" << std::flush;

		if(read_start_per_position.at(currentGene).count(first_level))
		{
			for(auto readID : read_start_per_position.at(currentGene).at(first_level))
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

		assert(activeAlleles_per_position.at(currentGene).count(first_level));
		assert(activeAlleles_per_position.at(currentGene).at(first_level).size() >= 1);

		size_t expectedStates_thisLevel =
					((currentGene_MSA_int_2_id.size() * currentGene_MSA_int_2_id.size() - currentGene_MSA_int_2_id.size()) / 2.0 + currentGene_MSA_int_2_id.size()) *
					possibleReadAssignmentStates_thisLevel.size() *
					(currentGene_activeAlleles_perPosition.at(first_level).size() * currentGene_activeAlleles_perPosition.at(first_level).size());
		statesByLevel.at(first_level).reserve(expectedStates_thisLevel);

		for(unsigned int copyFromI_1 = 0; copyFromI_1 < currentGene_MSA_int_2_id.size(); copyFromI_1++)
		{
			for(unsigned int copyFromI_2 = 0; copyFromI_2 < currentGene_MSA_int_2_id.size(); copyFromI_2++)
			{
				for(unsigned int readAssignmentStateI = 0; readAssignmentStateI < possibleReadAssignmentStates_thisLevel.size(); readAssignmentStateI++)
				{
					const std::string& thisReadAssignmentState = possibleReadAssignmentStates_thisLevel.at(readAssignmentStateI);
					size_t thisReadAssignmentState_index = readAssignmentState_2_index.at(thisReadAssignmentState);

					for(unsigned int h1_alleleIndex = 0; h1_alleleIndex < currentGene_activeAlleles_perPosition.at(first_level).size(); h1_alleleIndex++)
					{
						for(unsigned int h2_alleleIndex = 0; h2_alleleIndex < currentGene_activeAlleles_perPosition.at(first_level).size(); h2_alleleIndex++)
						{
							HMMstate s;
							s.copyingFrom = std::make_pair(copyFromI_1, copyFromI_2);
							s.haplotypes_alleles = std::make_pair(currentGene_activeAlleles_perPosition.at(first_level).at(h1_alleleIndex), currentGene_activeAlleles_perPosition.at(first_level).at(h2_alleleIndex));
							s.level = first_level;
							s.readAssignmentState = thisReadAssignmentState_index;
							statesByLevel.at(first_level).push_back(s);
							n_states++;

							size_t s_index = statesByLevel.at(first_level).size() - 1;
							level_readAssignmentState_2_states.at(first_level)[thisReadAssignmentState_index].insert(s_index);
						}
					}
				}
			}
		}

		if(read_stop_per_position.at(currentGene).count(first_level))
		{
			for(auto readID : read_stop_per_position.at(currentGene).at(first_level))
			{
				assert(runningReadIDs.count(readID));
				runningReadIDs.erase(readID);
				recompute_readAssingmentStates++;
			}
		}
	}

	initialProbabilities = computeInitialProbabilities();

	readAssignment_2_activeReads.clear();
	readAssignment_2_activeReads[readAssignmentTemplate] = {};
	for(const std::string& readAssignmentString : readAssignmentStates)
	{
		for(unsigned int readI = 0; readI < readAssignmentString.size(); readI++)
		{
			if(readAssignmentString.at(readI) != 'N')
			{
				std::string readID = readIndex_2_ID.at(readI);
				assert((readAssignmentString.at(readI) == '1') || (readAssignmentString.at(readI) == '2'));
				readAssignment_2_activeReads[readAssignmentString].insert(std::make_pair(readID, readAssignmentString.at(readI)));
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


	size_t n_jumps = 0;
	for(unsigned int levelI = 0; levelI < currentGene_geneLength; levelI++)
	{
		if((levelI % 10) == 0)
			std::cerr << "Round II: Level " << levelI << " / " << currentGene_geneLength << "\n" << std::flush;

		std::vector<double> emissionP = computeEmissionProbabilities(levelI);

		std::vector<std::map<size_t, double>> states_levelI_jumpFrom;

		if(levelI > 0)
		{
			std::vector<HMMtransition> transitions_to_next_level = computeLevelTransitions(levelI-1);

			states_levelI_jumpFrom.resize(statesByLevel.at(levelI).size());
			for(auto j : transitions_to_next_level)
			{
				assert(j.from_level == (levelI - 1));
				assert(states_levelI_jumpFrom.at(j.to_state).count(j.from_state) == 0);
				states_levelI_jumpFrom.at(j.to_state)[j.from_state] = j.P;
			}
		}

		size_t states = statesByLevel.at(levelI).size();
		size_t chunk_size = states / omp_get_num_threads();
		#pragma omp parallel
		{
			size_t thisThread = omp_get_thread_num();
			size_t first_stateI = thisThread * chunk_size;
			size_t last_stateI = (thisThread+1) * chunk_size - 1;
			if((thisThread == omp_get_num_threads()) && (last_stateI < (states-1)))
			{
				last_stateI = states - 1;
			}
			assert(last_stateI <= (states - 1));

			for(size_t stateI = first_stateI; stateI <= last_stateI ; stateI++)
			{
				HMMstate& s = statesByLevel.at(levelI).at(stateI);
				if(levelI == 0)
				{
					s.fw_p = initialProbabilities.at(stateI) * emissionP.at(stateI);
					s.viterbi_p = s.fw_p;
				}
				else
				{
					double running_viterbi_max_p = -1;
					size_t running_viterbi_max_whereFrom;

					s.fw_p = 0;
					assert(states_levelI_jumpFrom.at(stateI).size());
					for(auto jumpIntoThisState : states_levelI_jumpFrom.at(stateI))
					{
						s.fw_p += statesByLevel.at(levelI-1).at(jumpIntoThisState.first).fw_p * jumpIntoThisState.second;
						assert(s.fw_p >= 0);

						double thisJump_viterbi = statesByLevel.at(levelI-1).at(jumpIntoThisState.first).viterbi_p * jumpIntoThisState.second;
						assert(thisJump_viterbi >= 0);
						if(thisJump_viterbi > running_viterbi_max_p)
						{
							running_viterbi_max_p = thisJump_viterbi;
							running_viterbi_max_whereFrom = jumpIntoThisState.first;
						}
					}

					s.fw_p *= emissionP.at(stateI);
					assert(s.fw_p >= 0);

					running_viterbi_max_p *= emissionP.at(stateI);

					s.viterbi_p = running_viterbi_max_p;
					s.viterbi_p_whereFrom = running_viterbi_max_whereFrom;
				}
			}
		}
	}

	double forward_joint_p = 0;
	for(const auto& s: statesByLevel.at(statesByLevel.size() - 1))
	{
		forward_joint_p += s.fw_p;
	}

	std::cout << "Joint forward P: " << forward_joint_p << "\n" << std::flush;


	for(long long levelI = (currentGene_geneLength - 1); levelI >= 0; levelI--)
	{
		if((levelI % 10) == 0)
			std::cerr << "Round III (Backward): Level " << levelI << "\n" << std::flush;

		std::vector<double> emissionP = computeEmissionProbabilities(levelI);

		if(levelI == (currentGene_geneLength - 1))
		{
			for(size_t stateI = 0; stateI < statesByLevel.at(levelI).size(); stateI++)
			{
				HMMstate& s = statesByLevel.at(levelI).at(stateI);
				s.bw_p = emissionP.at(stateI);
			}
		}
		else
		{
			std::vector<HMMtransition> transitions_to_next_level = computeLevelTransitions(levelI);

			std::vector<std::map<size_t, double>> states_levelI_jumpInto;
			states_levelI_jumpInto.resize(statesByLevel.at(levelI).size());
			for(auto j : transitions_to_next_level)
			{
				assert((long long)j.from_level == levelI);
				assert(states_levelI_jumpInto.at(j.from_state).count(j.to_state) == 0);
				states_levelI_jumpInto.at(j.from_state)[j.to_state] = j.P;
			}

			for(size_t stateI = 0; stateI < statesByLevel.at(levelI).size(); stateI++)
			{
				HMMstate& s = statesByLevel.at(levelI).at(stateI);

				double backward_sum = 0;
				for(auto jumpsInto : states_levelI_jumpInto.at(stateI))
				{
					size_t nextLevel_jumpTarget = jumpsInto.first;
					double jumpP = jumpsInto.second;
					double nextLevel_jumpTarget_bw = statesByLevel.at(levelI+1).at(nextLevel_jumpTarget).bw_p;

					backward_sum += nextLevel_jumpTarget_bw * jumpP;
				}

				s.bw_p = (backward_sum * emissionP.at(stateI));

				if(levelI == 0)
				{
					s.bw_p *= initialProbabilities.at(stateI);
				}
			}
		}
	}

	double backward_joint_p = 0;
	for(const auto& s: statesByLevel.at(0))
	{
		backward_joint_p += s.bw_p;
	}

	std::cout << "Joint forward P: " << forward_joint_p << "\n" << std::flush;
	std::cout << "Joint backward P: " << backward_joint_p << "\n" << std::flush;


		/*
		unsigned int activeAlleles = activeAlleles_per_position.at(currentGene).at(i).size();
		long long n_states =
				std::pow(currentGene_MSA_ids.size(), 2) *
				std::pow(2, runningReadIDs.size()) *
				std::pow(activeAlleles, 2);

		states_per_position.at(i) = n_states;

		if(n_states > 1e6)
			std::cout << "\t" << i << "\t" << n_states << " [" << currentGene_MSA_ids.size() << " " << runningReadIDs.size() << " " << activeAlleles << "]" <<  "\n";

		*/

	std::vector<size_t> backtrace_Viterbi;
	{
		backtrace_Viterbi.reserve(currentGene_geneLength);
		double viterbi_maxP = -1;
		size_t viterbi_maxP_whichState;
		for(size_t stateI = 0; stateI < statesByLevel.at(currentGene_geneLength - 1).size(); stateI++)
		{
			const HMMstate& s = statesByLevel.at(currentGene_geneLength - 1).at(stateI);
			if(s.viterbi_p > viterbi_maxP)
			{
				viterbi_maxP = s.viterbi_p;
				viterbi_maxP_whichState = stateI;
			}
		}
		assert(viterbi_maxP > 0);
		backtrace_Viterbi.push_back(viterbi_maxP_whichState);

		size_t backtrace_runningState = viterbi_maxP_whichState;
		for(unsigned int levelI = (currentGene_geneLength - 1); levelI > 0; levelI--)
		{
			backtrace_runningState = statesByLevel.at(levelI).at(backtrace_runningState).viterbi_p_whereFrom;
			backtrace_Viterbi.push_back(backtrace_runningState);
		}
		assert(backtrace_Viterbi.size() == currentGene_geneLength);
		std::reverse(backtrace_Viterbi.begin(), backtrace_Viterbi.end());
	}
	std::string bt_h1;
	std::string bt_h2;
	{
		bt_h1.reserve(currentGene_geneLength);
		bt_h2.reserve(currentGene_geneLength);

		for(unsigned int levelI = 0; levelI < (currentGene_geneLength - 1); levelI++)
		{
			const HMMstate& s = statesByLevel.at(levelI).at(backtrace_Viterbi.at(levelI));
			bt_h1.append(s.haplotypes_alleles.first);
			bt_h2.append(s.haplotypes_alleles.second);
		}

		//assert(bt_h1.size() == currentGene_geneLength);
		//assert(bt_h2.size() == currentGene_geneLength);
	}

	std::ofstream tempOutputStream;
	tempOutputStream.open("tempOutput.txt", std::ios::out);
	if(! tempOutputStream.is_open())
	{
		throw std::runtime_error("Cannot open file for writing");
	}

	tempOutputStream << geneID << "\n";
	tempOutputStream << bt_h1 << "\n";
	tempOutputStream << bt_h2 << "\n";

	std::cout << currentGene << " done - " << n_states << " states -- " << n_jumps << " jumgs.\n" << std::flush;

	// auto states_min_max = std::minmax_element(states_per_position.begin(), states_per_position.end());
	// std::cout << "=====" << "\n" << "\tMin: " << *(states_min_max.first) << " - max: " << *(states_min_max.second) << "\n" << std::flush;

}

std::vector<double> fullLengthHMM::computeEmissionProbabilities(size_t level) const
{
	std::map<size_t, double> readAssignment2P;
	std::vector<double> emission_P;
	emission_P.reserve(statesByLevel.at(level).size());
	for(size_t stateI = 0; stateI < statesByLevel.at(level).size(); stateI++)
	{
		const HMMstate& s = statesByLevel.at(level).at(stateI);
		if(readAssignment2P.count(s.readAssignmentState))
		{
			emission_P.push_back(readAssignment2P.at(s.readAssignmentState));
		}
		else
		{
			double emissionP = 1;
			const std::string& readAssignmentString = readAssignmentStates.at(s.readAssignmentState);
			if(readAssignment_2_activeReads.count(readAssignmentString) == 0)
			{
				std::cerr << "Missing info for string: " << readAssignmentString << "\n" << std::flush;
			}
			assert(readAssignment_2_activeReads.count(readAssignmentString));
			for(auto activeReads : readAssignment_2_activeReads.at(readAssignmentString))
			{
				assert((activeReads.second == '1') || (activeReads.second == '2'));
				assert(read_genotypes_per_position.at(currentGene).count(activeReads.first));

				if(read_genotypes_per_position.at(currentGene).at(activeReads.first).count(level))
				{
					std::string readAllele = read_genotypes_per_position.at(currentGene).at(activeReads.first).at(level);
					std::string underlyingAllele = (activeReads.second == '1') ? s.haplotypes_alleles.first : s.haplotypes_alleles.second;
					emissionP *= ((readAllele == underlyingAllele) ? 0.99 : 0.01);
				}
			}
			readAssignment2P[s.readAssignmentState] = emissionP;
			emission_P.push_back(emissionP);
		}
	}

	return emission_P;
}

std::vector<double> fullLengthHMM::computeInitialProbabilities() const
{
	std::vector<double> forReturn;

	size_t hom_states = currentGene_MSA_ids_h1.size() + currentGene_MSA_ids_h2.size();
	double het_p = 0.9;
	double hom_p = 1 - het_p;
	double hom_p_h1 = hom_p * 0.5 * (1.0 / currentGene_MSA_ids_h1.size());
	double hom_p_h2 = hom_p * 0.5 * (1.0 / currentGene_MSA_ids_h2.size());

	size_t readAssignmentStates_at_level_0 = level_readAssignmentState_2_states.at(0).size();
	unsigned int n_activeAlleles_thisLevel = currentGene_activeAlleles_perPosition.at(0).size();
	assert(n_activeAlleles_thisLevel > 0);

	forReturn.reserve(statesByLevel.at(0).size());
	double sum_p = 0;
	for(size_t stateI = 0; stateI < statesByLevel.at(0).size(); stateI++)
	{
		const HMMstate& s = statesByLevel.at(0).at(stateI);

		double copyFromP;
		if(s.copyingFrom.first == s.copyingFrom.second)
		{
			if(currentGene_MSA_ids_h1.count(s.copyingFrom.first))
				copyFromP = hom_p_h1;
			else
				copyFromP = hom_p_h2;
		}
		else
		{
			if(s.copyingFrom.first < s.copyingFrom.second)
			{
				copyFromP = 2.0 / (statesByLevel.at(0).size() - hom_states);
			}
			else
			{
				copyFromP = 0;
			}
		}

		bool alleles_match_h1 = (s.haplotypes_alleles.first.length() == 1) && (s.haplotypes_alleles.first.at(0) == MSA_reference_sequences.at(currentGene).at(currentGene_MSA_int_2_id.at(s.copyingFrom.first)).at(0));
		bool alleles_match_h2 = (s.haplotypes_alleles.second.length() == 1) && (s.haplotypes_alleles.second.at(0) == MSA_reference_sequences.at(currentGene).at(currentGene_MSA_int_2_id.at(s.copyingFrom.second)).at(0));

		double allele_p =
				(alleles_match_h1 ? ((1 - currentGene_novel_allele_p) + (currentGene_novel_allele_p * (1.0 / n_activeAlleles_thisLevel))) : (currentGene_novel_allele_p * (1.0 / n_activeAlleles_thisLevel))) *
				(alleles_match_h2 ? ((1 - currentGene_novel_allele_p) + (currentGene_novel_allele_p * (1.0 / n_activeAlleles_thisLevel))) : (currentGene_novel_allele_p * (1.0 / n_activeAlleles_thisLevel)));

		double thisState_initialP = copyFromP * (1.0 / readAssignmentStates_at_level_0) * allele_p;
		forReturn.push_back(thisState_initialP);

		sum_p += thisState_initialP;
	}
	assert(abs(1 - sum_p) <= 1e-4);

	assert(forReturn.size() == statesByLevel.at(0).size());
	return forReturn;
}

std::vector<HMMtransition> fullLengthHMM::computeLevelTransitions(size_t first_level) const
{
	assert(currentGene.length());
	std::vector<HMMtransition> forReturn;

	double change_h1h2 = 0.1 * (1.0 / currentGene_geneLength); // probability to recombine into the other A1/A2 allele group initiate a haplotype recombination event to the other A1/A2 group
	double noChange_h1h2 = 1 - change_h1h2; // probability not to recombine into the other A1/A2 allele group

	double within_h1h2_changeTemplate = 1.0 / currentGene_geneLength; // probability to initiate a recombination event within the same A1/A2 group, conditional on not initiating a recombination event to the other A1/A2 allele group
	double within_h1h2_NoChangeTemplate = 1 - within_h1h2_changeTemplate; // probability not to initiate a recombination event, conditional on not initiating a recombination event to the other A1/A2 allele group

	double change_oneh1 = (currentGene_haplotypeResolved) ? (1.0 / currentGene_MSA_ids_h1.size()) : -1; // probability to jump to a particular haplotype within the A1 group upon recombination within A1
	double change_oneh2 = (currentGene_haplotypeResolved) ? (1.0 / currentGene_MSA_ids_h2.size()) : -1; // probability to jump to a particular haplotype within the A2 group upon recombination within A2
	double change_oneh = 1.0/MSA_reference_sequences.size(); // probability to jump to any haplotype upon recombination (in the absence of A1/A2 information)

	// no A1 / A2 allele group information
	double noHaplotypeResolution_remain_oneh = within_h1h2_NoChangeTemplate + within_h1h2_changeTemplate * change_oneh; // probability to stay on exactly the same haplotype
	double noHaplotypeResolution_change_oneh = within_h1h2_changeTemplate * change_oneh; // probability to change haplotype

	// with A1 / A2 allele group information
	double haplotypeResolution_remain_h1 = noChange_h1h2 * within_h1h2_NoChangeTemplate + noChange_h1h2 * within_h1h2_changeTemplate * change_oneh1; // probability to stay on exactly the same haplotype
	double haplotypeResolution_change_h1_withinGroup = noChange_h1h2 * within_h1h2_changeTemplate * change_oneh1; // probability to switch haplotype ot another one in the same A1/A2 allele group
	double haplotypeResolution_change_h1_outOfGroup = change_h1h2 * change_oneh2; // probability to switch haplotype to another of the other A1/A2 allele group

	double haplotypeResolution_remain_h2 = noChange_h1h2 * within_h1h2_NoChangeTemplate + noChange_h1h2 * within_h1h2_changeTemplate * change_oneh2; // probability to stay on exactly the same haplotype
	double haplotypeResolution_change_h2_withinGroup = noChange_h1h2 * within_h1h2_changeTemplate * change_oneh2; // probability to switch haplotype ot another one in the same A1/A2 allele group
	double haplotypeResolution_change_h2_outOfGroup = change_h1h2 * change_oneh1; // probability to switch haplotype to another of the other A1/A2 allele group


	std::set<size_t> diseappearing_read_IDs_indices;
	if(read_stop_per_position.at(currentGene).count(first_level))
	{
		for(auto readID : read_stop_per_position.at(currentGene).at(first_level))
		{
			size_t readID_index = readID_2_index.at(readID);
			diseappearing_read_IDs_indices.insert(readID_index);
		}
	}

	std::set<size_t> new_read_IDs_indices_nextLevel;
	if(read_start_per_position.at(currentGene).count(first_level+1))
	{
		for(auto readID : read_start_per_position.at(currentGene).at(first_level+1))
		{
			size_t readID_index = readID_2_index.at(readID);
			new_read_IDs_indices_nextLevel.insert(readID_index);
		}
	}

	forReturn.reserve(statesByLevel.at(first_level).size() * 15);
	size_t n_jumps = 0;
	for(size_t stateI = 0; stateI < statesByLevel.at(first_level).size(); stateI++)
	{
		const HMMstate& s = statesByLevel.at(first_level).at(stateI);

		std::set<std::string> forward_read_assignment_states = nextLevel_compatibleReadAssignments(readAssignmentStates.at(s.readAssignmentState), diseappearing_read_IDs_indices, new_read_IDs_indices_nextLevel);

		double _sum_jump_Ps = 0;

		double read_assignment_state_p = 1.0 / forward_read_assignment_states.size();
		unsigned int n_activeAlleles_nextLevel = currentGene_activeAlleles_perPosition.at(first_level+1).size();
		assert(n_activeAlleles_nextLevel > 0);

		for(std::string readAssignmentState_nextLevel : forward_read_assignment_states)
		{
			size_t readAssignmentState_nextLevel_index = readAssignmentState_2_index.at(readAssignmentState_nextLevel);
			if(level_readAssignmentState_2_states.at(first_level+1).count(readAssignmentState_nextLevel_index) == 0)
			{
				std::cerr << "Missing a read assingment state at the next level.\n";
				std::cerr << "\t" << "first_level" << ": " << first_level << "\n";
				std::cerr << "\t" << "first_level+1" << ": " << first_level+1 << "\n";
				std::cerr << "\t" << "readAssignmentState_nextLevel_index" << ": " << readAssignmentState_nextLevel_index << "\n";
				std::cerr << "\t" << "readAssignmentState_nextLevel_index" << ": " << readAssignmentState_nextLevel_index << "\n";
				std::cerr << std::flush;
			}
			std::set<size_t> nextLevel_jumpIntoStates = level_readAssignmentState_2_states.at(first_level+1).at(readAssignmentState_nextLevel_index);

			for(auto nextStateIndex : nextLevel_jumpIntoStates)
			{
				const HMMstate& next_s = statesByLevel.at(first_level+1).at(nextStateIndex);

				double haplotype_copy_p = 1;

				if(currentGene_haplotypeResolved)
				{
					if(s.copyingFrom.first != next_s.copyingFrom.first)
					{
						if(currentGene_MSA_ids_same_haploGroup.count(std::make_pair(s.copyingFrom.first, next_s.copyingFrom.first)))
						{
							haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ? haplotypeResolution_change_h1_withinGroup : haplotypeResolution_change_h2_withinGroup);
						}
						else
						{
							haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ?  haplotypeResolution_change_h1_outOfGroup : haplotypeResolution_change_h2_outOfGroup);
						}
					}
					else
					{
						haplotype_copy_p *=  (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ? haplotypeResolution_remain_h1 : haplotypeResolution_remain_h2);
					}


					if(s.copyingFrom.second != next_s.copyingFrom.second)
					{
						if(currentGene_MSA_ids_same_haploGroup.count(std::make_pair(s.copyingFrom.second, next_s.copyingFrom.second)))
						{
							haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.second) ? haplotypeResolution_change_h1_withinGroup : haplotypeResolution_change_h2_withinGroup);
						}
						else
						{
							haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.second) ?  haplotypeResolution_change_h1_outOfGroup : haplotypeResolution_change_h2_outOfGroup);
						}
					}
					else
					{
						haplotype_copy_p *=  (currentGene_MSA_ids_h1.count(s.copyingFrom.second) ? haplotypeResolution_remain_h1 : haplotypeResolution_remain_h2);
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

				bool alleles_match_h1 = (next_s.haplotypes_alleles.first.length() == 1) && (next_s.haplotypes_alleles.first.at(0) == MSA_reference_sequences.at(currentGene).at(currentGene_MSA_int_2_id.at(next_s.copyingFrom.first)).at(first_level+1));
				bool alleles_match_h2 = (next_s.haplotypes_alleles.second.length() == 1) && (next_s.haplotypes_alleles.second.at(0) == MSA_reference_sequences.at(currentGene).at(currentGene_MSA_int_2_id.at(next_s.copyingFrom.second)).at(first_level+1));

				double allele_p =
						(alleles_match_h1 ? ((1 - currentGene_novel_allele_p) + (currentGene_novel_allele_p * (1.0 / n_activeAlleles_nextLevel))) : (currentGene_novel_allele_p * (1.0 / n_activeAlleles_nextLevel))) *
						(alleles_match_h2 ? ((1 - currentGene_novel_allele_p) + (currentGene_novel_allele_p * (1.0 / n_activeAlleles_nextLevel))) : (currentGene_novel_allele_p * (1.0 / n_activeAlleles_nextLevel)));

				double p_jump = read_assignment_state_p * haplotype_copy_p * allele_p;

				HMMtransition thisTransition;
				thisTransition.from_level = first_level;
				thisTransition.from_state = stateI;
				thisTransition.to_state = nextStateIndex;
				thisTransition.P = p_jump;
				forReturn.push_back(thisTransition);

				_sum_jump_Ps += p_jump;
				n_jumps++;
			}
		}

		// std::cerr << "Level " << first_level << " state " << stateI << " sum of transition (J)  Ps: " << _sum_jump_Ps << "\n" << std::flush;
		assert(abs(1 - _sum_jump_Ps) <= 1e-4);
	}

	return forReturn;

}
