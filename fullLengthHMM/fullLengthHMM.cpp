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
#include <vector>

#include <algorithm>
#include <omp.h>

#include <cmath>
#include <sstream>

#include "../Utilities.h"

fullLengthHMM::fullLengthHMM(
	std::map<std::string, unsigned int> _gene_length,
	std::map<std::string, std::string> _reads_2_genes,
	//std::map<std::string, std::map<std::string, std::pair<unsigned int, unsigned int>>> _read_start_stop_positions,
	std::map<std::string, std::map<unsigned int, std::set<std::string>>> _thisGene_reads_start_per_position,
	std::map<std::string, std::map<unsigned int, std::set<std::string>>> _thisGene_reads_stop_per_position,
	std::map<std::string, std::map<std::string, std::map<unsigned int, std::string>>> _read_genotypes_per_position,
	std::map<std::string, std::map<unsigned int, std::set<std::string>>> _activeAlleles_per_position,
	std::map<std::string, std::map<std::string, std::string>> _MSA_reference_sequences,
	std::map<std::string, std::map<std::string, std::string>> _MSA_reference_sequences_whichHap
) :
	gene_length(_gene_length),
	all_reads_2_genes(_reads_2_genes),
	//all_reads_start_stop_positions(_read_start_stop_positions),
	all_reads_start_per_position(_thisGene_reads_start_per_position),
	all_reads_stop_per_position(_thisGene_reads_stop_per_position),
	read_genotypes_per_position(_read_genotypes_per_position),
	activeAlleles_per_position(_activeAlleles_per_position),
	MSA_reference_sequences(_MSA_reference_sequences),
	MSA_reference_sequences_whichHap(_MSA_reference_sequences_whichHap)
{

}

std::set<std::string> fullLengthHMM::nextLevel_compatibleReadAssignments(const std::string& thisReadAssignment, const std::set<size_t>& disappearingReadIDs_nextLevel, const std::set<size_t>& newReadIDs_nextLevel) const
{
	assert(! constrainedReadAssignmentStates);
	std::vector<std::string> readAssignmentStates = {thisReadAssignment};
	for(auto disappearingReadIDIdx : disappearingReadIDs_nextLevel)
	{
		readAssignmentStates.at(0).at(disappearingReadIDIdx) = 'N';
	}

	size_t naiveStateAllocation_size = std::pow(2, newReadIDs_nextLevel.size());
	if(naiveStateAllocation_size <= (5 * readAssignmentState_2_index.size()))
	{		
		readAssignmentStates.resize(naiveStateAllocation_size, readAssignmentTemplate);

		size_t vectorElements_filled = 1;
		for(const auto& newReadIdIdx : newReadIDs_nextLevel)
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

		std::vector<std::string> filtered_readAssignmentStates;
		filtered_readAssignmentStates.reserve(readAssignmentStates.size());
		for(auto oneReadAssignmentState : readAssignmentStates)
		{
			if(readAssignmentState_2_index.count(oneReadAssignmentState))
			{
				filtered_readAssignmentStates.push_back(oneReadAssignmentState);
			}
		}
		return std::set<std::string>(filtered_readAssignmentStates.begin(), filtered_readAssignmentStates.end());
		
	}
	else
	{
		std::set<std::string> forReturn;
		for(const auto& existingAssignmentStatesIt : readAssignmentState_2_index)
		{
			const std::string& thisExistingAssignmentState = existingAssignmentStatesIt.first;
			bool existingAssignmentCompatible = true;
			for(unsigned int posInAssignmentState = 0; posInAssignmentState < thisExistingAssignmentState.size(); posInAssignmentState++)
			{
				if(newReadIDs_nextLevel.count(posInAssignmentState))
				{
					if(thisExistingAssignmentState.at(posInAssignmentState) == 'N')
					{
						existingAssignmentCompatible = false;
						break;
					}
				}
				else
				{
					if(thisExistingAssignmentState.at(posInAssignmentState) != readAssignmentStates.at(0).at(posInAssignmentState))
					{
						existingAssignmentCompatible = false;
						break;
					}
				} 
			}
			if(existingAssignmentCompatible)
			{
				forReturn.insert(thisExistingAssignmentState);
			}
		}
		return forReturn;
	}
}

std::set<std::string> fullLengthHMM::previousLevel_compatibleReadAssignments(const std::string& thisReadAssignment, const std::set<size_t>& newReadIDs_previousLevel, const std::set<size_t>& disappearingReadIDs_thisLevel) const
{
	assert(! constrainedReadAssignmentStates);
	std::vector<std::string> readAssignmentStates = {thisReadAssignment};
	for(auto disappearingReadIDIdx : disappearingReadIDs_thisLevel)
	{
		readAssignmentStates.at(0).at(disappearingReadIDIdx) = 'N';
	}

	size_t naiveStateAllocation_size = std::pow(2, newReadIDs_previousLevel.size());
	if(naiveStateAllocation_size <= (5 * readAssignmentState_2_index.size()))
	{	
		readAssignmentStates.resize(std::pow(2, newReadIDs_previousLevel.size()), readAssignmentTemplate);

		size_t vectorElements_filled = 1;
		for(const auto& newReadIdIdx : newReadIDs_previousLevel)
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

		std::vector<std::string> filtered_readAssignmentStates;
		filtered_readAssignmentStates.reserve(readAssignmentStates.size());
		for(auto oneReadAssignmentState : readAssignmentStates)
		{
			if(readAssignmentState_2_index.count(oneReadAssignmentState))
			{
				filtered_readAssignmentStates.push_back(oneReadAssignmentState);
			}
		}

		return std::set<std::string>(filtered_readAssignmentStates.begin(), filtered_readAssignmentStates.end());
	}
	else
	{
		std::set<std::string> forReturn;
		for(const auto& existingAssignmentStatesIt : readAssignmentState_2_index)
		{
			const std::string& thisExistingAssignmentState = existingAssignmentStatesIt.first;
			bool existingAssignmentCompatible = true;
			for(unsigned int posInAssignmentState = 0; posInAssignmentState < thisExistingAssignmentState.size(); posInAssignmentState++)
			{
				if(newReadIDs_previousLevel.count(posInAssignmentState))
				{
					if(thisExistingAssignmentState.at(posInAssignmentState) == 'N')
					{
						existingAssignmentCompatible = false;
						break;
					}
				}
				else
				{
					if(thisExistingAssignmentState.at(posInAssignmentState) != readAssignmentStates.at(0).at(posInAssignmentState))
					{
						existingAssignmentCompatible = false;
						break;
					}
				} 
			}
			if(existingAssignmentCompatible)
			{
				forReturn.insert(thisExistingAssignmentState);
			}
		}			
		return forReturn;
	}
}

void fullLengthHMM::trimReadsToPolymorphicPositions(const std::string& geneID, std::set<std::string>& forRet_removedReads)
{
	forRet_removedReads.clear();
	std::set<std::string> readIDs;
	std::map<std::string, unsigned int> read_start_positions;
	std::map<std::string, unsigned int> read_stop_positions;

	for(auto read2geneEntry : all_reads_2_genes)
	{
		if(read2geneEntry.second == geneID)
		{
			readIDs.insert(read2geneEntry.first);
		}
	}

	for(const auto& readStartsPerPos : all_reads_start_per_position.at(geneID))
	{
		unsigned int levelI = readStartsPerPos.first;
		for(const std::string& readID : readStartsPerPos.second)
		{
			assert(readIDs.count(readID));
			assert(read_start_positions.count(readID) == 0);
			read_start_positions[readID] = levelI;
		}
	}

	for(const auto& readStopsPerPos : all_reads_stop_per_position.at(geneID))
	{
		unsigned int levelI = readStopsPerPos.first;
		for(const std::string& readID : readStopsPerPos.second)
		{
			assert(readIDs.count(readID));
			assert(read_stop_positions.count(readID) == 0);
			read_stop_positions[readID] = levelI;
		}
	}

	std::map<std::string, std::set<unsigned int>> readID_2_hetPos;
	std::map<unsigned int, std::set<std::string>> hetPos_2_readID;

	bool verbose = false;
	if(1 || verbose)
		std::cerr << "fullLengthHMM::trimReadsToPolymorphicPositions(..): Processing " << geneID << "\n" << std::flush;

	unsigned int reads_removed = 0;
	unsigned int reads_trimmed = 0;
	unsigned int reads_unchanged = 0;
	
	for(const std::string& readID : readIDs)
	{
		unsigned int positions_heterozygous = 0;
		int minPos_heterozygous;
		int maxPos_heterozygous;
		for(unsigned int readLevelI = read_start_positions.at(readID); readLevelI <= read_stop_positions.at(readID); readLevelI++)
		{
			if(activeAlleles_per_position.at(geneID).at(readLevelI).size() > 1)
			{
				if(read_genotypes_per_position.at(geneID).at(readID).count(readLevelI))
				{
					positions_heterozygous++;
					readID_2_hetPos[readID].insert(readLevelI);
					hetPos_2_readID[readLevelI].insert(readID);
					if(positions_heterozygous == 1)
						minPos_heterozygous = readLevelI;
					maxPos_heterozygous = readLevelI;
				}
			}
		}
		if(positions_heterozygous == 0)
		{
			forRet_removedReads.insert(readID);
			all_reads_2_genes.erase(readID);
			all_reads_start_per_position.at(geneID).at(read_start_positions.at(readID)).erase(readID);
			all_reads_stop_per_position.at(geneID).at(read_stop_positions.at(readID)).erase(readID);
			reads_removed++;
			if(verbose)
				std::cerr << "\tRead ID " << readID << ", positions_heterozygous = " << positions_heterozygous << ", ignore for further analyses." << std::flush;
		}
		else
		{
			unsigned int oldStartPos = read_start_positions.at(readID);
			unsigned int oldStopPos = read_stop_positions.at(readID);
			
			assert(minPos_heterozygous <= maxPos_heterozygous);
			all_reads_start_per_position.at(geneID).at(read_start_positions.at(readID)).erase(readID);
			all_reads_stop_per_position.at(geneID).at(read_stop_positions.at(readID)).erase(readID);
			all_reads_start_per_position.at(geneID)[minPos_heterozygous].insert(readID);
			all_reads_stop_per_position.at(geneID)[maxPos_heterozygous].insert(readID);
			if(verbose)
			{
				if((oldStartPos == minPos_heterozygous) && (oldStopPos == maxPos_heterozygous))
				{
					reads_unchanged++;
					std::cerr << "\tRead ID " << readID << ", positions_heterozygous = " << positions_heterozygous << ", coordinates unchanged.\n" << std::flush;
				}
				else
				{
					reads_trimmed++;
					std::cerr << "\tRead ID " << readID << ", positions_heterozygous = " << positions_heterozygous << ", trim to " << minPos_heterozygous << " - " << maxPos_heterozygous  << " (from " << oldStartPos << " - " << oldStopPos << ")" << "\n" << std::flush;					
				}
			}
		}
	}
	
	if(1 || verbose)
	{
		std::cerr << "\n\tSummary read stats:\n";
		std::cerr << "\t\t" << "reads_unchanged" << ": " << reads_unchanged << "\n";
		std::cerr << "\t\t" << "reads_trimmed" << ": " << reads_trimmed << "\n";
		std::cerr << "\t\t" << "reads_removed" << ": " << reads_removed << "\n";
		std::cerr << std::flush;
	}	

	if(verbose)
	{
		std::cerr << "\n\tHet pos linkage:\n";
	}
	std::set<unsigned int> lonelyHetPos;
	for(auto hetPos : hetPos_2_readID)
	{
		unsigned int hetPos_levelI = hetPos.first;
		std::map<unsigned int, std::set<std::string>> linkedHetPos;
		for(const std::string& readIDOverHetPos : hetPos.second)
		{
			for(unsigned int otherHetPos : readID_2_hetPos.at(readIDOverHetPos))
			{
				if(hetPos_levelI != otherHetPos)
				{
					linkedHetPos[otherHetPos].insert(readIDOverHetPos);
				}
			}
		}
		if(0 && verbose)
		{
			std::cerr << "\t\tPosition " << hetPos.first << " linked with: \n";
			for(auto linkedHetPosIt : linkedHetPos)
			{
				std::cerr << "\t\t\t" << linkedHetPosIt.first;
				for(auto readID : linkedHetPosIt.second)
				{
					std::cerr << " " << readID;
				}
				std::cerr << "\n" << std::flush;
			}
		}
		
		if(linkedHetPos.size() == 0)
		{
			lonelyHetPos.insert(hetPos_levelI);
		}
	}

	if(1 || verbose)
	{
		std::cerr << "\n\tlonelyHetPos.size(): " << lonelyHetPos.size() << "\n" << std::flush;
	}

}

std::vector<std::string> fullLengthHMM::computeReadAssignmentSets(const std::set<std::string>& runningReadIDs, const std::map<std::string, double>& oneReadP_h1, const std::map<std::string, std::map<std::string, double>>& readPair_differentHaplotypes_P) const
{
	//assert(readAssignmentTemplate.size());
	std::vector<std::string> readAssignmentStates;

	if(runningReadIDs.size() == 0)
	{
		readAssignmentStates.push_back(readAssignmentTemplate);
	}
	else
	{
		if((oneReadP_h1.size() == 0) && (readPair_differentHaplotypes_P.size() == 0))
		{
			readAssignmentStates.resize(std::pow(2, runningReadIDs.size()), readAssignmentTemplate);

			size_t vectorElements_filled = 1;

			for(auto readID : runningReadIDs)
			{
				if(readID_2_index.count(readID) == 0)
				{
					std::cerr << "Missing readID_2_index entry for " << readID << "\n" << std::flush;
				}
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
		else
		{
			size_t reserveNow = std::pow(2, runningReadIDs.size()/2);
			if(reserveNow > 10000)
			{
				reserveNow = 10000;
			}
			readAssignmentStates.reserve(reserveNow);
			readAssignmentStates.push_back(readAssignmentTemplate);

			for(const std::string& readID : runningReadIDs)
			{
				unsigned int readID_index = readID_2_index.at(readID);
				if(oneReadP_h1.size())
				{
					assert(oneReadP_h1.count(readID));
				}
				//assert(oneReadP_h1.count(readID) || readPair_differentHaplotypes_P.count(readID));
				size_t existingElements = readAssignmentStates.size();
				bool deleteOne = false;
				std::vector<bool> deleteExistingElements;
				deleteExistingElements.resize(existingElements, false);
				
				bool verbose = false;				
				if(verbose)
					std::cout << "computeReadAssignmentSets(..) for " << readID << "\n" << std::flush;
					
				for(unsigned int existingElementI = 0; existingElementI < existingElements; existingElementI++)
				{
					std::string readID_2_h1 = readAssignmentStates.at(existingElementI);
					std::string readID_2_h2 = readAssignmentStates.at(existingElementI);
					readID_2_h1.at(readID_index) = '1';
					readID_2_h2.at(readID_index) = '2';

					if(verbose)
						std::cout << "\tExtend existing element " << existingElementI << " / " << existingElements << ": " << readAssignmentStates.at(existingElementI) << "\n" << std::flush;
												
					bool add_readID_2_h1 = true;
					bool add_readID_2_h2 = true;

					if(oneReadP_h1.count(readID))
					{
						if(verbose)
							std::cout << "\t\t" << "oneReadP_h1.at(readID)" << ": " << oneReadP_h1.at(readID) << "\n" << std::flush;
												
						if(oneReadP_h1.at(readID) >= 0.99)
						{
							add_readID_2_h2 = false;
						}
						else if(oneReadP_h1.at(readID) <= 0.01)
						{
							add_readID_2_h1 = false;
						}
					}

					if(readPair_differentHaplotypes_P.count(readID))
					{
						if(verbose)
							std::cout << "\t\t" << "read pair analysis" << "\n" << std::flush;

						for(const std::string& readID2 : runningReadIDs)
						{
							if(verbose)
								std::cout << "\t\t\t" << "readID2: " << readID2 << "\n" << std::flush;
											
							if(readID2 == readID)
							{
								break;
							}

							if(readPair_differentHaplotypes_P.at(readID).count(readID2))
							{
								if(verbose)
									std::cout << "\t\t\t\t" << "readPair_differentHaplotypes_P.at(readID).at(readID2)" << ": " << readPair_differentHaplotypes_P.at(readID).at(readID2) << "\n" << std::flush;
																						
								unsigned int readID2_index = readID_2_index.at(readID2);
								char haplotypeAssignment_otherRead = readAssignmentStates.at(existingElementI).at(readID2_index);
								
								if(verbose)
									std::cout << "\t\t\t\t" << "haplotypeAssignment_otherRead" << ": " << haplotypeAssignment_otherRead << "\n" << std::flush;
																								
								assert((haplotypeAssignment_otherRead == '1') || (haplotypeAssignment_otherRead == '2'));
								if(readPair_differentHaplotypes_P.at(readID).at(readID2) >= 0.99)
								{
									if(haplotypeAssignment_otherRead == '1')
									{
										add_readID_2_h1 = false;
									} else
									{
										add_readID_2_h2 = false;
									}
								}
								else if(readPair_differentHaplotypes_P.at(readID).at(readID2) <= 0.01)
								{
									if(haplotypeAssignment_otherRead == '1')
									{
										add_readID_2_h2 = false;
									} else
									{
										add_readID_2_h1 = false;
									}
								}
							}
						}
					}
					
					if(readID == *(runningReadIDs.begin()))
					{
						assert(readAssignmentStates.size() == 1);
						assert(add_readID_2_h1 || add_readID_2_h2);
					}
					if(verbose)
					{
						std::cout << "\t" << "add_readID_2_h1" << ": " << add_readID_2_h1 << "\n" << std::flush;
						std::cout << "\t" << "add_readID_2_h2" << ": " << add_readID_2_h2 << "\n" << std::flush;
					}
					if(add_readID_2_h1 && add_readID_2_h2)
					{
						readAssignmentStates.at(existingElementI) = readID_2_h1;
						readAssignmentStates.push_back(readID_2_h2);
					}
					else if(add_readID_2_h1)
					{
						readAssignmentStates.at(existingElementI) = readID_2_h1;
					}
					else if(add_readID_2_h2)
					{
						readAssignmentStates.at(existingElementI) = readID_2_h2;
					}
					else
					{
						deleteExistingElements.at(existingElementI) = true;
						deleteOne = true;
					}
				}
				
				if(deleteOne)
				{
					std::vector<std::string> new_readAssignmentStates;					
					new_readAssignmentStates.reserve(readAssignmentStates.size());
					for(size_t assignmentStateI = 0; assignmentStateI < readAssignmentStates.size(); assignmentStateI++)
					{
						if(! deleteExistingElements.at(assignmentStateI))
						{
							new_readAssignmentStates.push_back(readAssignmentStates.at(assignmentStateI));
						}
					}
					readAssignmentStates = new_readAssignmentStates;
					assert(readAssignmentStates.size());
				}
			}
		}
	}

	return readAssignmentStates;
}

size_t fullLengthHMM::remainingEffectiveReadsPerGene(std::string geneID, const std::set<std::string>& useReadIDs)
{
	std::set<std::string> readIDs;

	for(auto read2geneEntry : all_reads_2_genes)
	{
		if((read2geneEntry.second == geneID) && (useReadIDs.count(read2geneEntry.first)))
		{
			readIDs.insert(read2geneEntry.first);
		}
	}

	return readIDs.size();
}

size_t fullLengthHMM::maxReadAssignmentStates(std::string geneID, const std::set<std::string>& useReadIDs, const std::map<std::string, double>& oneReadP_h1, const std::map<std::string, std::map<std::string, double>>& readPair_differentHaplotypes_P)
{
	_initInternalReadStates(geneID, useReadIDs, 0);
	std::set<std::string> runningReadIDs;
	size_t maxAssignmentStates;
	for(unsigned int first_level = 0; first_level < gene_length.at(geneID); first_level++)
	{
		assert(thisGene_reads_start_per_position.count(geneID));
		if(thisGene_reads_start_per_position.at(geneID).count(first_level))
		{
			for(auto readID : thisGene_reads_start_per_position.at(geneID).at(first_level))
			{
				assert(useReadIDs.count(readID));
				runningReadIDs.insert(readID);
			}
		}

		std::vector<std::string> possibleReadAssignmentStates_thisLevel = computeReadAssignmentSets(runningReadIDs, oneReadP_h1, readPair_differentHaplotypes_P);

		if((first_level == 0) || (possibleReadAssignmentStates_thisLevel.size() > maxAssignmentStates))
		{
			maxAssignmentStates = possibleReadAssignmentStates_thisLevel.size();
		}

		if(thisGene_reads_stop_per_position.at(geneID).count(first_level))
		{
			for(auto readID : thisGene_reads_stop_per_position.at(geneID).at(first_level))
			{
				assert(runningReadIDs.count(readID));
				runningReadIDs.erase(readID);
			}
		}
	}

	return maxAssignmentStates;
}

std::set<std::string> fullLengthHMM::_initInternalReadStates(std::string geneID, const std::set<std::string>& useReadIDs, const std::vector<std::string>* forConstraint_readAssignmentStates_readIDs)
{
	std::set<std::string> readIDs;

	thisGene_reads_start_per_position.clear();
	thisGene_reads_stop_per_position.clear();
	thisGene_reads_start_per_position[geneID] = std::map<unsigned int, std::set<std::string>>();
	thisGene_reads_stop_per_position[geneID] = std::map<unsigned int, std::set<std::string>>();
	for(auto read2geneEntry : all_reads_2_genes)
	{
		if((read2geneEntry.second == geneID) && (useReadIDs.count(read2geneEntry.first)))
		{
			readIDs.insert(read2geneEntry.first);
		}
	}
	for(const auto& positionIterator : all_reads_start_per_position.at(geneID))
	{
		for(const auto& readIDIterator : positionIterator.second)
		{
			if(readIDs.count(readIDIterator))
			{
				thisGene_reads_start_per_position[geneID][positionIterator.first].insert(readIDIterator);
			}
		}
	}
	for(const auto& positionIterator : all_reads_stop_per_position.at(geneID))
	{
		for(const auto& readIDIterator : positionIterator.second)
		{
			if(readIDs.count(readIDIterator))
			{
				thisGene_reads_stop_per_position[geneID][positionIterator.first].insert(readIDIterator);
			}
		}
	}

	readID_2_index.clear();
	readIndex_2_ID.clear();
	std::vector<std::string> readIDs_vector;

	if(forConstraint_readAssignmentStates_readIDs != 0)
	{
		if(!(forConstraint_readAssignmentStates_readIDs->size() == readIDs.size()))
		{
			std::cerr << "forConstraint_readAssignmentStates_readIDs->size()" << ": " << forConstraint_readAssignmentStates_readIDs->size() << "\n";
			std::cerr << "readIDs.size()" << ": " << readIDs.size() << "\n";
			std::cerr << "\n" << std::flush;
		}
		assert(forConstraint_readAssignmentStates_readIDs->size() == readIDs.size());
		for(auto readID : *forConstraint_readAssignmentStates_readIDs)
		{
			assert(readIDs.count(readID));
		}
		readIDs_vector = *forConstraint_readAssignmentStates_readIDs;
	}
	else
	{
		readIDs_vector = std::vector<std::string>(readIDs.begin(), readIDs.end());
	}

	readIndex_2_ID = readIDs_vector;
	for(unsigned int i = 0; i < readIDs_vector.size(); i++)
	{
		readID_2_index[readIDs_vector.at(i)] = i;
	}

	readAssignmentTemplate.clear();
	//assert(readIDs.size());
	readAssignmentTemplate.resize(readIDs.size(), 'N');
	assert(readAssignmentTemplate.size() == readIDs.size());
	
	return readIDs;
}

double fullLengthHMM::makeInference(std::string geneID, bool outputToFilestreams, std::ofstream& output_fasta, std::ofstream& output_graphLevels, std::string outputPrefix_furtherOutput, const std::set<std::string>& useReadIDs, std::map<std::string, double>& forRet_oneReadP_h1, std::map<std::pair<std::string, std::string>, double>& forRet_readPair_differentHaplotypes_P, std::map<unsigned int, std::map<std::pair<std::string, std::string>, double>>& forRet_genotypes_P, std::map<unsigned int, std::pair<std::map<std::string, double>, std::map<std::string, double>>>& forRet_allele_by_haplotype_P, std::vector<std::vector<std::string>>* forRet_samples_readAssignmentStates, std::vector<std::string>* forRet_readAssignmentStates_readIDs, size_t generateHaplotypeSamples, const std::map<std::string, double>* forConstraint_oneReadP_h1, const std::map<std::string, std::map<std::string, double>>* forConstraint_readPair_differentHaplotypes_P, const std::vector<std::vector<std::string>>* forConstraint_readAssignmentStates, const std::vector<std::string>* forConstraint_readAssignmentStates_readIDs)
{
	forRet_readPair_differentHaplotypes_P.clear();
	forRet_oneReadP_h1.clear();
	forRet_genotypes_P.clear();
	forRet_allele_by_haplotype_P.clear();

	std::cout << "makeInference" << std::flush;

	omp_set_num_threads(16);
	
	assert(gene_length.count(geneID));
	currentGene = geneID;

	currentGene_haplotypeResolved = true;
	unsigned int MSA_h1_n = 0;
	unsigned int MSA_h2_n = 0;

	currentGene_geneLength = gene_length.at(geneID);
	std::cout << "Now making HMM-based inference for gene '" << geneID << "' --length " << currentGene_geneLength << "\n" << std::flush;

	currentGene_novel_allele_p = 3.0/currentGene_geneLength;

	currentGene_MSA_ids.clear();
	currentGene_MSA_id_2_int.clear();
	currentGene_MSA_int_2_id.clear();
	currentGene_MSA_int_2_HaplotypeID.clear();

	for(auto MSAentry : MSA_reference_sequences.at(currentGene))
	{
		currentGene_MSA_ids.insert(MSAentry.first);
		currentGene_MSA_id_2_int[MSAentry.first] = currentGene_MSA_ids.size() - 1;
		currentGene_MSA_int_2_id.push_back(MSAentry.first);
		currentGene_MSA_int_2_HaplotypeID.push_back(MSA_reference_sequences_whichHap.at(currentGene).at(MSAentry.first));
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
	currentGene_activeAlleles_perPosition_allelesLength1_chars.clear();
	currentGene_activeAlleles_perPosition_allelesLength1_chars.resize(currentGene_geneLength);
	
	for(unsigned int i = 0; i < currentGene_geneLength; i++)
	{
		std::vector<std::string> activeAlleles_vec(activeAlleles_per_position.at(currentGene).at(i).begin(), activeAlleles_per_position.at(currentGene).at(i).end());
		currentGene_activeAlleles_perPosition.at(i) = activeAlleles_vec;
		for(const std::string& a : activeAlleles_vec)
		{
			if(a.length() == 1)
			{
				currentGene_activeAlleles_perPosition_allelesLength1_chars.at(i).insert(a.at(0));
			}
		}
	}

	statesByLevel.clear();
	statesByLevel.resize(currentGene_geneLength);

	level_readAssignmentState_2_states.clear();
	level_readAssignmentState_2_states.resize(currentGene_geneLength);

	std::vector<long long> states_per_position;
	states_per_position.resize(currentGene_geneLength, 0);

	std::set<std::string> readIDs = _initInternalReadStates(geneID, useReadIDs, forConstraint_readAssignmentStates_readIDs);

	// assert(readID_2_index.count("HLAA_h0_A*02:90_14_568_1:0:0_2:0:0_63"));

	readAssignmentStates.clear();
	readAssignmentState_2_index.clear();

	assert(
			((forConstraint_readAssignmentStates == 0) && (forConstraint_readAssignmentStates_readIDs == 0)) ||
			((forConstraint_readAssignmentStates != 0) && (forConstraint_readAssignmentStates_readIDs != 0))
	);

	constrainedReadAssignmentStates = (forConstraint_readAssignmentStates != 0);
	constrainedReadAssignmentStates_transitions_forward.clear();
	constrainedReadAssignmentStates_transitions_backward.clear();

	std::cout << "Gene " << currentGene << ", starting inference with " << readIDs.size() << " reads.\n" << std::flush;

	{
		tr_change_h1h2 = 0.1 * (1.0 / currentGene_geneLength); // probability to recombine into the other A1/A2 allele group initiate a haplotype recombination event to the other A1/A2 group
		tr_noChange_h1h2 = 1 - tr_change_h1h2; // probability not to recombine into the other A1/A2 allele group

		tr_within_h1h2_changeTemplate = 1.0 / currentGene_geneLength; // probability to initiate a recombination event within the same A1/A2 group, conditional on not initiating a recombination event to the other A1/A2 allele group
		tr_within_h1h2_NoChangeTemplate = 1 - tr_within_h1h2_changeTemplate; // probability not to initiate a recombination event, conditional on not initiating a recombination event to the other A1/A2 allele group

		tr_change_oneh1 = (currentGene_haplotypeResolved) ? (1.0 / currentGene_MSA_ids_h1.size()) : -1; // probability to jump to a particular haplotype within the A1 group upon recombination within A1
		tr_change_oneh2 = (currentGene_haplotypeResolved) ? (1.0 / currentGene_MSA_ids_h2.size()) : -1; // probability to jump to a particular haplotype within the A2 group upon recombination within A2
		tr_change_oneh = 1.0/MSA_reference_sequences.at(currentGene).size(); // probability to jump to any haplotype upon recombination (in the absence of A1/A2 information)

		// no A1 / A2 allele group information
		tr_noHaplotypeResolution_remain_oneh = tr_within_h1h2_NoChangeTemplate + tr_within_h1h2_changeTemplate * tr_change_oneh; // probability to stay on exactly the same haplotype
		tr_noHaplotypeResolution_change_oneh = tr_within_h1h2_changeTemplate * tr_change_oneh; // probability to change haplotype

		// with A1 / A2 allele group information
		tr_haplotypeResolution_remain_h1 = tr_noChange_h1h2 * tr_within_h1h2_NoChangeTemplate + tr_noChange_h1h2 * tr_within_h1h2_changeTemplate * tr_change_oneh1; // probability to stay on exactly the same haplotype
		tr_haplotypeResolution_change_h1_withinGroup = tr_noChange_h1h2 * tr_within_h1h2_changeTemplate * tr_change_oneh1; // probability to switch haplotype ot another one in the same A1/A2 allele group
		tr_haplotypeResolution_change_h1_outOfGroup = tr_change_h1h2 * tr_change_oneh2; // probability to switch haplotype to another of the other A1/A2 allele group

		tr_haplotypeResolution_remain_h2 = tr_noChange_h1h2 * tr_within_h1h2_NoChangeTemplate + tr_noChange_h1h2 * tr_within_h1h2_changeTemplate * tr_change_oneh2; // probability to stay on exactly the same haplotype
		tr_haplotypeResolution_change_h2_withinGroup = tr_noChange_h1h2 * tr_within_h1h2_changeTemplate * tr_change_oneh2; // probability to switch haplotype ot another one in the same A1/A2 allele group
		tr_haplotypeResolution_change_h2_outOfGroup = tr_change_h1h2 * tr_change_oneh1; // probability to switch haplotype to another of the other A1/A2 allele group
	}

	std::map<std::string, double> oneReadP_h1_empty;
	std::map<std::string, std::map<std::string, double>>  readPair_differentHaplotypes_P_empty;

	const std::map<std::string, double>* oneReadP_h1 = (forConstraint_oneReadP_h1 == 0) ? &oneReadP_h1_empty : forConstraint_oneReadP_h1;
	const std::map<std::string, std::map<std::string, double>>*  readPair_differentHaplotypes_P = (forConstraint_readPair_differentHaplotypes_P == 0) ? &readPair_differentHaplotypes_P_empty : forConstraint_readPair_differentHaplotypes_P;

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
				double tr_haplotypeResolution_change_h2_outOfGroup = change_h1h2 * change_oneh1;

				if(s.copyingFrom.first != next_s.copyingFrom.first)
				{
					if(currentGene_MSA_ids_same_haploGroup.count(std::make_pair(s.copyingFrom.first, next_s.copyingFrom.first)))
					{
						haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ? haplotypeResolution_change_h1_withinGroup : haplotypeResolution_change_h2_withinGroup);
					}
					else
					{
						haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ?  haplotypeResolution_change_h1_outOfGroup : tr_haplotypeResolution_change_h2_outOfGroup);
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
						double tr_haplotypeResolution_change_h2_outOfGroup = change_h1h2 * change_oneh1;

						if(s.copyingFrom.first != next_s.copyingFrom.first)
						{
							if(currentGene_MSA_ids_same_haploGroup.count(std::make_pair(s.copyingFrom.first, next_s.copyingFrom.first)))
							{
								haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ? haplotypeResolution_change_h1_withinGroup : haplotypeResolution_change_h2_withinGroup);
							}
							else
							{
								haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ?  haplotypeResolution_change_h1_outOfGroup : tr_haplotypeResolution_change_h2_outOfGroup);
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
								haplotype_copy_p *= (currentGene_MSA_ids_h1.count(s.copyingFrom.second) ?  haplotypeResolution_change_h1_outOfGroup : tr_haplotypeResolution_change_h2_outOfGroup);
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

	if(forConstraint_readAssignmentStates != 0)
	{
		for(const auto& oneReadAssignmentPath : *forConstraint_readAssignmentStates)
		{
			if(!(oneReadAssignmentPath.size() == currentGene_geneLength))
			{
				std::cerr << "oneReadAssignmentPath.size()" << ": " << oneReadAssignmentPath.size() << "\n";
				std::cerr << "currentGene_geneLength" << ": " << currentGene_geneLength << "\n";
				std::cerr << "\n" << std::flush;
			}			
			assert(oneReadAssignmentPath.size() == currentGene_geneLength);
		}
	}

	std::vector<unsigned int> readAssingmentStates;
	std::vector<std::vector<size_t>> constrainedReadAssignmentStates_numerical_byLevel;
	std::set<std::string> runningReadIDs;
	size_t n_states_total = 0;
	int recompute_readAssingmentStates = 0;
	constrainedReadAssignmentStates_numerical_byLevel.resize(currentGene_geneLength);
	for(unsigned int first_level = 0; first_level < currentGene_geneLength; first_level++)
	{
		if((first_level % 1000) == 0)
			std::cerr << "Round I: Level " << first_level << " / " << currentGene_geneLength << "\n" << std::flush;

		if(thisGene_reads_start_per_position.at(currentGene).count(first_level))
		{
			for(auto readID : thisGene_reads_start_per_position.at(currentGene).at(first_level))
			{
				assert(useReadIDs.count(readID));
				runningReadIDs.insert(readID);
				recompute_readAssingmentStates++;
			}
		}

		size_t n_possibleReadAssignmentStates_thisLevel;
		std::vector<std::string> possibleReadAssignmentStates_thisLevel;
		if(forConstraint_readAssignmentStates == 0)
		{
			possibleReadAssignmentStates_thisLevel = computeReadAssignmentSets(runningReadIDs, *oneReadP_h1, *readPair_differentHaplotypes_P);
			for(size_t readAssignmentStateI = 0; readAssignmentStateI < possibleReadAssignmentStates_thisLevel.size(); readAssignmentStateI++)
			{
				const std::string& thisReadAssignmentState = possibleReadAssignmentStates_thisLevel.at(readAssignmentStateI);
				if(readAssignmentState_2_index.count(thisReadAssignmentState) == 0)
				{
					readAssignmentStates.push_back(thisReadAssignmentState);
					readAssignmentState_2_index[thisReadAssignmentState] = readAssignmentStates.size() - 1;

				}
			}
			n_possibleReadAssignmentStates_thisLevel = possibleReadAssignmentStates_thisLevel.size();
		}
		else
		{
			n_possibleReadAssignmentStates_thisLevel = forConstraint_readAssignmentStates->size();
			constrainedReadAssignmentStates_numerical_byLevel.at(first_level).reserve(n_possibleReadAssignmentStates_thisLevel);
			possibleReadAssignmentStates_thisLevel.reserve(n_possibleReadAssignmentStates_thisLevel);
			for(const auto& oneReadAssignmentPath : *forConstraint_readAssignmentStates)
			{
				possibleReadAssignmentStates_thisLevel.push_back(oneReadAssignmentPath.at(first_level));
				readAssignmentStates.push_back(oneReadAssignmentPath.at(first_level));
				constrainedReadAssignmentStates_numerical_byLevel.at(first_level).push_back(readAssignmentStates.size() - 1);
			}
		}
		assert(activeAlleles_per_position.at(currentGene).count(first_level));
		assert(activeAlleles_per_position.at(currentGene).at(first_level).size() >= 1);

		size_t expectedStates_thisLevel =
					((currentGene_MSA_int_2_id.size() * currentGene_MSA_int_2_id.size() - currentGene_MSA_int_2_id.size()) / 2.0 + currentGene_MSA_int_2_id.size()) *
					n_possibleReadAssignmentStates_thisLevel *
					(currentGene_activeAlleles_perPosition.at(first_level).size() * currentGene_activeAlleles_perPosition.at(first_level).size());
		statesByLevel.at(first_level).reserve(expectedStates_thisLevel);

		for(unsigned int copyFromI_1 = 0; copyFromI_1 < currentGene_MSA_int_2_id.size(); copyFromI_1++)
		{
			for(unsigned int copyFromI_2 = 0; copyFromI_2 < currentGene_MSA_int_2_id.size(); copyFromI_2++)
			{
				for(unsigned int readAssignmentStateI = 0; readAssignmentStateI < n_possibleReadAssignmentStates_thisLevel; readAssignmentStateI++)
				{
					size_t thisReadAssignmentState_index;
					if(forConstraint_readAssignmentStates == 0)
					{
						const std::string& thisReadAssignmentState = possibleReadAssignmentStates_thisLevel.at(readAssignmentStateI);
						thisReadAssignmentState_index = readAssignmentState_2_index.at(thisReadAssignmentState);
					}
					else
					{
						thisReadAssignmentState_index = constrainedReadAssignmentStates_numerical_byLevel.at(first_level).at(readAssignmentStateI);
					}

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
							n_states_total++;

							size_t s_index = statesByLevel.at(first_level).size() - 1;
							level_readAssignmentState_2_states.at(first_level)[thisReadAssignmentState_index].insert(s_index);
						}
					}
				}
			}
		}

		if(thisGene_reads_stop_per_position.at(currentGene).count(first_level))
		{
			for(auto readID : thisGene_reads_stop_per_position.at(currentGene).at(first_level))
			{
				assert(runningReadIDs.count(readID));
				runningReadIDs.erase(readID);
				recompute_readAssingmentStates++;
			}
		}
	}

	if(forConstraint_readAssignmentStates != 0)
	{
		for(size_t readAssignmentPathI = 0; readAssignmentPathI < forConstraint_readAssignmentStates->size(); readAssignmentPathI++)
		{
			const auto& oneReadAssignmentPath = forConstraint_readAssignmentStates->at(readAssignmentPathI);
			for(size_t levelI = 0; levelI < (currentGene_geneLength-1); levelI++)
			{
				size_t thisReadAssignmentState_index = constrainedReadAssignmentStates_numerical_byLevel.at(levelI).at(readAssignmentPathI);
				size_t nextReadAssignmentState_index = constrainedReadAssignmentStates_numerical_byLevel.at(levelI+1).at(readAssignmentPathI);
				assert(readAssignmentStates.at(thisReadAssignmentState_index) == forConstraint_readAssignmentStates->at(readAssignmentPathI).at(levelI));
				assert(constrainedReadAssignmentStates_transitions_forward.count(thisReadAssignmentState_index) == 0);
				assert(constrainedReadAssignmentStates_transitions_backward.count(nextReadAssignmentState_index) == 0);
				constrainedReadAssignmentStates_transitions_forward[thisReadAssignmentState_index] = nextReadAssignmentState_index;
				constrainedReadAssignmentStates_transitions_backward[nextReadAssignmentState_index] = thisReadAssignmentState_index;
			}
		}
		
		for(unsigned int levelI = 0; levelI < currentGene_geneLength; levelI++)
		{
			for(const HMMstate& s : statesByLevel.at(levelI))
			{
				size_t readAssignmentState = s.readAssignmentState;				
				if((levelI + 1) < currentGene_geneLength) 
				{
					assert(constrainedReadAssignmentStates_transitions_forward.count(readAssignmentState));
				}				
				if(levelI > 0)
				{
					assert(constrainedReadAssignmentStates_transitions_backward.count(readAssignmentState));					
				}
			}
		}



		
		/*
		constrainedReadAssignmentStates_transitions.resize(currentGene_geneLength);
		constrainedReadAssignmentStates_transitions_backward.resize(currentGene_geneLength);
		for(size_t levelI = 0; levelI < currentGene_geneLength; levelI++)
		{
			size_t readAssignmentStatesAtLevel = level_readAssignmentState_2_states.at(levelI).size();
			constrainedReadAssignmentStates_transitions.reserve(readAssignmentStatesAtLevel);
			constrainedReadAssignmentStates_transitions_backward.reserve(readAssignmentStatesAtLevel);
			if((levelI + 1) < currentGene_geneLength)
			{
				for(size_t stateI = 0; stateI < readAssignmentStatesAtLevel; stateI++)
				{
					const std::string& thisAssignmentState = forConstraint_readAssignmentStates->at(levelI).at(levelI);
					const std::string& nextAssignmentState = forConstraint_readAssignmentStates->at(levelI).at(levelI);
				}
			}
		}
		*/
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
				// std::cout << "For assignment string " << readAssignmentString << ", set " << readID << " to " << readAssignmentString.at(readI) << "\n" << std::flush;
			}
		}
	}
	
	// assert( 1== 0);
	
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

	bool paranoid = true;
	if(paranoid)
	{
		std::cerr << "Check that computed transition probabilites in forward and backward direction agree...\n" << std::flush;
		for(unsigned int levelI = 0; levelI < currentGene_geneLength; levelI++)
		{
			// std::cerr << "\t\tCheck level " << levelI << "\n";
			if((levelI % 1000) == 0)
				std::cerr << "Round (Paranoia): Level " << levelI << " / " << currentGene_geneLength << "\n" << std::flush;

			std::map<std::pair<size_t, size_t>, double> map_forward;
			std::map<std::pair<size_t, size_t>, double> map_backward;
			if(levelI > 0)
			{
				std::vector<HMMtransition> transitions_to_next_level = computeLevelTransitions(levelI-1);
				for(auto transition : transitions_to_next_level)
				{
					map_forward[make_pair(transition.from_state, transition.to_state)] = transition.P;
				}
				std::vector<HMMtransition> transitions_to_next_level_backward = computeLevelTransitions_backward(levelI);
				for(auto transition : transitions_to_next_level_backward)
				{
					map_backward[make_pair(transition.from_state, transition.to_state)] = transition.P;
				}
				assert(transitions_to_next_level.size() == transitions_to_next_level_backward.size());
				//std::cerr << "\t\t\t" << transitions_to_next_level.size() << " transisitions agree.\n";
				for(auto statePair : map_forward)
				{
					std::pair<size_t, size_t> k = statePair.first;
					assert(map_forward.at(k) == map_backward.at(k));
				}
			}
		}
		std::cerr << "\t... success.\n" << std::flush;
		// assert(1 == 0);
	}


	unsigned int normalization_interval = 100;
	double normalization_forward_sum_log = 0;
	double normalization_viterbi_sum_log = 0;
	double normalization_backward_sum_log = 0;

	size_t n_jumps = 0;
	for(unsigned int levelI = 0; levelI < currentGene_geneLength; levelI++)
	{
		if((levelI % 1000) == 0)
			std::cerr << "Round II: Level " << levelI << " / " << currentGene_geneLength << "\n" << std::flush;

		std::vector<double> emissionP = computeEmissionProbabilities(levelI);

		std::vector<std::map<size_t, double>> states_levelI_jumpFrom;
		size_t n_states = statesByLevel.at(levelI).size();

		// std::cerr << "\t" << "A" << "\n" << std::flush;

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

		// std::cerr << "\t" << "B" << "\n" << std::flush;

		#pragma omp parallel for
		for(size_t stateI = 0; stateI < n_states ; stateI++)
		{
			HMMstate& s = statesByLevel.at(levelI).at(stateI);
			if(levelI == 0)
			{
				s.fw_p = initialProbabilities.at(stateI) * emissionP.at(stateI);
				s.viterbi_p = s.fw_p;
			}
			else
			{
				double running_viterbi_max_p = 0;
				size_t running_viterbi_max_whereFrom;

				s.fw_p = 0;
				//assert(states_levelI_jumpFrom.at(stateI).size());
				
				for(auto jumpIntoThisState : states_levelI_jumpFrom.at(stateI))
				{
					s.fw_p += statesByLevel.at(levelI-1).at(jumpIntoThisState.first).fw_p * jumpIntoThisState.second;
					assert(s.fw_p >= 0);

					double thisJump_viterbi = statesByLevel.at(levelI-1).at(jumpIntoThisState.first).viterbi_p * jumpIntoThisState.second;
					if(!(thisJump_viterbi >= 0))
					{
						#pragma omp critical
						{
							std::cerr << "jumpIntoThisState.first" << ": " << jumpIntoThisState.first << "\n";
							std::cerr << "thisJump_viterbi" << ": " << thisJump_viterbi << "\n";
							std::cerr << "levelI" << ": " << levelI << "\n";
							std::cerr << "statesByLevel.at(levelI-1).at(jumpIntoThisState.first).viterbi_p" << ": " << statesByLevel.at(levelI-1).at(jumpIntoThisState.first).viterbi_p << "\n";
							std::cerr << "jumpIntoThisState.second" << ": " << jumpIntoThisState.second << "\n";
							std::cerr << std::flush;
						}
					}
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

		// std::cerr << "\t" << "C" << "\n" << std::flush;

		if((levelI % normalization_interval) == 0)
		{
			double sum_forward_level = 0;
			double max_forward_viterbi = 0;
			for(size_t stateI = 0; stateI < n_states ; stateI++)
			{
				sum_forward_level += statesByLevel.at(levelI).at(stateI).fw_p;
				if((stateI == 0) || (statesByLevel.at(levelI).at(stateI).viterbi_p > max_forward_viterbi))
				{
					max_forward_viterbi = statesByLevel.at(levelI).at(stateI).viterbi_p;
				}
			}
			assert(sum_forward_level > 0);
			assert(max_forward_viterbi > 0);

			#pragma omp parallel for
			for(size_t stateI = 0; stateI < n_states ; stateI++)
			{
				statesByLevel.at(levelI).at(stateI).fw_p = statesByLevel.at(levelI).at(stateI).fw_p / sum_forward_level;
				statesByLevel.at(levelI).at(stateI).viterbi_p = statesByLevel.at(levelI).at(stateI).viterbi_p / max_forward_viterbi;
				assert((statesByLevel.at(levelI).at(stateI).fw_p >= 0) && (statesByLevel.at(levelI).at(stateI).fw_p <= 1));
				assert((statesByLevel.at(levelI).at(stateI).viterbi_p >= 0) && (statesByLevel.at(levelI).at(stateI).viterbi_p <= 1));
			}
			assert(sum_forward_level > 0);
			assert(max_forward_viterbi > 0);

			normalization_forward_sum_log += log(sum_forward_level);
			normalization_viterbi_sum_log += log(max_forward_viterbi);
		}

		// std::cerr << "\t" << "D" << "\n" << std::flush;
	}

	double forward_joint_p_normalized = 0;
	for(const auto& s: statesByLevel.at(statesByLevel.size() - 1))
	{
		forward_joint_p_normalized += s.fw_p;
	}
	assert(forward_joint_p_normalized > 0);

	double forward_joint_log_p = normalization_forward_sum_log + log(forward_joint_p_normalized);
	std::cout << "Joint forward log P: " << forward_joint_log_p << "\n" << std::flush;

	bool doBackward = false;
	if(doBackward)
	{
		for(long long levelI = (currentGene_geneLength - 1); levelI >= 0; levelI--)
		{
			if((levelI % 1000) == 0)
				std::cerr << "Round III (Backward): Level " << levelI << "\n" << std::flush;

			std::vector<double> emissionP = computeEmissionProbabilities(levelI);

			size_t n_states = statesByLevel.at(levelI).size();
			if(levelI == (currentGene_geneLength - 1))
			{
				for(size_t stateI = 0; stateI < n_states; stateI++)
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

				#pragma omp parallel for
				for(size_t stateI = 0; stateI < n_states; stateI++)
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

			if((levelI > 0) && (levelI % normalization_interval) == 0)
			{
				double sum_backward_level = 0;
				for(size_t stateI = 0; stateI < n_states ; stateI++)
				{
					sum_backward_level += statesByLevel.at(levelI).at(stateI).bw_p;
				}
				assert(sum_backward_level > 0);

				#pragma omp parallel for
				for(size_t stateI = 0; stateI < n_states ; stateI++)
				{
					statesByLevel.at(levelI).at(stateI).bw_p = statesByLevel.at(levelI).at(stateI).bw_p / sum_backward_level;
					assert((statesByLevel.at(levelI).at(stateI).bw_p >= 0) && (statesByLevel.at(levelI).at(stateI).bw_p <= 1));
				}

				normalization_backward_sum_log += log(sum_backward_level);
			}
		}

		double backward_joint_p_normalized = 0;
		for(const auto& s: statesByLevel.at(0))
		{
			backward_joint_p_normalized += s.bw_p;
		}
		double backward_joint_log_p = normalization_backward_sum_log + log(backward_joint_p_normalized);

		std::cout << "Joint forward log P: " << forward_joint_log_p << "\n" << std::flush;
		std::cout << "Joint backward log P: " << backward_joint_log_p << "\n" << std::flush;

		double log_diff = forward_joint_log_p - backward_joint_log_p;
		double pror_diff = exp(log_diff);

		std::cout << "log_diff: " << log_diff << "\n" << std::flush;
		std::cout << "pror_diff: " << pror_diff << "\n" << std::flush;

		assert(abs(1 - pror_diff) <= 1e-2);
	}

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
	
	std::vector<std::string> bt_h1_separated;
	std::vector<std::string> bt_h2_separated;
	
	std::vector<std::string> bt_graphLevels;
	

	{
		bt_h1.reserve(currentGene_geneLength);
		bt_h2.reserve(currentGene_geneLength);

		bt_h1_separated.reserve(currentGene_geneLength);
		bt_h2_separated.reserve(currentGene_geneLength);

		bt_graphLevels.reserve(currentGene_geneLength);
		
		for(unsigned int levelI = 0; levelI < (currentGene_geneLength - 1); levelI++)
		{
			const HMMstate& s = statesByLevel.at(levelI).at(backtrace_Viterbi.at(levelI));
			bt_h1.append(s.haplotypes_alleles.first);
			bt_h2.append(s.haplotypes_alleles.second);

			bt_h1_separated.push_back(s.haplotypes_alleles.first);
			bt_h2_separated.push_back(s.haplotypes_alleles.second);
			
			bt_graphLevels.push_back(std::to_string(levelI));
		}

		//assert(bt_h1.size() == currentGene_geneLength);
		//assert(bt_h2.size() == currentGene_geneLength);
	}

	if(outputToFilestreams)
	{
		output_fasta << ">" << geneID << "-H1\n";
		output_fasta << bt_h1 << "\n";
		output_fasta << ">" << geneID << "-H2\n";
		output_fasta << bt_h2 << "\n";
		output_fasta << std::flush;

		output_graphLevels << ">" << geneID << "-H1\n";
		output_graphLevels << Utilities::join(bt_h1_separated, ";") << "\n";
		output_graphLevels << ">" << geneID << "-H2\n";
		output_graphLevels << Utilities::join(bt_h2_separated, ";") << "\n";
		output_graphLevels << ">" << geneID << "-Levels\n";
		output_graphLevels << Utilities::join(bt_graphLevels, ";") << "\n";
		output_graphLevels << std::flush;
	}
	
	std::cout << currentGene << " done -- " << n_states_total << " states -- " << n_jumps << " jumps.\n" << std::flush;

	std::string outputFn_sampledCompleteHaplotypes = outputPrefix_furtherOutput + ".sampledCompleteHaplotypes";
	std::ofstream outputStream_sampledCompleteHaplotypes;
	outputStream_sampledCompleteHaplotypes.open(outputFn_sampledCompleteHaplotypes.c_str());
	assert(outputStream_sampledCompleteHaplotypes.is_open());

	std::string outputFn_sampledGenotypesAndCopiedFrom = outputPrefix_furtherOutput + ".genotypesAndCopiedFromPerPosition";
	std::ofstream outputStream_sampledGenotypesAndCopiedFrom;
	outputStream_sampledGenotypesAndCopiedFrom.open(outputFn_sampledGenotypesAndCopiedFrom.c_str());
	assert(outputStream_sampledGenotypesAndCopiedFrom.is_open());

	std::string outputFn_readID2Haplotype = outputPrefix_furtherOutput + ".readID2HaplotypeAssignment";
	std::ofstream outputStream_readID2Haplotype;
	outputStream_readID2Haplotype.open(outputFn_readID2Haplotype.c_str());
	assert(outputStream_readID2Haplotype.is_open());

	std::string outputFn_readIDSameHaplotype = outputPrefix_furtherOutput + ".readIDSameHaplotype";
	std::ofstream outputStream_readIDSameHaplotype;
	outputStream_readIDSameHaplotype.open(outputFn_readIDSameHaplotype.c_str());
	assert(outputStream_readIDSameHaplotype.is_open());

	std::vector<std::map<std::string, unsigned int>> samples_genotypes_by_position;
	std::vector<std::map<std::string, unsigned int>> samples_copyingFrom_by_position;
	std::vector<std::map<std::string, unsigned int>> samples_copyingFromPlusAllele_by_position;
	std::map<std::string, std::map<std::string, unsigned int>> readID_2_haplotype;
	samples_genotypes_by_position.resize(currentGene_geneLength);
	samples_copyingFrom_by_position.resize(currentGene_geneLength);
	samples_copyingFromPlusAllele_by_position.resize(currentGene_geneLength);

	std::map<std::pair<std::string, std::string>, std::map<unsigned int, std::pair<unsigned int, unsigned int>>> readPairs_haplotypeIdentical;

	std::set<std::pair<std::string, std::string>> overlappingReadIDPairs;

	if(forRet_samples_readAssignmentStates != 0)
	{
		forRet_samples_readAssignmentStates->clear();
		forRet_samples_readAssignmentStates->resize(generateHaplotypeSamples);
	}

	if(forRet_readAssignmentStates_readIDs != 0)
	{
		forRet_readAssignmentStates_readIDs->clear();
		*forRet_readAssignmentStates_readIDs = readIndex_2_ID;
	}

	for(size_t sampleI = 0; sampleI < generateHaplotypeSamples; sampleI++)
	{
		// std::cout << "Sample " << sampleI << " / " << generateHaplotypeSamples << "\n" << std::flush;
		
		std::vector<size_t> thisSample;
		thisSample.reserve(currentGene_geneLength);
		size_t runningState_lastLevel;
		for(long long levelI = (currentGene_geneLength - 1); levelI >= 0; levelI--)
		{
			std::vector<double> currentLevel_selectFrom;
			currentLevel_selectFrom.reserve(statesByLevel.at(levelI).size());
			if(levelI == (currentGene_geneLength - 1))
			{
				for(const auto& s: statesByLevel.at(levelI))
				{
					currentLevel_selectFrom.push_back(s.fw_p);
				}
				size_t thisLevel_chosenState = Utilities::chooseFromVector(currentLevel_selectFrom);
				thisSample.push_back(thisLevel_chosenState);
				runningState_lastLevel = thisLevel_chosenState;
			}
			else
			{
				std::vector<HMMtransition> transitions_from = computeLevelTransitions_backward(levelI+1, runningState_lastLevel);
				for(const auto& s : transitions_from)
				{
					assert(s.from_level == levelI);
					assert(s.to_state == runningState_lastLevel);
					double weight = s.P * statesByLevel.at(levelI).at(s.from_state).fw_p;
					currentLevel_selectFrom.push_back(weight);
				}
				size_t thisLevel_chosenState = Utilities::chooseFromVector(currentLevel_selectFrom);
				thisSample.push_back(transitions_from.at(thisLevel_chosenState).from_state);
				runningState_lastLevel = transitions_from.at(thisLevel_chosenState).from_state;
			}
		}
		assert(thisSample.size() == currentGene_geneLength);
		std::reverse(thisSample.begin(), thisSample.end());

		std::string h1;
		std::string h2;
		{
			h1.reserve(currentGene_geneLength);
			h1.reserve(currentGene_geneLength);
			if(forRet_samples_readAssignmentStates != 0)
			{
				forRet_samples_readAssignmentStates->at(sampleI).resize(currentGene_geneLength);
			}

			std::map<std::string, std::string> thisSample_readID_2_haplotype;
			for(unsigned int levelI = 0; levelI < currentGene_geneLength; levelI++)
			{
				const HMMstate& s = statesByLevel.at(levelI).at(thisSample.at(levelI));
				h1.append(s.haplotypes_alleles.first);
				h2.append(s.haplotypes_alleles.second);

				std::string genotype = s.haplotypes_alleles.first + "|" + s.haplotypes_alleles.second;
				if(samples_genotypes_by_position.at(levelI).count(genotype) == 0)
					samples_genotypes_by_position.at(levelI)[genotype] = 0;

				samples_genotypes_by_position.at(levelI).at(genotype)++;

				std::string copyingFrom = currentGene_MSA_int_2_id.at(s.copyingFrom.first) + "|" + currentGene_MSA_int_2_id.at(s.copyingFrom.second);
				if(samples_copyingFrom_by_position.at(levelI).count(copyingFrom) == 0)
					samples_copyingFrom_by_position.at(levelI)[copyingFrom] = 0;

				samples_copyingFrom_by_position.at(levelI).at(copyingFrom)++;

				double oneSampleWeight = 1.0/double(generateHaplotypeSamples);
				std::pair<std::string, std::string> gt_pair = makeGt(s.haplotypes_alleles.first, s.haplotypes_alleles.second);
				if(forRet_genotypes_P[levelI].count(gt_pair) == 0)
				{
					forRet_genotypes_P[levelI][gt_pair] = 0;
				}
				forRet_genotypes_P.at(levelI).at(gt_pair) += oneSampleWeight;

				if(forRet_allele_by_haplotype_P[levelI].first.count(s.haplotypes_alleles.first) == 0)
				{
					forRet_allele_by_haplotype_P[levelI].first[s.haplotypes_alleles.first] = 0;
				}
				forRet_allele_by_haplotype_P.at(levelI).first.at(s.haplotypes_alleles.first) += oneSampleWeight;

				if(forRet_allele_by_haplotype_P[levelI].second.count(s.haplotypes_alleles.second) == 0)
				{
					forRet_allele_by_haplotype_P[levelI].second[s.haplotypes_alleles.second] = 0;
				}
				forRet_allele_by_haplotype_P.at(levelI).second.at(s.haplotypes_alleles.second) += oneSampleWeight;

				std::string copyingFromPlusAlleles = currentGene_MSA_int_2_id.at(s.copyingFrom.first) + "|" + currentGene_MSA_int_2_id.at(s.copyingFrom.second) + "/" + s.haplotypes_alleles.first + "|" + s.haplotypes_alleles.second;
				if(samples_copyingFromPlusAllele_by_position.at(levelI).count(copyingFromPlusAlleles) == 0)
					samples_copyingFromPlusAllele_by_position.at(levelI)[copyingFromPlusAlleles] = 0;

				samples_copyingFromPlusAllele_by_position.at(levelI).at(copyingFromPlusAlleles)++;


				size_t readAssignmentState_i = s.readAssignmentState;
				std::string readAssignmentState_str = readAssignmentStates.at(readAssignmentState_i);
				if(forRet_samples_readAssignmentStates != 0)
				{
					forRet_samples_readAssignmentStates->at(sampleI).at(levelI) = readAssignmentState_str;
				}
				if(readAssignment_2_activeReads.count(readAssignmentState_str) == 0)
				{
					std::cerr << "Missing info for string: " << readAssignmentState_str << "\n" << std::flush;
				}
				assert(readAssignment_2_activeReads.count(readAssignmentState_str));
				for(const auto activeRead : readAssignment_2_activeReads.at(readAssignmentState_str))
				{
					assert((activeRead.second == '1') || (activeRead.second == '2'));
					assert(read_genotypes_per_position.at(currentGene).count(activeRead.first));

					std::string readID = activeRead.first;
					std::string thisAssignment_thisRead_haplotype = (activeRead.second == '1') ? "H1" : "H2";

					if(thisSample_readID_2_haplotype.count(readID))
					{
						assert(thisSample_readID_2_haplotype.at(readID) == thisAssignment_thisRead_haplotype);
					}
					else
					{
						thisSample_readID_2_haplotype[readID] = thisAssignment_thisRead_haplotype;
					}


					for(const auto activeRead2 : readAssignment_2_activeReads.at(readAssignmentState_str))
					{
						std::string readID2 = activeRead2.first;
						std::string thisAssignment_thisRead2_haplotype = (activeRead2.second == '1') ? "H1" : "H2";

						if(readID < readID2)
						{
							std::pair<std::string, std::string> readPairID = std::make_pair(readID, readID2);
							overlappingReadIDPairs.insert(readPairID);

							if(readPairs_haplotypeIdentical[readPairID].count(levelI) == 0)
							{
								readPairs_haplotypeIdentical[readPairID][levelI] = std::make_pair(0,0);
							}

							if(thisAssignment_thisRead_haplotype == thisAssignment_thisRead2_haplotype)
							{
								readPairs_haplotypeIdentical[readPairID][levelI].first++;
							}
							else
							{
								readPairs_haplotypeIdentical[readPairID][levelI].second++;
							}
						}
					}
				}
			}
			outputStream_sampledCompleteHaplotypes << h1 << "\n" << h2 << "\n\n";

			for(auto readID2Haplotype : thisSample_readID_2_haplotype)
			{
				if(readID_2_haplotype[readID2Haplotype.first].count(readID2Haplotype.second) == 0)
				{
					readID_2_haplotype[readID2Haplotype.first][readID2Haplotype.second] = 0;
				}
				readID_2_haplotype.at(readID2Haplotype.first).at(readID2Haplotype.second)++;
			}
		}
	}

	auto printSortedVector = [](const std::vector<std::pair<std::string, double>>& v) -> std::string {
		std::string forReturn;
		for(auto e : v)
		{
			forReturn = forReturn + e.first + ":" + std::to_string(e.second) + " ";
		}
		return forReturn;
	};
	for(unsigned int levelI = 0; levelI < currentGene_geneLength; levelI++)
	{
		std::map<std::string, unsigned int> level_read_genotypes;
		for(auto oneRead : read_genotypes_per_position.at(currentGene))
		{
			if(oneRead.second.count(levelI))
			{
				std::string genotype = oneRead.second.at(levelI);
				if(level_read_genotypes.count(genotype) == 0)
				{
					level_read_genotypes[genotype] = 0;
				}
				level_read_genotypes.at(genotype)++;
			}
		}

		outputStream_sampledGenotypesAndCopiedFrom
			<< levelI
			<< "\t"
			<< printSortedVector(Utilities::map2Freq_sorted(samples_genotypes_by_position.at(levelI)))
			<< "\t" << printSortedVector(Utilities::map2Freq_sorted(samples_copyingFrom_by_position.at(levelI)))
			<< "\t" << printSortedVector(Utilities::map2Freq_sorted(samples_copyingFromPlusAllele_by_position.at(levelI)))
			<< "\t" << Utilities::join(currentGene_activeAlleles_perPosition.at(levelI), ";")
			<< "\t" << Utilities::JoinMapUInt2Str(level_read_genotypes)
			<< "\n";
	}

	for(auto readData : readID_2_haplotype)
	{
		unsigned int n_h1 = 0;
		if(readData.second.count("H1"))
			n_h1 = readData.second.at("H1");

		unsigned int n_h2 = 0;
		if(readData.second.count("H2"))
			n_h2 = readData.second.at("H2");

		assert((n_h1 + n_h2) > 0);
		double p_h1 = double(n_h1)/double(n_h1 + n_h2);

		forRet_oneReadP_h1[readData.first] = p_h1;

		outputStream_readID2Haplotype << readData.first << "\t" << printSortedVector(Utilities::map2Freq_sorted(readID_2_haplotype.at(readData.first))) << "\n";
	}

	for(const auto& readPair : overlappingReadIDPairs)
	{
		std::pair<unsigned int, unsigned int> allLevelData;
		bool firstIteration = true;
		for(const auto& levelData : readPairs_haplotypeIdentical.at(readPair))
		{
			unsigned int levelI = levelData.first;
			if(firstIteration)
			{
				allLevelData = levelData.second;
			}
			else
			{
				assert(allLevelData.first == levelData.second.first);
				assert(allLevelData.second == levelData.second.second);
			}
			firstIteration = false;
		}
		assert(firstIteration == false);
		outputStream_readIDSameHaplotype << readPair.first << " " << readPair.second << " " << allLevelData.first << " " << (allLevelData.first + allLevelData.second) << "\n";

		forRet_readPair_differentHaplotypes_P[readPair] = 1 - double(allLevelData.first)/double(allLevelData.first + allLevelData.second);
	}

	// auto states_min_max = std::minmax_element(states_per_position.begin(), states_per_position.end());
	// std::cout << "=====" << "\n" << "\tMin: " << *(states_min_max.first) << " - max: " << *(states_min_max.second) << "\n" << std::flush;

	std::cout << "Inference for " << geneID << " done.\n" << std::flush;

	return forward_joint_log_p;
}

std::vector<double> fullLengthHMM::computeEmissionProbabilities(size_t level) const
{
	std::map<std::pair<std::string, size_t>, double> readAssignment2P;
	std::vector<double> emission_P;
	emission_P.reserve(statesByLevel.at(level).size());
	for(size_t stateI = 0; stateI < statesByLevel.at(level).size(); stateI++)
	{
		const HMMstate& s = statesByLevel.at(level).at(stateI);
		std::string underlyingGenotype = s.haplotypes_alleles.first + "|" + s.haplotypes_alleles.second;
		if(readAssignment2P.count(std::make_pair(underlyingGenotype, s.readAssignmentState)))
		{
			emission_P.push_back(readAssignment2P.at(std::make_pair(underlyingGenotype, s.readAssignmentState)));
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
			size_t reads_h1 = 0;
			size_t reads_h2 = 0;
			size_t matches = 0;
			size_t mismatches = 0;
			std::string alleles_reads_h1;
			std::string alleles_reads_h2;

			for(auto activeReads : readAssignment_2_activeReads.at(readAssignmentString))
			{
				assert((activeReads.second == '1') || (activeReads.second == '2'));
				assert(read_genotypes_per_position.at(currentGene).count(activeReads.first));

				if(read_genotypes_per_position.at(currentGene).at(activeReads.first).count(level))
				{
					std::string readAllele = read_genotypes_per_position.at(currentGene).at(activeReads.first).at(level);
					std::string underlyingAllele = (activeReads.second == '1') ? s.haplotypes_alleles.first : s.haplotypes_alleles.second;
					emissionP *= ((readAllele == underlyingAllele) ? 0.98 : 0.02);

					if(activeReads.second == '1')
					{
						reads_h1++;
						alleles_reads_h1.append(readAllele);
					}
					else
					{
						reads_h2++;
						alleles_reads_h2.append(readAllele);
					}

					if(readAllele == underlyingAllele)
					{
						matches++;
					}
					else
					{
						mismatches++;
					}
					// std::cerr << "readAllele: " << readAllele << "; underlyingAllele: " << underlyingAllele << "; emissionP = " << ((readAllele == underlyingAllele) ? 0.99 : 0.01) << "\n";
				}
			}
			readAssignment2P[std::make_pair(underlyingGenotype, s.readAssignmentState)] = emissionP;
			emission_P.push_back(emissionP);

			
			if(0 && (level >= 1375) && (level <= 1375))
			{
				std::cout << "\t" << "State " << stateI << " at level " << level << "\n";
				std::cout << "\t\t" << "readAssignment_2_activeReads.at(readAssignmentString).size() " << ": " << readAssignment_2_activeReads.at(readAssignmentString).size() << "\n"; 
				std::cout << "\t\t" << "haplotypes_alleles" << ": " << s.haplotypes_alleles.first << " / " << s.haplotypes_alleles.second << "\n";
				std::cout << "\t\t" << "copyingFrom" << ": " << s.copyingFrom.first << " / " << s.copyingFrom.second << "\n";
				std::cout << "\t\t" << "readAssignmentState" << ": " << s.readAssignmentState << "\n";
				std::cout << "\t\t" << "emissionP" << ": " << emissionP<< "\n";	
				std::cout << "\t\t" << "Underlying genotype " <<  underlyingGenotype << ", read assignment state " << s.readAssignmentState << ": " << reads_h1 << " reads on h1, " << reads_h2 << " reads on h2, " << matches << " matches, " << mismatches << " mismatches, emission P: " << emissionP << "\n" << std::flush;
			}
		}
	}

	//assert(2 == 3);
	return emission_P;
}

std::vector<double> fullLengthHMM::computeInitialProbabilities() const
{
	std::vector<double> forReturn;

	size_t hom_states = (currentGene_haplotypeResolved) ? (pow(currentGene_MSA_ids_h1.size(), 2) + pow(currentGene_MSA_ids_h2.size(), 2)) : currentGene_MSA_ids.size();
	size_t het_states = pow(currentGene_MSA_ids.size(), 2) - hom_states;

	double het_haploGroups_p = 0.9;
	double hom_haploGroups_p = 1 - het_haploGroups_p;

	double hom_haploGroups_withinH1 = hom_haploGroups_p * 0.5 * (1.0 / pow(currentGene_MSA_ids_h1.size(), 2));
	double hom_haploGroups_withinH2 = hom_haploGroups_p * 0.5 * (1.0 / pow(currentGene_MSA_ids_h2.size(), 2));

	bool verbose = false;
	
	if(! currentGene_haplotypeResolved)
	{
		if(verbose)
			std::cerr << "(! currentGene_haplotypeResolved) active!\n" << std::flush;
	
		assert(currentGene_MSA_ids.size());
		hom_haploGroups_withinH1 = hom_haploGroups_p * 0.5 * (1.0 / (double)currentGene_MSA_ids.size());
		hom_haploGroups_withinH2 = hom_haploGroups_withinH1;
	} 
	
	if(verbose)
	{
		std::cerr << "currentGene_haplotypeResolved" << ": " << currentGene_haplotypeResolved << "\n";
		std::cerr << "currentGene_MSA_ids_h1.size()" << ": " << currentGene_MSA_ids_h1.size() << "\n";
		std::cerr << "currentGene_MSA_ids_h2.size()" << ": " << currentGene_MSA_ids_h2.size() << "\n";
		std::cerr << "currentGene_MSA_ids.size()" << ": " << currentGene_MSA_ids.size() << "\n";
		std::cerr << "hom_states" << ": " << hom_states << "\n";
		std::cerr << "het_states" << ": " << het_states << "\n";
		std::cerr << "het_haploGroups_p" << ": " << het_haploGroups_p << "\n";
		std::cerr << "hom_haploGroups_p" << ": " << hom_haploGroups_p << "\n";
		std::cerr << "hom_haploGroups_withinH1" << ": " << hom_haploGroups_withinH1 << "\n";
		std::cerr << "hom_haploGroups_withinH2" << ": " << hom_haploGroups_withinH2 << "\n";
		std::cerr << std::flush;
	}
	
	size_t readAssignmentStates_at_level_0 = level_readAssignmentState_2_states.at(0).size();
	unsigned int n_activeAlleles_thisLevel = currentGene_activeAlleles_perPosition.at(0).size();
	assert(n_activeAlleles_thisLevel > 0);

	forReturn.reserve(statesByLevel.at(0).size());
	double sum_p = 0;
	size_t hom_states_check = 0;
	for(size_t stateI = 0; stateI < statesByLevel.at(0).size(); stateI++)
	{
		const HMMstate& s = statesByLevel.at(0).at(stateI);

		double copyFromP;

		const std::string& copyFrom_h1_string = currentGene_MSA_int_2_HaplotypeID.at(s.copyingFrom.first);
		const std::string& copyFrom_h2_string = currentGene_MSA_int_2_HaplotypeID.at(s.copyingFrom.second);
		
		const std::string& copyFrom_h1_allele = currentGene_MSA_int_2_id.at(s.copyingFrom.first);
		const std::string& copyFrom_h2_allele = currentGene_MSA_int_2_id.at(s.copyingFrom.second);

		if(currentGene_haplotypeResolved)
		{
			if(copyFrom_h1_string == copyFrom_h2_string)
			{
				//hom_states_check++;
				if(currentGene_MSA_ids_h1.count(s.copyingFrom.first))
					copyFromP = hom_haploGroups_withinH1;
				else
					copyFromP = hom_haploGroups_withinH2;
			}
			else
			{
				if(s.copyingFrom.first < s.copyingFrom.second)
				{
					copyFromP = 2.0 * het_haploGroups_p * (1.0/double(het_states));
				}
				else
				{
					copyFromP = 0;
				}
			}
		}
		else
		{
			if(copyFrom_h1_allele == copyFrom_h2_allele)
			{
				copyFromP = 1.0/(double)pow(currentGene_MSA_ids.size(), 2);
			}
			else
			{
				if(s.copyingFrom.first < s.copyingFrom.second)
				{
					copyFromP = 2.0/(double)pow(currentGene_MSA_ids.size(), 2);
				}
				else
				{
					copyFromP = 0;
				}
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

		if(verbose)
		{
			std::cerr << "Initial probabilities state " << stateI << "\n";
			std::cerr << "\t" << "s.copyingFrom.first" << ": " << s.copyingFrom.first << " (" << copyFrom_h1_string << ", " << copyFrom_h1_allele << ")\n";
			std::cerr << "\t" << "s.copyingFrom.second" << ": " << s.copyingFrom.second << " (" << copyFrom_h2_string << ", " << copyFrom_h2_allele << ")\n";
			std::cerr << "\t" << "copyFromP" << ": " << copyFromP << "\n";
			std::cerr << "\t" << "alleles_match_h1" << ": " << alleles_match_h1 << "\n";
			std::cerr << "\t" << "alleles_match_h2" << ": " << alleles_match_h2 << "\n";
			std::cerr << "\t" << "allele_p" << ": " << allele_p << "\n";
			std::cerr << "\t" << "thisState_initialP" << ": " << thisState_initialP << "\n";
			std::cerr << "\n" << std::flush;
		}

	}
	/*
	std::cerr << "hom_states" << ": " << hom_states << "\n";
	std::cerr << "hom_states_check" << ": " << hom_states_check << "\n";
	std::cerr << "sum_p" << ": " << sum_p << "\n";
	std::cerr << std::flush;
	*/
	//assert(hom_states == hom_states_check);
	if(!(abs(1 - sum_p) <= 1e-4))
	{
		std::cerr << "sum_p: " << sum_p << "\n" << std::flush;
		std::cerr << "(abs(1 - sum_p) <= 1e-4): " << abs(1 - sum_p) << "\n" << std::flush;
	}
	assert(abs(1 - sum_p) <= 1e-4);

	assert(forReturn.size() == statesByLevel.at(0).size());
	return forReturn;
}

std::vector<HMMtransition> fullLengthHMM::computeLevelTransitions_backward(size_t next_level, long long limitToToState) const
{
	assert(currentGene.length());
	std::vector<HMMtransition> forReturn;

	// std::cerr << "Next_level: " << next_level << "\n";
	// std::cerr << "currentGene.length(): " << currentGene.length() << "\n" << std::flush;

	assert(next_level > 0);
	assert(next_level < currentGene_geneLength);

	std::set<size_t> new_read_IDs_indices_previousLevel;
	if(thisGene_reads_stop_per_position.at(currentGene).count(next_level-1))
	{
		for(auto readID : thisGene_reads_stop_per_position.at(currentGene).at(next_level-1))
		{
			size_t readID_index = readID_2_index.at(readID);
			new_read_IDs_indices_previousLevel.insert(readID_index);
		}
	}

	std::set<size_t> diseappearing_read_IDs_thisLevel;
	if(thisGene_reads_start_per_position.at(currentGene).count(next_level))
	{
		for(auto readID : thisGene_reads_start_per_position.at(currentGene).at(next_level))
		{
			size_t readID_index = readID_2_index.at(readID);
			diseappearing_read_IDs_thisLevel.insert(readID_index);
		}
	}
	
	std::set<size_t> diseappearing_read_IDs_previousLevel;
	if(thisGene_reads_stop_per_position.at(currentGene).count(next_level-1))
	{
		for(auto readID : thisGene_reads_stop_per_position.at(currentGene).at(next_level-1))
		{
			size_t readID_index = readID_2_index.at(readID);
			diseappearing_read_IDs_previousLevel.insert(readID_index);
		}
	}

	std::set<size_t> new_read_IDs_indices_thisLevel;
	if(thisGene_reads_start_per_position.at(currentGene).count(next_level))
	{
		for(auto readID : thisGene_reads_start_per_position.at(currentGene).at(next_level))
		{
			size_t readID_index = readID_2_index.at(readID);
			new_read_IDs_indices_thisLevel.insert(readID_index);
		}
	}	


	if(limitToToState == -1)
	{
		forReturn.reserve(statesByLevel.at(next_level-1).size() * 15);
	}
	else
	{
		forReturn.reserve(15);
	}

	size_t n_jumps = 0;
	std::map<size_t, std::set<size_t>> backward_read_assignment_states;
	std::map<size_t, std::set<size_t>> forward_read_assignment_states;
	// double forward_readState_transisition_p = diseappearing_read_IDs_thisLevel.size() ? 1.0/pow(2, (double)diseappearing_read_IDs_thisLevel.size()) : 1;

	const std::set<char>& activeAlleles_length1_thisLevel = currentGene_activeAlleles_perPosition_allelesLength1_chars.at(next_level);
	unsigned int n_activeAlleles_thisLevel = currentGene_activeAlleles_perPosition.at(next_level).size();
	assert(n_activeAlleles_thisLevel > 0);
	
	auto computeOneToState = [&](size_t stateI) -> void {	
		const HMMstate& next_s = statesByLevel.at(next_level).at(stateI);

		if(backward_read_assignment_states.count(next_s.readAssignmentState) == 0)
		{
			backward_read_assignment_states[next_s.readAssignmentState] = std::set<size_t>();
			if(constrainedReadAssignmentStates)
			{
				assert(constrainedReadAssignmentStates_transitions_backward.count(next_s.readAssignmentState));
				std::cout << "a" << std::flush;
				size_t transition_back = constrainedReadAssignmentStates_transitions_backward.at(next_s.readAssignmentState);
				std::cout << "b" << std::flush;
				backward_read_assignment_states[next_s.readAssignmentState] = {transition_back};
			}
			else
			{
				std::set<std::string> backwardAssignmentStates = previousLevel_compatibleReadAssignments(readAssignmentStates.at(next_s.readAssignmentState), new_read_IDs_indices_previousLevel, diseappearing_read_IDs_thisLevel);;
				for(auto const& oneAssingmentState : backwardAssignmentStates)
				{
					backward_read_assignment_states[next_s.readAssignmentState].insert(readAssignmentState_2_index.at(oneAssingmentState));
				}
			}
		}

		if((! constrainedReadAssignmentStates) && (forward_read_assignment_states.count(next_s.readAssignmentState) == 0))
		{
			// if(constrainedReadAssignmentStates)
			// { 
				// assert(constrainedReadAssignmentStates_transitions_forward.count(next_s.readAssignmentState));
				// forward_read_assignment_states[next_s.readAssignmentState] = {constrainedReadAssignmentStates_transitions_forward.at(next_s.readAssignmentState)};
			// }
			// else
			// { 
				std::set<std::string> forwardAssignmentStates = nextLevel_compatibleReadAssignments(readAssignmentStates.at(next_s.readAssignmentState), diseappearing_read_IDs_previousLevel, new_read_IDs_indices_thisLevel);
				for(auto const& oneAssingmentState : forwardAssignmentStates)
				{
					forward_read_assignment_states[next_s.readAssignmentState].insert(readAssignmentState_2_index.at(oneAssingmentState));
				}
			// }
		}
		
		//double _sum_jump_Ps = 0;

		// double read_assignment_state_p = forward_readState_transisition_p;
		double read_assignment_state_p = (constrainedReadAssignmentStates) ? 1 : (1.0/double(forward_read_assignment_states.at(next_s.readAssignmentState).size()));

		assert(backward_read_assignment_states.count(next_s.readAssignmentState));
		for(size_t readAssignmentState_previousLevel_index : backward_read_assignment_states.at(next_s.readAssignmentState))
		{
			// size_t readAssignmentState_previousLevel_index = readAssignmentState_2_index.at(readAssignmentState_previousLevel);
			if(level_readAssignmentState_2_states.at(next_level-1).count(readAssignmentState_previousLevel_index) == 0)
			{
				std::cerr << "Missing a read assingment state at the previous level.\n";
				std::cerr << "\t" << "next_level" << ": " << next_level << "\n";
				std::cerr << "\t" << "next_level-1" << ": " << next_level-1 << "\n";
				std::cerr << "\t" << "readAssignmentState_previousLevel_index" << ": " << readAssignmentState_previousLevel_index << "\n";
				std::cerr << std::flush;
			}
			
			assert(level_readAssignmentState_2_states.at(next_level-1).count(readAssignmentState_previousLevel_index));
			std::set<size_t> previousLevel_jumpFromStates = level_readAssignmentState_2_states.at(next_level-1).at(readAssignmentState_previousLevel_index);

			for(auto previousStateIndex : previousLevel_jumpFromStates)
			{				
				const HMMstate& previous_s = statesByLevel.at(next_level-1).at(previousStateIndex);

				double haplotype_copy_p = 1;

				if(currentGene_haplotypeResolved)
				{
					if(previous_s.copyingFrom.first != next_s.copyingFrom.first)
					{
						if(currentGene_MSA_ids_same_haploGroup.count(std::make_pair(previous_s.copyingFrom.first, next_s.copyingFrom.first)))
						{
							haplotype_copy_p *= (currentGene_MSA_ids_h1.count(previous_s.copyingFrom.first) ? tr_haplotypeResolution_change_h1_withinGroup : tr_haplotypeResolution_change_h2_withinGroup);
						}
						else
						{
							haplotype_copy_p *= (currentGene_MSA_ids_h1.count(previous_s.copyingFrom.first) ?  tr_haplotypeResolution_change_h1_outOfGroup : tr_haplotypeResolution_change_h2_outOfGroup);
						}
					}
					else
					{
						haplotype_copy_p *=  (currentGene_MSA_ids_h1.count(previous_s.copyingFrom.first) ? tr_haplotypeResolution_remain_h1 : tr_haplotypeResolution_remain_h2);
					}


					if(previous_s.copyingFrom.second != next_s.copyingFrom.second)
					{
						if(currentGene_MSA_ids_same_haploGroup.count(std::make_pair(previous_s.copyingFrom.second, next_s.copyingFrom.second)))
						{
							haplotype_copy_p *= (currentGene_MSA_ids_h1.count(previous_s.copyingFrom.second) ? tr_haplotypeResolution_change_h1_withinGroup : tr_haplotypeResolution_change_h2_withinGroup);
						}
						else
						{
							haplotype_copy_p *= (currentGene_MSA_ids_h1.count(previous_s.copyingFrom.second) ?  tr_haplotypeResolution_change_h1_outOfGroup : tr_haplotypeResolution_change_h2_outOfGroup);
						}
					}
					else
					{
						haplotype_copy_p *=  (currentGene_MSA_ids_h1.count(previous_s.copyingFrom.second) ? tr_haplotypeResolution_remain_h1 : tr_haplotypeResolution_remain_h2);
					}
				}
				else
				{
					bool isSwitch_h1 = (previous_s.copyingFrom.first != next_s.copyingFrom.first);
					bool isSwitch_h2 = (previous_s.copyingFrom.second != next_s.copyingFrom.second);
					haplotype_copy_p =
								(isSwitch_h1 ? tr_noHaplotypeResolution_change_oneh : tr_noHaplotypeResolution_remain_oneh) *
								(isSwitch_h2 ? tr_noHaplotypeResolution_change_oneh : tr_noHaplotypeResolution_remain_oneh);
				}

				const std::string& activeAllele_h1 = next_s.haplotypes_alleles.first;
				const std::string& activeAllele_h2 = next_s.haplotypes_alleles.second;
				assert(MSA_reference_sequences.count(currentGene));
				assert(MSA_reference_sequences.at(currentGene).count(currentGene_MSA_int_2_id.at(next_s.copyingFrom.first)));
				assert(MSA_reference_sequences.at(currentGene).count(currentGene_MSA_int_2_id.at(next_s.copyingFrom.second)));
				const char& copiedAlleleFrom_h1 = MSA_reference_sequences.at(currentGene).at(currentGene_MSA_int_2_id.at(next_s.copyingFrom.first)).at(next_level);
				const char& copiedAlleleFrom_h2 = MSA_reference_sequences.at(currentGene).at(currentGene_MSA_int_2_id.at(next_s.copyingFrom.second)).at(next_level);
				
				bool copiedAlleleFrom_h1_exists_in_activeAlleles = (activeAlleles_length1_thisLevel.count(copiedAlleleFrom_h1) > 0);
				bool copiedAlleleFrom_h2_exists_in_activeAlleles = (activeAlleles_length1_thisLevel.count(copiedAlleleFrom_h2) > 0);
				
				bool alleles_match_h1 = (activeAllele_h1.length() == 1) && (activeAllele_h1.at(0) == copiedAlleleFrom_h1);
				bool alleles_match_h2 = (activeAllele_h2.length() == 1) && (activeAllele_h2.at(0) == copiedAlleleFrom_h2);
				
				double local_currentGene_novel_allele_p_h1 = currentGene_novel_allele_p;
				double local_currentGene_novel_allele_p_h2 = currentGene_novel_allele_p;
				
				if(! copiedAlleleFrom_h1_exists_in_activeAlleles)
					local_currentGene_novel_allele_p_h1 = 1;
				if(! copiedAlleleFrom_h2_exists_in_activeAlleles)
					local_currentGene_novel_allele_p_h2 = 1;
				
				double allele_p =
						(alleles_match_h1 ? ((1 - local_currentGene_novel_allele_p_h1) + (local_currentGene_novel_allele_p_h1 * (1.0 / n_activeAlleles_thisLevel))) : (local_currentGene_novel_allele_p_h1 * (1.0 / n_activeAlleles_thisLevel))) *
						(alleles_match_h2 ? ((1 - local_currentGene_novel_allele_p_h2) + (local_currentGene_novel_allele_p_h2 * (1.0 / n_activeAlleles_thisLevel))) : (local_currentGene_novel_allele_p_h2 * (1.0 / n_activeAlleles_thisLevel)));

				double p_jump = read_assignment_state_p * haplotype_copy_p * allele_p;

				HMMtransition thisTransition;
				thisTransition.from_level = next_level - 1;
				thisTransition.from_state = previousStateIndex;
				thisTransition.to_state = stateI;
				thisTransition.P = p_jump;
				forReturn.push_back(thisTransition);

				// _sum_jump_Ps += p_jump;
				n_jumps++;
			}
		}

		// std::cerr << "Level " << first_level << " state " << stateI << " sum of transition (J)  Ps: " << _sum_jump_Ps << "\n" << std::flush;
		// assert(abs(1 - _sum_jump_Ps) <= 1e-4);
	};

	if(limitToToState == -1)
	{
		for(size_t stateI = 0; stateI < statesByLevel.at(next_level).size(); stateI++)
		{
			computeOneToState(stateI);
		}
	}
	else
	{
		assert((limitToToState >= 0) && (limitToToState < statesByLevel.at(next_level).size()));
		computeOneToState(limitToToState);
	}

	return forReturn;
}

std::vector<HMMtransition> fullLengthHMM::computeLevelTransitions(size_t first_level, long long limitToFromState) const
{
	assert(currentGene.length());
	std::vector<HMMtransition> forReturn;

	// std::cerr << "\t" << "limitToFromState" << ": " << limitToFromState << "\n" << std::flush;
	// std::cerr << "\t" << "first_level" << ": " << first_level << "\n" << std::flush;

	std::set<size_t> diseappearing_read_IDs_indices;
	std::set<size_t> new_read_IDs_indices_nextLevel;

	if(! constrainedReadAssignmentStates)
	{
		if(thisGene_reads_stop_per_position.at(currentGene).count(first_level))
		{
			for(auto readID : thisGene_reads_stop_per_position.at(currentGene).at(first_level))
			{
				size_t readID_index = readID_2_index.at(readID);
				diseappearing_read_IDs_indices.insert(readID_index);
			}
		}

		if(thisGene_reads_start_per_position.at(currentGene).count(first_level+1))
		{
			for(auto readID : thisGene_reads_start_per_position.at(currentGene).at(first_level+1))
			{
				size_t readID_index = readID_2_index.at(readID);
				new_read_IDs_indices_nextLevel.insert(readID_index);
			}
		}
	}

	if(limitToFromState == -1)
	{
		forReturn.reserve(statesByLevel.at(first_level).size() * 15);
	}
	else
	{
		forReturn.reserve(15);
	}
	

	size_t n_jumps = 0;
	std::map<size_t, std::set<size_t>> forward_read_assignment_states;

	// double forward_readState_transisition_p = new_read_IDs_indices_nextLevel.size() ? 1.0/pow(2, (double)new_read_IDs_indices_nextLevel.size()) : 1;

	unsigned int n_activeAlleles_nextLevel = currentGene_activeAlleles_perPosition.at(first_level+1).size();
	assert(n_activeAlleles_nextLevel > 0);
	
	const std::set<char>& activeAlleles_length1_nextLevel = currentGene_activeAlleles_perPosition_allelesLength1_chars.at(first_level+1);
	
	auto computeOneFromState = [&](size_t stateI) -> void {
		
		bool verbose = ((first_level == 0) && (stateI == 0));
		
		std::stringstream debOut;
		
		const HMMstate& s = statesByLevel.at(first_level).at(stateI);

		if(forward_read_assignment_states.count(s.readAssignmentState) == 0)
		{
			if(constrainedReadAssignmentStates)
			{
				forward_read_assignment_states[s.readAssignmentState] = {constrainedReadAssignmentStates_transitions_forward.at(s.readAssignmentState)};
			}
			else
			{
				std::set<std::string> forwardAssignmentStates = nextLevel_compatibleReadAssignments(readAssignmentStates.at(s.readAssignmentState), diseappearing_read_IDs_indices, new_read_IDs_indices_nextLevel);
				for(auto const& oneAssingmentState : forwardAssignmentStates)
				{
					forward_read_assignment_states[s.readAssignmentState].insert(readAssignmentState_2_index.at(oneAssingmentState));
				}
			}
		}

		if(verbose)
		{
			debOut << "Computing outgoing transitions from state " << stateI << "\n";
			debOut << "\t" << "forward_read_assignment_states.at(s.readAssignmentState).size()" << ": " << forward_read_assignment_states.at(s.readAssignmentState).size() << "\n";	
			debOut << "\t" << "tr_noHaplotypeResolution_change_oneh" << ": " << tr_noHaplotypeResolution_change_oneh << "\n";	
			debOut << "\t" << "tr_noHaplotypeResolution_remain_oneh" << ": " << tr_noHaplotypeResolution_remain_oneh << "\n";	
			debOut << "\t" << "tr_change_oneh" << ": " << tr_change_oneh << "\n";	
		}
		
		double _sum_jump_Ps = 0;


		double read_assignment_state_p = 1.0 / forward_read_assignment_states.at(s.readAssignmentState).size();
		/*
		if(!(abs(read_assignment_state_p - forward_readState_transisition_p) <= 1e-5))
		{
			debOut << "(abs(read_assignment_state_p - forward_readState_transisition_p) <= 1e-5)!" << "\n";
			debOut << "new_read_IDs_indices_nextLevel.size(): " << new_read_IDs_indices_nextLevel.size() << "\n";
			debOut << "forward_read_assignment_states.at(s.readAssignmentState).size(): " << forward_read_assignment_states.at(s.readAssignmentState).size() << "\n";
			debOut << "read_assignment_state_p: " << read_assignment_state_p << "\n";
			debOut << "forward_readState_transisition_p: " << forward_readState_transisition_p << "\n";
		}
		assert(abs(read_assignment_state_p - forward_readState_transisition_p) <= 1e-5);
		*/
		
		if(verbose)
		{
			debOut << "\t" << "n_activeAlleles_nextLevel" << ": " << n_activeAlleles_nextLevel << "\n";	
		}
		
		for(size_t readAssignmentState_nextLevel_index : forward_read_assignment_states.at(s.readAssignmentState))
		{
			if(level_readAssignmentState_2_states.at(first_level+1).count(readAssignmentState_nextLevel_index) == 0)
			{
				debOut<< "Missing a read assingment state at the next level.\n";
				debOut << "\t" << "first_level" << ": " << first_level << "\n";
				debOut << "\t" << "first_level+1" << ": " << first_level+1 << "\n";
				debOut << "\t" << "readAssignmentState_nextLevel_index" << ": " << readAssignmentState_nextLevel_index << "\n";
				debOut << "\t" << "readAssignmentState_nextLevel_index" << ": " << readAssignmentState_nextLevel_index << "\n";
				debOut << std::flush;
			}
			std::set<size_t> nextLevel_jumpIntoStates = level_readAssignmentState_2_states.at(first_level+1).at(readAssignmentState_nextLevel_index);

			if(verbose)
			{
				debOut << "\t\t" << "Read assignment state " << readAssignmentState_nextLevel_index << "\n";	
			}
			
			for(auto nextStateIndex : nextLevel_jumpIntoStates)
			{
				const HMMstate& next_s = statesByLevel.at(first_level+1).at(nextStateIndex);

				if(verbose)
				{
					debOut << "\t\t\t" << "State " << nextStateIndex << "\n";	
				}
				
				double haplotype_copy_p = 1;

				if(verbose)
				{
					debOut << "\t\t\t\t" << "Transition computation" << "\n";
				}
				
				if(currentGene_haplotypeResolved)
				{
					if(verbose)
					{
						debOut << "\t\t\t\t\t" << "currentGene_haplotypeResolved = 1" << "\n";
						debOut << "\t\t\t\t\t" << "s.copyingFrom.first" << ": " << s.copyingFrom.first << "\n";
						debOut << "\t\t\t\t\t" << "next_s.copyingFrom.first" << ": " << next_s.copyingFrom.first << "\n";
					}
					
					double haplotype_copy_p_h1 = 1;
					if(s.copyingFrom.first != next_s.copyingFrom.first)
					{
						if(verbose)
							debOut << "\t\t\t\t\t\t" << "s.copyingFrom.first != next_s.copyingFrom.first\n";	
						
						if(currentGene_MSA_ids_same_haploGroup.count(std::make_pair(s.copyingFrom.first, next_s.copyingFrom.first)))
						{
							if(verbose)
								debOut << "\t\t\t\t\t\t\t" << "currentGene_MSA_ids_same_haploGroup = 1\n";	
							
							haplotype_copy_p_h1 *= (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ? tr_haplotypeResolution_change_h1_withinGroup : tr_haplotypeResolution_change_h2_withinGroup);
						}
						else
						{
							if(verbose)
								debOut << "\t\t\t\t\t\t\t" << "currentGene_MSA_ids_same_haploGroup = 0\n";
							
							haplotype_copy_p_h1 *= (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ?  tr_haplotypeResolution_change_h1_outOfGroup : tr_haplotypeResolution_change_h2_outOfGroup);
						}
					}
					else
					{
						if(verbose)
							debOut << "\t\t\t\t\t\t" << "s.copyingFrom.first == next_s.copyingFrom.first\n";												
						
						haplotype_copy_p_h1 *=  (currentGene_MSA_ids_h1.count(s.copyingFrom.first) ? tr_haplotypeResolution_remain_h1 : tr_haplotypeResolution_remain_h2);
						
						if(verbose)
							debOut << "\t\t\t\t\t\t" << "currentGene_MSA_ids_h1.count(s.copyingFrom.first" << ": " << currentGene_MSA_ids_h1.count(s.copyingFrom.first) << "\n";																		
					}
					if(verbose)
					{
						debOut << "\t\t\t\t\t\t" << "haplotype_copy_p_h1" << ": " << haplotype_copy_p_h1 << "\n";
					}

					double haplotype_copy_p_h2 = 1;
					if(s.copyingFrom.second != next_s.copyingFrom.second)
					{
						if(verbose)
							debOut << "\t\t\t\t\t\t" << "s.copyingFrom.second != next_s.copyingFrom.second\n";	
												
						if(currentGene_MSA_ids_same_haploGroup.count(std::make_pair(s.copyingFrom.second, next_s.copyingFrom.second)))
						{
							if(verbose)
								debOut << "\t\t\t\t\t\t\t" << "currentGene_MSA_ids_same_haploGroup = 1\n";	
														
							haplotype_copy_p_h2 *= (currentGene_MSA_ids_h1.count(s.copyingFrom.second) ? tr_haplotypeResolution_change_h1_withinGroup : tr_haplotypeResolution_change_h2_withinGroup);
						}
						else
						{
							if(verbose)
								debOut << "\t\t\t\t\t\t\t" << "currentGene_MSA_ids_same_haploGroup = 0\n";							
							haplotype_copy_p_h2 *= (currentGene_MSA_ids_h1.count(s.copyingFrom.second) ?  tr_haplotypeResolution_change_h1_outOfGroup : tr_haplotypeResolution_change_h2_outOfGroup);
						}
					}
					else
					{
						if(verbose)
							debOut << "\t\t\t\t\t\t" << "s.copyingFrom.second == next_s.copyingFrom.second\n";												
												
						haplotype_copy_p_h2 *=  (currentGene_MSA_ids_h1.count(s.copyingFrom.second) ? tr_haplotypeResolution_remain_h1 : tr_haplotypeResolution_remain_h2);
						
						if(verbose)
							debOut << "\t\t\t\t\t\t" << "currentGene_MSA_ids_h1.count(s.copyingFrom.second)" << ": " << currentGene_MSA_ids_h1.count(s.copyingFrom.second) << "\n";																								
					}
					
					if(verbose)
					{
						debOut << "\t\t\t\t\t\t" << "haplotype_copy_p_h2" << ": " << haplotype_copy_p_h2 << "\n";
					}
					
					haplotype_copy_p = (haplotype_copy_p_h1 * haplotype_copy_p_h2);
					
					if(verbose)
					{
						debOut << "\t\t\t\t\t" << "haplotype_copy_p" << ": " << haplotype_copy_p << "\n";
						debOut << "\t\t\t\t\t\t" << "haplotype_copy_p_h1" << ": " << haplotype_copy_p_h1 << "\n";
						debOut << "\t\t\t\t\t\t" << "haplotype_copy_p_h2" << ": " << haplotype_copy_p_h2 << "\n";
					}					
				}
				else
				{
					bool isSwitch_h1 = (s.copyingFrom.first != next_s.copyingFrom.first);
					bool isSwitch_h2 = (s.copyingFrom.second != next_s.copyingFrom.second);
					haplotype_copy_p =
								(isSwitch_h1 ? tr_noHaplotypeResolution_change_oneh : tr_noHaplotypeResolution_remain_oneh) *
								(isSwitch_h2 ? tr_noHaplotypeResolution_change_oneh : tr_noHaplotypeResolution_remain_oneh);
				}

				const std::string& activeAllele_h1 = next_s.haplotypes_alleles.first;
				const std::string& activeAllele_h2 = next_s.haplotypes_alleles.second;
				const char& copiedAlleleFrom_h1 = MSA_reference_sequences.at(currentGene).at(currentGene_MSA_int_2_id.at(next_s.copyingFrom.first)).at(first_level+1);
				const char& copiedAlleleFrom_h2 = MSA_reference_sequences.at(currentGene).at(currentGene_MSA_int_2_id.at(next_s.copyingFrom.second)).at(first_level+1);
				
				bool copiedAlleleFrom_h1_exists_in_activeAlleles = (activeAlleles_length1_nextLevel.count(copiedAlleleFrom_h1) > 0);
				bool copiedAlleleFrom_h2_exists_in_activeAlleles = (activeAlleles_length1_nextLevel.count(copiedAlleleFrom_h2) > 0);
				
				bool alleles_match_h1 = (activeAllele_h1.length() == 1) && (activeAllele_h1.at(0) == copiedAlleleFrom_h1);
				bool alleles_match_h2 = (activeAllele_h2.length() == 1) && (activeAllele_h2.at(0) == copiedAlleleFrom_h2);
				
				double local_currentGene_novel_allele_p_h1 = currentGene_novel_allele_p;
				double local_currentGene_novel_allele_p_h2 = currentGene_novel_allele_p;
				
				if(! copiedAlleleFrom_h1_exists_in_activeAlleles)
					local_currentGene_novel_allele_p_h1 = 1;
				if(! copiedAlleleFrom_h2_exists_in_activeAlleles)
					local_currentGene_novel_allele_p_h2 = 1;
				
				double allele_p =
						(alleles_match_h1 ? ((1 - local_currentGene_novel_allele_p_h1) + (local_currentGene_novel_allele_p_h1 * (1.0 / n_activeAlleles_nextLevel))) : (local_currentGene_novel_allele_p_h1 * (1.0 / n_activeAlleles_nextLevel))) *
						(alleles_match_h2 ? ((1 - local_currentGene_novel_allele_p_h2) + (local_currentGene_novel_allele_p_h2 * (1.0 / n_activeAlleles_nextLevel))) : (local_currentGene_novel_allele_p_h2 * (1.0 / n_activeAlleles_nextLevel)));

				if(verbose)
				{
					debOut << "\t\t\t\t" << "read_assignment_state_p " << read_assignment_state_p << "\n";	
					debOut << "\t\t\t\t" << "allele_p " << allele_p << "\n";	
					debOut << "\t\t\t\t\t" << "next_s.haplotypes_alleles.first " << next_s.haplotypes_alleles.first << "\n";	
					debOut << "\t\t\t\t\t" << "next_s.haplotypes_alleles.second " << next_s.haplotypes_alleles.second << "\n";	
					debOut << "\t\t\t\t\t" << "MSA_reference_sequences.at(currentGene).at(currentGene_MSA_int_2_id.at(next_s.copyingFrom.first)).at(first_level+1)) " << MSA_reference_sequences.at(currentGene).at(currentGene_MSA_int_2_id.at(next_s.copyingFrom.first)).at(first_level+1) << "\n";	
					debOut << "\t\t\t\t\t" << "MSA_reference_sequences.at(currentGene).at(currentGene_MSA_int_2_id.at(next_s.copyingFrom.second)).at(first_level+1)) " << MSA_reference_sequences.at(currentGene).at(currentGene_MSA_int_2_id.at(next_s.copyingFrom.second)).at(first_level+1) << "\n";	
					debOut << "\t\t\t\t\t" << "This state (s.copyingFrom.first / s.copyingFrom.second) " << s.copyingFrom.first << " / " << s.copyingFrom.second << "\n";	
					debOut << "\t\t\t\t\t" << "Next state (next_s.copyingFrom.first / next_s.copyingFrom.second) " << next_s.copyingFrom.first << " / " << next_s.copyingFrom.second << "\n";	
					debOut << "\t\t\t\t\t" << "alleles_match_h1 " << alleles_match_h1 << "\n";	
					debOut << "\t\t\t\t\t" << "alleles_match_h2 " << alleles_match_h2 << "\n";	
				}
				
				double p_jump = read_assignment_state_p * haplotype_copy_p * allele_p;

				if(verbose)
				{
					debOut << "\t\t\t\t" << "p_jump " << p_jump << "\n";	
				}
				
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
		if(!(abs(1 - _sum_jump_Ps) <= 1e-4))
		{
			std::cerr << "_sum_jump_Ps: " << _sum_jump_Ps << "\n";
			std::cerr << debOut.str() << "\n";
			std::cerr << "first_level " << first_level << "\n";
			std::cerr << "stateI " << stateI << "\n";
			std::cerr << std::flush;
		}
		assert(abs(1 - _sum_jump_Ps) <= 1e-4);
	};

	if(limitToFromState == -1)
	{
		for(size_t stateI = 0; stateI < statesByLevel.at(first_level).size(); stateI++)
		{
			//std::cerr << "\t\t" << "stateI" << ": " << stateI << "\n" << std::flush;
			computeOneFromState(stateI);
			//std::cerr << "\t\t" << "forReturn.size()" << ": " << forReturn.size() << "\n" << std::flush;
		}
	}
	else
	{
		assert((limitToFromState >= 0) && (limitToFromState < statesByLevel.at(first_level).size()));
		computeOneFromState(limitToFromState);
	}

	return forReturn;
}

std::pair<std::string, std::string> fullLengthHMM::makeGt(std::string allele1, std::string allele2)
{
	if(allele1 < allele2)
	{
		return std::make_pair(allele1, allele2);
	}
	else
	{
		return std::make_pair(allele2, allele1);
	}
}

void fullLengthHMM::removeActiveAllele(const std::string& geneID, unsigned int position, const std::string& allele)
{
	assert(activeAlleles_per_position.at(geneID).at(position).count(allele));
	activeAlleles_per_position.at(geneID).at(position).erase(allele);
	assert(activeAlleles_per_position.at(geneID).size() >= 1);
}

std::map<unsigned int, std::set<std::string>> fullLengthHMM::getActiveAllelesForGene(const std::string& geneID)
{
	return activeAlleles_per_position.at(geneID);
}

