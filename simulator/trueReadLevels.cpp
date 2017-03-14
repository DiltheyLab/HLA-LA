/*
 * trueReadLevels.cpp
 *
 *  Created on: 01.10.2015
 *      Author: AlexanderDilthey
 */

#include "trueReadLevels.h"

#include <fstream>
#include <iostream>
#include <assert.h>
#include "../Utilities.h"
#include <algorithm>

namespace simulator {

std::pair<size_t, size_t> trueReadLevels::evaluateAlignment(const mapper::reads::verboseSeedChainPair& alignment, const mapper::reads::oneReadPair& rP, mapper::aligner::extensionAligner* eA)
{
	std::pair<size_t, size_t> forReturn;
	forReturn.first = 0;
	forReturn.second = 0;

	if(!trueLevels.count(alignment.readID))
	{
		std::cerr << "! trueLevels.count(alignment.readID)" << "\n";
		std::cerr << "alignment.readID" << ": " << alignment.readID << "\n";
		std::cerr << "Existing IDs:\n";
		int rI = 0;
		for(std::map<std::string, std::pair<std::vector<int>, std::vector<int>>>::iterator readIDit = trueLevels.begin(); readIDit != trueLevels.end(); readIDit++)
		{
			std::cerr << "\t" << readIDit->first << "\n";
			rI++;
			if(rI > 5)
				break;
		}
		std::cout << std::flush;
	}
	assert(trueLevels.count(alignment.readID));
	assert(trueLevels.count(alignment.readID));
	assert(true_fullAlignment_coordinates_edgePath.count(alignment.readID));

	auto evaluateLevels = [&](bool first) -> void {
		const mapper::reads::verboseSeedChain& seedChain = ( (first) ? alignment.chains.first : alignment.chains.second );
		const mapper::reads::oneRead& read = ( (first) ? rP.reads.first : rP.reads.second );

		const std::string& aligned_sequence = ( (first) ? alignment.chains.first.sequence_aligned : alignment.chains.second.sequence_aligned );
		const std::vector<int>& aligned_graph_levels = ( (first) ? alignment.chains.first.graph_aligned_levels : alignment.chains.second.graph_aligned_levels );
		assert(aligned_sequence.size() == aligned_graph_levels.size());

		std::vector<int> graph_aligned_levels_withSequenceCharacters;
		for(unsigned int i = 0; i < aligned_sequence.length(); i++)
		{
			if(aligned_sequence.at(i) != '_')
			{
				graph_aligned_levels_withSequenceCharacters.push_back(aligned_graph_levels.at(i));
			}
		}

		std::vector<int> true_levels = ( (first) ? trueLevels.at(alignment.readID).first : trueLevels.at(alignment.readID).second );
		assert(true_levels.size() == graph_aligned_levels_withSequenceCharacters.size());

		if((first && alignment.chains.first.reverse) || ((!first) && alignment.chains.second.reverse))
		{
			std::reverse(true_levels.begin(), true_levels.end());
		}

		int thisAlignment_L = 0;
		int thisAlignment_Q = 0;
		for(unsigned int i = 0; i < true_levels.size(); i++)
		{
			forReturn.first++;
			thisAlignment_L++;
			if(graph_aligned_levels_withSequenceCharacters.at(i) == true_levels.at(i))
			{
				forReturn.second++;
				thisAlignment_Q++;
			}
		}

		double prop_OK = (double)thisAlignment_Q/(double)thisAlignment_L;
		//if(thisAlignment_L != thisAlignment_Q)
		if(prop_OK < 0.9)
		{
			double Score = eA->scoreOneAlignment(seedChain, read);

			std::vector<int> true_fa_coordinates_edgePath = ( (first) ? true_fullAlignment_coordinates_edgePath.at(alignment.readID).first : true_fullAlignment_coordinates_edgePath.at(alignment.readID).second );
			std::string true_fa_underlyingEdgelabels = ( (first) ? true_fullAlignment_underlyingEdgeLabels.at(alignment.readID).first : true_fullAlignment_underlyingEdgeLabels.at(alignment.readID).second );
			std::string true_fa_sequence = ( (first) ? true_fullAlignment_sequence.at(alignment.readID).first : true_fullAlignment_sequence.at(alignment.readID).second );

			mapper::reads::oneRead originalRead_alignmentOrientation = read;
			if((first && alignment.chains.first.reverse) || ((!first) && alignment.chains.second.reverse))
			{
				std::reverse(true_fa_coordinates_edgePath.begin(), true_fa_coordinates_edgePath.end());
				std::reverse(true_fa_underlyingEdgelabels.begin(), true_fa_underlyingEdgelabels.end());
				std::reverse(true_fa_sequence.begin(), true_fa_sequence.end());
				originalRead_alignmentOrientation.invert();
			}

			std::cout << "Evaluation " << alignment.readID << " -- levels correct: " << thisAlignment_Q << " of " << thisAlignment_L << "\n";
			std::cout << ((first) ? "\tFirst" : "\tSecond") << "\n";
			std::cout << "\t\tReverse: " << ((first) ? alignment.chains.first.reverse : alignment.chains.second.reverse ) << "\n";
			std::cout << "\t\t" << Utilities::join(Utilities::ItoStr(true_levels), ", ") << "\n";
			std::cout << "\t\t" << Utilities::join(Utilities::ItoStr(graph_aligned_levels_withSequenceCharacters), ", ") << "\n";
			std::cout << "\t\t Inferred alignment:\n";
			std::cout << "\t\t\t" << Utilities::join(Utilities::ItoStr(seedChain.graph_aligned_levels), ", ") << "\n";
			std::cout << "\t\t\t" << Utilities::join(Utilities::BtoStr(seedChain.is_from_BWAseed), "") << "\n";
			std::cout << "\t\t\t" << seedChain.graph_aligned << "\n";
			std::cout << "\t\t\t" << seedChain.sequence_aligned << "\n";

			std::cout << "\t\t" << "Score" << ": " << Score << "\n";

			assert(true_fa_coordinates_edgePath.size() == true_fa_underlyingEdgelabels.size());
			assert(true_fa_coordinates_edgePath.size() == true_fa_sequence.size());

			mapper::reads::verboseSeedChain seedChain_originalAlignment;
			seedChain_originalAlignment.reverse = false;
			int lastPrintedLevel = -1;
			assert(true_fa_sequence.size() == true_fa_coordinates_edgePath.size());
			for(unsigned int lI = 0; lI < true_fa_coordinates_edgePath.size(); lI++)
			{
				int thisLevel = true_fa_coordinates_edgePath.at(lI);
				if(thisLevel != -1)
				{
					if((lI != 0) && ((lastPrintedLevel + 1) != thisLevel))
					{
						while((lastPrintedLevel + 1) != thisLevel)
						{
							lastPrintedLevel++;
							assert(lastPrintedLevel < thisLevel);
							seedChain_originalAlignment.graph_aligned_levels.push_back(lastPrintedLevel);
							seedChain_originalAlignment.graph_aligned.push_back('_');
							seedChain_originalAlignment.sequence_aligned.push_back('_');
						}
					}
				}

				seedChain_originalAlignment.graph_aligned_levels.push_back(true_fa_coordinates_edgePath.at(lI));
				seedChain_originalAlignment.graph_aligned.push_back(true_fa_underlyingEdgelabels.at(lI));
				seedChain_originalAlignment.sequence_aligned.push_back(true_fa_sequence.at(lI));

				if(thisLevel != -1)
				{
					lastPrintedLevel = thisLevel;
				}
			}

			std::string alignedSequence_noGaps = Utilities::removeGaps(seedChain.sequence_aligned);
			assert(true_fa_sequence.size() <= seedChain_originalAlignment.sequence_aligned.size());

			std::string original_original_Sequence_noGaps = Utilities::removeGaps(true_fa_sequence);
			if(!(original_original_Sequence_noGaps == alignedSequence_noGaps))
			{
				std::cerr << "original_original_Sequence_noGaps" << ": " << original_original_Sequence_noGaps << "\n";
				std::cerr << "alignedSequence_noGaps" << ": " << alignedSequence_noGaps << "\n";
				std::cerr << "\n" << std::flush;
			}
			assert(original_original_Sequence_noGaps == alignedSequence_noGaps);



			std::string originalSequence_noGaps = Utilities::removeGaps(seedChain_originalAlignment.sequence_aligned);
			if(!(originalSequence_noGaps == alignedSequence_noGaps))
			{
				std::cerr << "true_fa_sequence" << ": " << true_fa_sequence << "\n";
				std::cerr << "original_original_Sequence_noGaps" << ": " << original_original_Sequence_noGaps << "\n";
				std::cerr << "originalSequence_noGaps" << ": " << originalSequence_noGaps << "\n";
				std::cerr << "alignedSequence_noGaps" << ": " << alignedSequence_noGaps << "\n";
				std::cerr << "\n" << std::flush;
			}
			assert(originalSequence_noGaps == alignedSequence_noGaps);

			seedChain.checkLevelContiguity();
			seedChain_originalAlignment.checkLevelContiguity();

			std::cout << "\t\t Original alignment:\n";
			std::cout << "\t\t\t" << Utilities::join(Utilities::ItoStr(seedChain_originalAlignment.graph_aligned_levels), ", ") << "\n";
			std::cout << "\t\t\t" << seedChain_originalAlignment.graph_aligned << "\n";
			std::cout << "\t\t\t" << seedChain_originalAlignment.sequence_aligned << "\n";

			double Score_original = eA->scoreOneAlignment(seedChain_originalAlignment, originalRead_alignmentOrientation);
			std::cout << "\t\t" << "Score" << ": " << Score_original << "\n";


			std::cout << std::flush;
		}
	};

	evaluateLevels(true);
	evaluateLevels(false);

	total_and_correct.first += forReturn.first;
	total_and_correct.second += forReturn.second;

	return forReturn;
}

std::pair<size_t, size_t> trueReadLevels::get_total_and_correct()
{
	return total_and_correct;
}

trueReadLevels::trueReadLevels(std::string R1_levels, std::string R2_levels) {

	total_and_correct.first = 0;
	total_and_correct.second = 0;

	auto read_true_levels = [&](std::string file, bool first) -> void {
		std::ifstream levelsStream;
		levelsStream.open(file.c_str());
		assert(levelsStream.is_open());

		std::string line;
		while(levelsStream.good())
		{
			std::getline(levelsStream, line);
			Utilities::eraseNL(line);
			if(line.length() == 0)
				continue;

			assert(line.substr(0, 1) == "@");

			std::string readID = line.substr(1);
			assert(levelsStream.good());

			std::string levels_str;
			std::getline(levelsStream, levels_str);
			Utilities::eraseNL(levels_str);

			std::vector<std::string> levels_vector = Utilities::split(levels_str, " ");


			std::vector<int> levels;
			levels.reserve ( levels_vector.size() );
			for(unsigned int i = 0; i < levels_vector.size(); i++)
			{
				std::string level = levels_vector.at(i);
				levels.push_back( Utilities::StrtoI(level) );
			}

			if(first)
			{
				trueLevels[readID].first = levels;
			}
			else
			{
				trueLevels[readID].second = levels;
			}


			std::string edgeLabels_str;
			std::getline(levelsStream, edgeLabels_str);
			Utilities::eraseNL(edgeLabels_str);
			if(first)
			{
				true_underlyingEdgeLabels[readID].first = edgeLabels_str;
			}
			else
			{
				true_underlyingEdgeLabels[readID].second = edgeLabels_str;
			}

			std::string fA_coordinates_edgePath_str;
			std::getline(levelsStream, fA_coordinates_edgePath_str);
			Utilities::eraseNL(fA_coordinates_edgePath_str);
			std::vector<std::string> fA_coordinates_edgePath_vector = Utilities::split(fA_coordinates_edgePath_str, " ");
			std::vector<int> fA_coordinates_edgePath;
			fA_coordinates_edgePath.reserve ( fA_coordinates_edgePath_vector.size() );
			for(unsigned int i = 0; i < fA_coordinates_edgePath_vector.size(); i++)
			{
				std::string level = fA_coordinates_edgePath_vector.at(i);
				fA_coordinates_edgePath.push_back( Utilities::StrtoI(level) );
			}

			if(first)
			{
				true_fullAlignment_coordinates_edgePath[readID].first = fA_coordinates_edgePath;
			}
			else
			{
				true_fullAlignment_coordinates_edgePath[readID].second = fA_coordinates_edgePath;
			}

			std::string fa_underlyingEdgeLabels_str;
			std::getline(levelsStream, fa_underlyingEdgeLabels_str);
			Utilities::eraseNL(fa_underlyingEdgeLabels_str);
			if(first)
			{
				true_fullAlignment_underlyingEdgeLabels[readID].first = fa_underlyingEdgeLabels_str;
			}
			else
			{
				true_fullAlignment_underlyingEdgeLabels[readID].second = fa_underlyingEdgeLabels_str;
			}

			std::string fa_sequence_str;
			std::getline(levelsStream, fa_sequence_str);
			Utilities::eraseNL(fa_sequence_str);
			if(first)
			{
				true_fullAlignment_sequence[readID].first = fa_sequence_str;
			}
			else
			{
				true_fullAlignment_sequence[readID].second = fa_sequence_str;
			}

			assert(fA_coordinates_edgePath.size() == fa_underlyingEdgeLabels_str.size());
			assert(fA_coordinates_edgePath.size() == fa_sequence_str.size());

		}

	};

	read_true_levels(R1_levels, true);
	read_true_levels(R2_levels, false);



}

trueReadLevels::~trueReadLevels() {
	// TODO Auto-generated destructor stub
}

} /* namespace simulator */
