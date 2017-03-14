/*
 * verboseSeedChain.cpp
 *
 *  Created on: 23.09.2015
 *      Author: AlexanderDilthey
 */

#include "verboseSeedChain.h"

namespace mapper {
namespace reads {

verboseSeedChain::verboseSeedChain() {
	mapQ = -1;
	removed_columns_noGap_restriction = -1;
	improvement_through_bt = -1;
}

verboseSeedChain::~verboseSeedChain() {

}

void verboseSeedChain::extendWithOtherSeedChain(const verboseSeedChain& otherChain, bool left)
{
	assert(otherChain.reverse == reverse);
	if(left)
	{
		assert((otherChain.sequence_end + 1) == sequence_begin);
		sequence_begin = otherChain.sequence_begin;
		graph_aligned_levels.insert(graph_aligned_levels.begin(), otherChain.graph_aligned_levels.begin(), otherChain.graph_aligned_levels.end());
		graph_aligned_edges.insert(graph_aligned_edges.begin(), otherChain.graph_aligned_edges.begin(), otherChain.graph_aligned_edges.end());
		graph_aligned.insert(graph_aligned.begin(), otherChain.graph_aligned.begin(), otherChain.graph_aligned.end());
		sequence_aligned.insert(sequence_aligned.begin(), otherChain.sequence_aligned.begin(), otherChain.sequence_aligned.end());
		is_from_BWAseed.insert(is_from_BWAseed.begin(), otherChain.is_from_BWAseed.begin(), otherChain.is_from_BWAseed.end());
	}
	else
	{
		assert(otherChain.sequence_begin == (sequence_end + 1));
		sequence_end = otherChain.sequence_end;
		graph_aligned_levels.insert(graph_aligned_levels.end(), otherChain.graph_aligned_levels.begin(), otherChain.graph_aligned_levels.end());
		graph_aligned_edges.insert(graph_aligned_edges.end(), otherChain.graph_aligned_edges.begin(), otherChain.graph_aligned_edges.end());
		graph_aligned.insert(graph_aligned.end(), otherChain.graph_aligned.begin(), otherChain.graph_aligned.end());
		sequence_aligned.insert(sequence_aligned.end(), otherChain.sequence_aligned.begin(), otherChain.sequence_aligned.end());
		is_from_BWAseed.insert(is_from_BWAseed.end(), otherChain.is_from_BWAseed.begin(), otherChain.is_from_BWAseed.end());
	}
}

void verboseSeedChain::checkChainConcordanceWithSequence(const std::string& sequence) const
{
	assert(graph_aligned_edges.size() == graph_aligned_levels.size());
	assert(graph_aligned_edges.size() == graph_aligned.size());
	assert(graph_aligned_edges.size() == sequence_aligned.size());

	assert(sequence_begin <= sequence_end);
	assert(sequence_begin >= 0);
	assert(sequence_end < (int)sequence.length());

	std::string seedChain_sequence_aligned_noGaps = filter(sequence_aligned, [](unsigned char c){return (c != '_');});
	assert(seedChain_sequence_aligned_noGaps.find('_') == std::string::npos);
	assert(sequence.find('_') == std::string::npos);

	std::string sequence_coveredByChain = sequence.substr(sequence_begin, sequence_end - sequence_begin + 1);
	if(reverse)
	{
//		sequence_coveredByChain = Utilities::seq_reverse_complement(sequence_coveredByChain);
	}

	if(!(sequence_coveredByChain == seedChain_sequence_aligned_noGaps))
	{
		std::cerr << "Chain concordance error!\n" << std::flush;
		std::cerr << "\t"<< "sequence_coveredByChain" << ": " << sequence_coveredByChain << "\n";
		std::cerr << "\t"<< "seedChain_sequence_aligned_noGaps" << ": " << seedChain_sequence_aligned_noGaps << "\n";
		print();
	}

	assert(sequence_coveredByChain == seedChain_sequence_aligned_noGaps);
}

void verboseSeedChain::extendToFullSequenceLength(const std::string& sequence)
{
	int missing_characters_left = sequence_begin;
	int missing_characters_right = sequence.length() - sequence_end - 1;

	if(missing_characters_left)
	{
		std::string insertSequence = sequence.substr(0, missing_characters_left);

		std::vector<int> insertGraphLevels;
		insertGraphLevels.resize(missing_characters_left, -1);

		std::vector<Edge*> insertAlignedEdges;
		insertAlignedEdges.resize(missing_characters_left, 0);

		std::string insertGraphAligned;
		insertGraphAligned.resize(missing_characters_left, '_');

		std::vector<bool> insert_is_from_BWAseed;
		insert_is_from_BWAseed.resize(missing_characters_left, 0);

		graph_aligned_edges.insert(graph_aligned_edges.begin(), insertAlignedEdges.begin(), insertAlignedEdges.end());
		graph_aligned_levels.insert(graph_aligned_levels.begin(), insertGraphLevels.begin(), insertGraphLevels.end());
		graph_aligned.insert(graph_aligned.begin(), insertGraphAligned.begin(), insertGraphAligned.end());
		sequence_aligned.insert(sequence_aligned.begin(), insertSequence.begin(), insertSequence.end());
		is_from_BWAseed.insert(is_from_BWAseed.begin(), insert_is_from_BWAseed.begin(), insert_is_from_BWAseed.end());

		sequence_begin = 0;
	}

	if(missing_characters_right)
	{
		std::string insertSequence = sequence.substr(sequence.length() - missing_characters_right, missing_characters_right);
		assert((int)insertSequence.length() == missing_characters_right);

		std::vector<int> insertGraphLevels;
		insertGraphLevels.resize(missing_characters_right, -1);

		std::vector<Edge*> insertAlignedEdges;
		insertAlignedEdges.resize(missing_characters_right, 0);

		std::string insertGraphAligned;
		insertGraphAligned.resize(missing_characters_right, '_');

		std::vector<bool> insert_is_from_BWAseed;
		insert_is_from_BWAseed.resize(missing_characters_right, 0);

		graph_aligned_edges.insert(graph_aligned_edges.end(), insertAlignedEdges.begin(), insertAlignedEdges.end());
		graph_aligned_levels.insert(graph_aligned_levels.end(), insertGraphLevels.begin(), insertGraphLevels.end());
		graph_aligned.insert(graph_aligned.end(), insertGraphAligned.begin(), insertGraphAligned.end());
		sequence_aligned.insert(sequence_aligned.end(), insertSequence.begin(), insertSequence.end());
		is_from_BWAseed.insert(is_from_BWAseed.end(), insert_is_from_BWAseed.begin(), insert_is_from_BWAseed.end());

		sequence_end = sequence.length() - 1;
	}

	checkChainConcordanceWithSequence(sequence);
}

std::set<std::string> verboseSeedChain::getSegments(const std::vector<std::string>& levels) const
{
	std::set<std::string> segments;
	if(levels.size())
	{
		for(int level : graph_aligned_levels)
		{
			if(level != -1)
			{
				assert(level >= 0);
				assert(level < (int)levels.size());
				std::string levelID = levels.at(level);
				std::vector<std::string> parts = Utilities::split(levelID, "_");
				assert(parts.size() >= 3);
				std::vector<std::string> parts_for_print;
				for(unsigned int i = 1; i < (parts.size() - 1); i++)
				{
					parts_for_print.push_back(parts.at(i));
				}
				std::string id_for_print = Utilities::join(parts_for_print, "_");
				segments.insert(id_for_print);
			}
		}
	}

	return segments;
}
void verboseSeedChain::print(const std::vector<std::string>& levels) const
{
	std::set<std::string> segments = getSegments(levels);
	if(levels.size())
	{
		for(int level : graph_aligned_levels)
		{
			if(level != -1)
			{
				assert(level >= 0);
				assert(level < (int)levels.size());
				std::string levelID = levels.at(level);
				std::vector<std::string> parts = Utilities::split(levelID, "_");
				assert(parts.size() >= 3);
				std::vector<std::string> parts_for_print;
				for(unsigned int i = 1; i < (parts.size() - 1); i++)
				{
					parts_for_print.push_back(parts.at(i));
				}
				std::string id_for_print = Utilities::join(parts_for_print, "_");
				segments.insert(id_for_print);
			}
		}
	}
	std::vector<std::string> segments_vector(segments.begin(), segments.end());
	
	std::cout << "verboseSeedChain:\n";
	std::cout << "mapQ:" << mapQ << "\n";
	std::cout << "Segments:" << Utilities::join(segments_vector, ", ") << "\n";
	std::cout << "graph_aligned_levels       " << ": " << Utilities::join(Utilities::ItoStr(graph_aligned_levels), ", ") << "\n";
	std::cout << "graph_aligned              " << ": " << graph_aligned << "\n";
	std::cout << "sequence_aligned           " << ": " << sequence_aligned << "\n";
	std::cout << "sequence_begin" << ": " << sequence_begin << "\n";
	std::cout << "sequence_end " << ": " << sequence_end << "\n";
	std::cout << "reverse                    " << ": " << reverse << "\n";
	std::cout << "\n" << std::flush;	
}
void verboseSeedChain::print() const
{
	print({});
}

} /* namespace reads */
} /* namespace mapper */
