/*
 * BAMalignment.cpp
 *
 *  Created on: 22.09.2015
 *      Author: AlexanderDilthey
 */

#include "PRGContigBAMAlignment.h"

#include <assert.h>
#include <iostream>
#include "../../Utilities.h"

namespace mapper {
namespace reads {

PRGContigBAMAlignment::PRGContigBAMAlignment() {
	// TODO Auto-generated constructor stub

}

PRGContigBAMAlignment::~PRGContigBAMAlignment() {
	// TODO Auto-generated destructor stub
}

void PRGContigBAMAlignment::removeSequenceCharacters(bool left, int n)
{
	if(left)
	{
		int nonGap_characters = 0;
		int newStartCoordinate_left = -1;
		for(unsigned int i = 0; i < sequence_aligned.length(); i++)
		{
			if(sequence_aligned.at(i) != '_')
			{
				nonGap_characters++;
				if(nonGap_characters == n)
				{
					newStartCoordinate_left = i + 1;
					break;
				}
			}
		}
		assert(newStartCoordinate_left != -1);
		assert(newStartCoordinate_left < (int)(sequence_aligned.length()));

		graph_aligned_levels = std::vector<int>(graph_aligned_levels.begin() + newStartCoordinate_left, graph_aligned_levels.end());
		graph_aligned = std::string(graph_aligned.begin() + newStartCoordinate_left, graph_aligned.end());
		sequence_aligned = std::string(sequence_aligned.begin() + newStartCoordinate_left, sequence_aligned.end());
		sequence_aligned_startInRaw += n;
	}
	else
	{
		int nonGap_characters = 0;
		int newStartCoordinate_right = -1;
		for(int i = sequence_aligned.length() - 1; i >= 0; i--)
		{
			if(sequence_aligned.at(i) != '_')
			{
				nonGap_characters++;
				if(nonGap_characters == n)
				{
					newStartCoordinate_right = i - 1;
					break;
				}
			}
		}
		assert(newStartCoordinate_right != -1);
		assert(newStartCoordinate_right >= 0);

		int remainingCharacters = newStartCoordinate_right + 1;
		graph_aligned_levels = std::vector<int>(graph_aligned_levels.begin(), graph_aligned_levels.begin() + remainingCharacters);
		graph_aligned = std::string(graph_aligned.begin(), graph_aligned.begin() + remainingCharacters);
		sequence_aligned = std::string(sequence_aligned.begin(), sequence_aligned.begin() + remainingCharacters);
		sequence_aligned_stopInRaw -= n;
	}
}

void PRGContigBAMAlignment::print()
{
	std::cout << "PRGContigBAMAlignment:\n";
	std::cout << "graph_aligned_levels       " << ": " << Utilities::join(Utilities::ItoStr(graph_aligned_levels), ", ") << "\n";
	std::cout << "graph_aligned              " << ": " << graph_aligned << "\n";
	std::cout << "sequence_aligned           " << ": " << sequence_aligned << "\n";
	std::cout << "sequence_aligned_startInRaw" << ": " << sequence_aligned_startInRaw << "\n";
	std::cout << "sequence_aligned_stopInRaw " << ": " << sequence_aligned_stopInRaw << "\n";
	std::cout << "reverse                    " << ": " << reverse << "\n";
	std::cout << "\n" << std::flush;
}

void PRGContigBAMAlignment::checkAlignmentConcordanceWithSequence(const std::string& sequence) const
{
	assert(sequence_aligned_startInRaw <= sequence_aligned_stopInRaw);
	assert(sequence_aligned_startInRaw >= 0);
	assert(sequence_aligned_stopInRaw < (int)sequence.length());

	assert(sequence.find('_') == std::string::npos);
	std::string subSequence = sequence.substr(sequence_aligned_startInRaw, sequence_aligned_stopInRaw - sequence_aligned_startInRaw + 1);
//	if(reverse)
//	{
//		subSequence = Utilities::seq_reverse_complement(subSequence);
//	}

	std::string sequence_aligned_noGaps = filter(sequence_aligned, [](unsigned char c){return (c != '_');});

	if(!(sequence_aligned_noGaps == subSequence))
	{
		std::cerr << "Mismatch in PRGContigBAMAlignment::checkAlignmentConcordanceWithSequence!\n";
		std::cerr << "\t" << "sequence" << ": " << sequence << "\n";
		std::cerr << "\t" << "sequence_aligned" << ": " << sequence_aligned << "\n";
		std::cerr << "\t" << "sequence_aligned_noGaps" << ": " << sequence_aligned_noGaps << "\n";
		std::cerr << "\t" << "subSequence" << ": " << subSequence << "\n";
		std::cerr << "\t" << "sequence_aligned_startInRaw" << ": " << sequence_aligned_startInRaw << "\n";
		std::cerr << "\t" << "sequence_aligned_stopInRaw" << ": " << sequence_aligned_stopInRaw << "\n";
		std::cerr << "\n" << std::flush;
	}
	
	assert(sequence_aligned_noGaps == subSequence);
}

} /* namespace reads */
} /* namespace mapper */
