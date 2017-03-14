/*
 * verboseSeedChain.h
 *
 *  Created on: 23.09.2015
 *      Author: AlexanderDilthey
 */

#ifndef MAPPER_READS_VERBOSESEEDCHAIN_H_
#define MAPPER_READS_VERBOSESEEDCHAIN_H_

#include "../../Graph/Edge.h"
#include "../../Utilities.h"

#include <string>
#include <vector>
#include <utility>
#include <assert.h>

namespace mapper {
namespace reads {

class verboseSeedChain {
public:
	verboseSeedChain();
	virtual ~verboseSeedChain();

	int sequence_begin;
	int sequence_end;
	bool reverse;
	int removed_columns_noGap_restriction;
	double improvement_through_bt;

	std::vector<bool> is_from_BWAseed;
	std::vector<Edge*> graph_aligned_edges;
	std::vector<int> graph_aligned_levels;
	std::string graph_aligned;
	std::string sequence_aligned;

	double mapQ;
	std::string mapQ_perPosition;

	void extendWithOtherSeedChain(const verboseSeedChain& otherChain, bool left);

	void print() const;
	void print(const std::vector<std::string>& levels) const;

	void checkChainConcordanceWithSequence(const std::string& sequence) const;

	void extendToFullSequenceLength(const std::string& sequence);

	std::set<std::string> getSegments(const std::vector<std::string>& levels) const;

	std::pair<int, int> quality()
	{
		assert(graph_aligned.length() == sequence_aligned.length());
		int L = graph_aligned.length();
		int M = 0;
		for(int i = 0; i < L; i++)
		{
			if(graph_aligned.at(i) == sequence_aligned.at(i))
			{
				M++;
			}
		}
		return make_pair(L, M);
	}

	std::pair<int, int> qualityLevels(const std::vector<int>& originalLevels, const std::string& originalSequence)
	{
		assert(originalLevels.size() == originalSequence.size());

		assert(sequence_begin >= 0);
		assert(sequence_begin <= sequence_end);
		unsigned int sequenceL = sequence_end - sequence_begin + 1;

		std::vector<int> originalLevels_definedSequence_thisAlignment = std::vector<int>(originalLevels.begin() + sequence_begin, originalLevels.begin() + sequence_begin + sequenceL);
		assert(originalLevels_definedSequence_thisAlignment.size() == sequenceL);

		std::vector<int> graph_aligned_levels_noSequenceGaps;
		{
			for(unsigned int i = 0; i < sequence_aligned.length(); i++)
			{
				if(sequence_aligned.at(i) != '_')
				{
					graph_aligned_levels_noSequenceGaps.push_back(graph_aligned_levels.at(i));
				}
			}
		}

		assert(graph_aligned_levels_noSequenceGaps.size() == sequenceL);

		int L = graph_aligned_levels_noSequenceGaps.size();
		int M = 0;
		for(int i = 0; i < L; i++)
		{
			if(graph_aligned_levels_noSequenceGaps.at(i) == (int)originalLevels_definedSequence_thisAlignment.at(i))
			{
				M++;
			}
		}

		return make_pair(L, M);
	}

	Node* alignment_firstNode() const
	{
		assert(graph_aligned_levels.size() == graph_aligned_edges.size());
		// int firstLevel = -1;
		for(unsigned int i = 0; i < graph_aligned_levels.size(); i++)
		{
			int level = graph_aligned_levels.at(i);
			if(level != -1)
			{
				return graph_aligned_edges.at(i)->From;
			}
		}
		return 0;
	}

	int alignment_firstLevel() const
	{
		int firstLevel = -1;
		for(unsigned int i = 0; i < graph_aligned_levels.size(); i++)
		{
			int level = graph_aligned_levels.at(i);
			if(level != -1)
			{
				firstLevel = level;
				break;
			}
		}
		return firstLevel;
	}

	std::vector<int> alignment_firstLevels(int nMaxLevels) const
	{
		assert(nMaxLevels > 0);
		std::vector<int> forReturn;
		forReturn.reserve(nMaxLevels);
		for(unsigned int i = 0; i < graph_aligned_levels.size(); i++)
		{
			int level = graph_aligned_levels.at(i);
			if(level != -1)
			{
				forReturn.push_back(level);
				if((int)forReturn.size() >= nMaxLevels)
				{
					break;
				}
			}
		}
		return forReturn;
	}

	int alignment_lastLevel() const
	{
		int lastLevel = -1;
		// std::cout << "graph_aligned_levels.size(): " << graph_aligned_levels.size()  << "\n" << std::flush;

		if(graph_aligned_levels.size() > 0)
		{
			for(int i = ((int)graph_aligned_levels.size() - 1); i >= 0; i--)
			{
				// std::cout << "i: " << i << "/" << graph_aligned_levels.size() << "\n" << std::flush;

				assert(i >= 0);
				assert(i < (int)graph_aligned_levels.size());
				int level = graph_aligned_levels.at(i);
				if(level != -1)
				{
					lastLevel = level;
					break;
				}
			}
		}
		else
		{
			// std::cerr << "Warning: call alignment_lastLevel() on unmapped (?) read!\n" << std::flush;
		}
		return lastLevel;
	}

	Node* alignment_lastNode() const
	{
		assert(graph_aligned_edges.size() == graph_aligned_levels.size());
		if(graph_aligned_levels.size() > 0)
		{
			for(int i = ((int)graph_aligned_levels.size() - 1); i >= 0; i--)
			{
				// std::cout << "i: " << i << "/" << graph_aligned_levels.size() << "\n" << std::flush;

				assert(i >= 0);
				assert(i < (int)graph_aligned_levels.size());
				int level = graph_aligned_levels.at(i);
				if(level != -1)
				{
					return graph_aligned_edges.at(i)->To;
				}
			}
		}
		return 0;
	}

	std::vector<int> alignment_lastLevels(int nMaxLevels) const
	{
		std::vector<int> forReturn;
		forReturn.reserve(nMaxLevels);
		if(graph_aligned_levels.size() > 0)
		{
			for(int i = ((int)graph_aligned_levels.size() - 1); i >= 0; i--)
			{
				assert(i >= 0);
				assert(i < (int)graph_aligned_levels.size());
				int level = graph_aligned_levels.at(i);
				if(level != -1)
				{
					forReturn.push_back(level);
					if((int)forReturn.size() >= nMaxLevels)
					{
						break;
					}
				}
			}
		}
		return forReturn;
	}

	std::map<int, int> alignment_end_originalSequenceAnchors(int nMaxLevels, const std::vector<std::map<int, int>>& graphLevel_2_underlyingSequencePositions) const
	{
		std::vector<int> levels_for_anchors = alignment_lastLevels(nMaxLevels);
		std::map<int, int> forReturn;
		for(unsigned int levelI = 0; levelI < levels_for_anchors.size(); levelI++)
		{
			int level = levels_for_anchors.at(levelI);
			assert(level >= 0);
			if(graphLevel_2_underlyingSequencePositions.at(level).size())
			{
				for(std::map<int, int>::const_iterator originalSequenceIt = graphLevel_2_underlyingSequencePositions.at(level).begin(); originalSequenceIt != graphLevel_2_underlyingSequencePositions.at(level).end(); originalSequenceIt++)
				{
					int originalSequenceID = originalSequenceIt->first;
					int positionAlongSequenceID = originalSequenceIt->second;
					assert(positionAlongSequenceID >= 0);
					if(forReturn.count(originalSequenceID) == 0)
					{
						forReturn[originalSequenceID] = positionAlongSequenceID;
					}
				}
			}
		}

		return forReturn;
	}

	std::map<int, int> alignment_begin_originalSequenceAnchors(int nMaxLevels, const std::vector<std::map<int, int>>& graphLevel_2_underlyingSequencePositions) const
	{
		std::vector<int> levels_for_anchors = alignment_firstLevels(nMaxLevels);
		std::map<int, int> forReturn;
		for(unsigned int levelI = 0; levelI < levels_for_anchors.size(); levelI++)
		{
			int level = levels_for_anchors.at(levelI);
			assert(level >= 0);
			if(graphLevel_2_underlyingSequencePositions.at(level).size())
			{
				for(std::map<int, int>::const_iterator originalSequenceIt = graphLevel_2_underlyingSequencePositions.at(level).begin(); originalSequenceIt != graphLevel_2_underlyingSequencePositions.at(level).end(); originalSequenceIt++)
				{
					int originalSequenceID = originalSequenceIt->first;
					int positionAlongSequenceID = originalSequenceIt->second;
					assert(positionAlongSequenceID >= 0);
					if(forReturn.count(originalSequenceID) == 0)
					{
						forReturn[originalSequenceID] = positionAlongSequenceID;
					}
				}
			}
		}

		return forReturn;
	}

	void checkLevelContiguity() const
	{
		assert(graph_aligned.size() == graph_aligned_levels.size());
		assert(sequence_aligned.size() == graph_aligned_levels.size());

		if(graph_aligned_levels.size() > 1)
		{
			bool error = false;
			int last_examined_level = -1;
			for(unsigned int lI = 0; lI < graph_aligned_levels.size(); lI++)
			{
				int thisLevel = graph_aligned_levels.at(lI);

				if(thisLevel != -1)
				{
					if(last_examined_level != -1)
					{
						if(!((last_examined_level + 1) == thisLevel))
						{
							error = true;
						}
					}
					last_examined_level = thisLevel;
				}
			}

			if(error)
			{
				std::cerr << "Chain level contiguity error!\n";
				std::cerr << Utilities::join(Utilities::ItoStr(graph_aligned_levels), ",") << "\n" << std::flush;
			}
			assert(error == false);
		}
	}
};

class verboseSeedChainPair
{
public:
	std::string readID;
	std::pair<verboseSeedChain, verboseSeedChain> chains;
	double mapQ;

	verboseSeedChainPair()
	{
		mapQ = -1;
	}

	void print() const
	{
		print({});
	}
	
	void print(const std::vector<std::string>& levels) const
	{
		std::cout << "verboseSeedChainPair:" << "\n";
		std::cout << "mapQ:" << mapQ << "\n";
		std::cout << "Chain I\n================================\n";
		chains.first.print(levels);
		std::cout << "Chain II\n================================\n";
		chains.second.print(levels);
		std::cout << "\n";		
	}
};

} /* namespace reads */
} /* namespace mapper */

#endif /* MAPPER_READS_VERBOSESEEDCHAIN_H_ */
