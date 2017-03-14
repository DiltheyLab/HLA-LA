/*
 * extensionAligner.cpp
 *
 *  Created on: 26.09.2015
 *      Author: AlexanderDilthey
 */

#include "extensionAligner.h"

#include "alignmentContext.h"

#include "../reads/oneReadPair.h"
#include <algorithm>
#include <assert.h>
#include <omp.h>

namespace mapper {
namespace aligner {

extensionAligner::extensionAligner(Graph* g_) : alignerBase(g_) {
	threads_initialized = 0;
	paranoid = true;

	if(!g->get_haveComputedGapEdgePaths())
	{
		g->computeGapEdgePaths();
	}

	init_for_threads(1);
}

extensionAligner::~extensionAligner() {
	// TODO Auto-generated destructor stub
}

void extensionAligner::init_for_threads(unsigned int threads)
{
	if(threads_initialized != threads)
	{
		rng_seeds.resize(threads);
		srand(time(NULL));
		for(unsigned int tI = 0; tI < threads; tI++)
		{
			rng_seeds.at(tI) = rand();
		}

		threads_initialized = threads;
	}
}

 
double extensionAligner::scoreOneAlignment(const mapper::reads::verboseSeedChain& alignment, const mapper::reads::oneRead& underlyingRead) const
{
	int indexIntoOriginalReadData = -1;

	bool conservativeReadQualities = true;

	double rate_deletions = log(0.001);
	double rate_insertions = log(0.001);
	double rate_match_mismatch = log(1 - exp(rate_deletions) - exp(rate_insertions));
	double combined_log_likelihood = 0;

	int totalMismatches = 0;

	bool verbose = false;

	if(verbose) std::cout << "SCORE\n";

	for(unsigned int cI = 0; cI < alignment.sequence_aligned.length(); cI++)
	{
		std::string sequenceCharacter = alignment.sequence_aligned.substr(cI, 1);
		std::string graphCharacter = alignment.graph_aligned.substr(cI, 1);
		int graphLevel = alignment.graph_aligned_levels.at(cI);

		if(sequenceCharacter != "_")
		{
			indexIntoOriginalReadData++;
			int indexIntoOriginalReadData_correctlyAligned = indexIntoOriginalReadData;
			if(alignment.reverse)
			{
				indexIntoOriginalReadData_correctlyAligned = underlyingRead.sequence.length() - indexIntoOriginalReadData_correctlyAligned - 1;
			}
			assert(indexIntoOriginalReadData_correctlyAligned >= 0);
			assert(indexIntoOriginalReadData_correctlyAligned < (int)underlyingRead.sequence.length());;

			if(underlyingRead.name != "SOFTCLIPPED")
			{
				std::string underlyingReadCharacter = underlyingRead.sequence.substr(indexIntoOriginalReadData_correctlyAligned, 1);
				if(alignment.reverse)
				{
					underlyingReadCharacter = Utilities::seq_reverse_complement(underlyingReadCharacter);
				}
				if(!(underlyingReadCharacter == sequenceCharacter))
				{
					std::cerr << "!(underlyingReadCharacter == sequenceCharacter)" << "\n";
					std::cerr << "underlyingReadCharacter" << ": " << underlyingReadCharacter << "\n";
					std::cerr << "sequenceCharacter" << ": " << sequenceCharacter << "\n";
					std::cerr << "cI" << ": " << cI << "/" << alignment.sequence_aligned.length() << "\n";
					std::cerr << "alignment.reverse" << ": " << alignment.reverse << "\n";
					std::cerr << "alignment.sequence_aligned" << ":\t" << alignment.sequence_aligned << "\n";
					std::cerr << "alignment.graph_aligned" << ":\t" << alignment.graph_aligned << "\n";
					std::cerr << "underlyingRead.sequence" << ":\t" << underlyingRead.sequence << "\n";
					std::cerr << "underlyingRead.name" << ":\t" << underlyingRead.name << "\n";

					std::cerr << std::flush;
				}
				assert(underlyingReadCharacter == sequenceCharacter);
			}

			if(graphCharacter == "_")
			{
				// sequence non gap, graph gap -- insertion
				combined_log_likelihood += (rate_insertions + log(1.0/4.0));
				totalMismatches++;
			}
			else
			{
				// two well-defined characters
				combined_log_likelihood += rate_match_mismatch;
				char qualityCharacter = underlyingRead.quality.at(indexIntoOriginalReadData_correctlyAligned);
				double pCorrect = Utilities::PhredToPCorrect(qualityCharacter);
				if(conservativeReadQualities && (pCorrect > 0.999))
				{
					pCorrect = 0.999;
				}
				if(pCorrect == 0)
				{
					pCorrect = 0.00001;
				}
				assert((pCorrect > 0) && (pCorrect <= 1));
				if(sequenceCharacter == graphCharacter)
				{
					combined_log_likelihood += log(pCorrect);
				}
				else
				{
					double pIncorrect = 1 - pCorrect;
					pIncorrect *= (1.0/3.0);
					assert((pIncorrect > 0) && (pIncorrect < 1));
					combined_log_likelihood += log(pIncorrect);
					totalMismatches++;
				}
			}

		}
		else
		{
			assert(sequenceCharacter == "_");
			if(graphCharacter == "_")
			{
				// sequence gap, graph gap - likelihood 1
				assert(graphLevel != -1);
			}
			else
			{
				// sequence gap, graph non gap - deletion in sequence
				combined_log_likelihood += rate_deletions;
				totalMismatches++;
			}
		}

		if(verbose) std::cout << "\t" << cI << "\t" << combined_log_likelihood << "\n";
	}

	double l = exp(combined_log_likelihood);
	if(0 && !((l > 0) && (l <= 1)))
	{
		std::cerr << "! ((l > 0) && (l <= 1))" << "\n";
		std::cerr << "l" << ": " << l << "\n";  
		std::cerr << "combined_log_likelihood" << ": " << combined_log_likelihood << "\n";
		std::cerr << "\n" << std::flush;
		
	}
	assert((l >= 0) && (l <= 1));
	return combined_log_likelihood;
}



mapper::reads::verboseSeedChain extensionAligner::extendSeedChain(const std::string& sequence, const mapper::reads::verboseSeedChain& seedChain) const
{
	assert(seedChain.sequence_begin <= seedChain.sequence_end);
	assert(seedChain.sequence_begin >= 0);
	if(!(seedChain.sequence_end < (int)sequence.length()))
	{
		std::cerr << "seedChain.sequence_end" << ": " << seedChain.sequence_end << "\n";
		std::cerr << "sequence.length()" << ": " << sequence.length() << "\n";
		std::cerr << std::flush;
	}
	assert(seedChain.sequence_end < (int)sequence.length());

	if(paranoid)
	{
		seedChain.checkChainConcordanceWithSequence(sequence);
	}

	mapper::reads::verboseSeedChain forReturn = seedChain;

	// left extension
	assert(seedChain.graph_aligned_edges.size());

	VirtualNWTable_Unique* vNW = 0;
	alignmentContext* aC = 0;
	auto prepare_vNW = [&]() -> void {
		if(vNW == 0)
		{
			aC = new alignmentContext();
			aC->sequence = &sequence;
			aC->aligner = this;
			vNW = new VirtualNWTable_Unique(aC);
		}
	};

	if(seedChain.sequence_begin != 0)
	{
		Edge* existingChain_firstEdge = seedChain.graph_aligned_edges.front();
		Node* firstNode = existingChain_firstEdge->From;
		if(firstNode->level > 0)
		{
			prepare_vNW();

			std::vector<mapper::reads::verboseSeedChain> extensionChains = fullNeedleman_diagonal_extension_gapJumper(
//			std::vector<mapper::reads::verboseSeedChain> extensionChains = fullNeedleman_diagonal_extension(
					sequence,
					seedChain.sequence_begin,
					firstNode->level,
					nodesPerLevel_ordered_rev.at(firstNode->level).at(firstNode),
					0,
					0,
					-16,
					vNW,
					false,
					false,
					true
			);

			for(unsigned int i = 0; i < extensionChains.size(); i++)
			{
				extensionChains.at(i).reverse = seedChain.reverse;
				//std::cout << "Left extension #" << i << ":\n" << std::flush;
				//extensionChain.at(i).print();
				extensionChains.at(i).checkChainConcordanceWithSequence(sequence);
				extensionChains.at(i).is_from_BWAseed.resize(extensionChains.at(i).graph_aligned_levels.size(), 0);
			}

			std::sort(extensionChains.begin(), extensionChains.end(), [](const mapper::reads::verboseSeedChain& a, const mapper::reads::verboseSeedChain& b){
				return (a.graph_aligned_levels.size() < b.graph_aligned_levels.size());
			});

			if(extensionChains.size() > 0)
			{
				if(extensionChains.size() > 1)
				{
					assert(extensionChains.back().graph_aligned_levels.size() >= extensionChains.at(extensionChains.size() - 2).graph_aligned_levels.size());
				}

				assert(extensionChains.at(0).is_from_BWAseed.size() == extensionChains.at(0).graph_aligned_levels.size());
				forReturn.extendWithOtherSeedChain(extensionChains.at(0), true);
				assert(forReturn.is_from_BWAseed.size() == forReturn.graph_aligned_levels.size());
			}
		}
	}

	// right extension
	if(seedChain.sequence_end != (int)(sequence.length() - 1))
	{
		Edge* existingChain_lastEdge = seedChain.graph_aligned_edges.back();
		Node* lastNode = existingChain_lastEdge->To;
		if(lastNode->level < (g->NodesPerLevel.size()-1))
		{
			prepare_vNW();

			std::vector<mapper::reads::verboseSeedChain> extensionChains = fullNeedleman_diagonal_extension_gapJumper(
//			std::vector<mapper::reads::verboseSeedChain> extensionChains = fullNeedleman_diagonal_extension(
					sequence,
					seedChain.sequence_end + 1,
					lastNode->level,
					nodesPerLevel_ordered_rev.at(lastNode->level).at(lastNode),
					g->NodesPerLevel.size() - 1,
					sequence.length(),
					-16,
					vNW,
					true,
					false,
					true
			);

			for(unsigned int i = 0; i < extensionChains.size(); i++)
			{
				extensionChains.at(i).reverse = seedChain.reverse;
				//std::cout << "Right extension #" << i << ":\n" << std::flush;
				//extensionChain.at(i).print();
				extensionChains.at(i).checkChainConcordanceWithSequence(sequence);
				extensionChains.at(i).is_from_BWAseed.resize(extensionChains.at(i).graph_aligned_levels.size(), 0);
			}

			std::sort(extensionChains.begin(), extensionChains.end(), [](const mapper::reads::verboseSeedChain& a, const mapper::reads::verboseSeedChain& b){
				return (a.graph_aligned_levels.size() < b.graph_aligned_levels.size());
			});

			if(extensionChains.size() > 0)
			{
				if(extensionChains.size() > 1)
				{
					assert(extensionChains.back().graph_aligned_levels.size() >= extensionChains.at(extensionChains.size() - 2).graph_aligned_levels.size());
				}

				assert(extensionChains.at(0).is_from_BWAseed.size() == extensionChains.at(0).graph_aligned_levels.size());
				forReturn.extendWithOtherSeedChain(extensionChains.back(), false);
				assert(forReturn.is_from_BWAseed.size() == forReturn.graph_aligned_levels.size());
			}
		}
	}

	forReturn.checkChainConcordanceWithSequence(sequence);

	forReturn.extendToFullSequenceLength(sequence);

	if(vNW != 0)
	{
		assert(aC != 0);
		delete(vNW);
		delete(aC);
	}

	return forReturn;
}

std::vector<mapper::reads::verboseSeedChain> extensionAligner::fullNeedleman_diagonal_extension_gapJumper(const std::string& sequence, int start_sequence, int startLevel_graph, int startZ_graph, int maxLevel_graph, int maxPosition_sequence, int diagonal_stop_threshold, VirtualNWTable_Unique* blockedPathsTable, bool directionPositive, bool returnGlobalScore, bool preferSequenceCompleAlignments) const
{
	assert(omp_get_num_threads() == (int)threads_initialized);

	assert(! ((returnGlobalScore) && (preferSequenceCompleAlignments)));

	assert(g != 0);
	if(directionPositive)
	{
		if(!((maxLevel_graph == -1) || (maxLevel_graph > startLevel_graph)))
		{
			std::cerr << "maxLevel_graph: " << maxLevel_graph << "\n";
			std::cerr << "startLevel_graph: " << startLevel_graph << "\n" << std::flush;
		}
		if(!((maxPosition_sequence == -1) || (maxPosition_sequence > start_sequence)))
		{
			std::cerr << "maxPosition_sequence: " << maxPosition_sequence << "\n";
			std::cerr << "start_sequence: " << start_sequence << "\n" << std::flush;
		}
		assert((maxLevel_graph == -1) || (maxLevel_graph > startLevel_graph));
		assert((maxPosition_sequence == -1) || (maxPosition_sequence > start_sequence));
	}
	else
	{
		assert((maxLevel_graph == -1) || (maxLevel_graph < startLevel_graph));
		assert((maxPosition_sequence == -1) || (maxPosition_sequence < start_sequence));
	}

	double minusInfinity = -1 * numeric_limits<double>::max();

	if(returnGlobalScore)
	{
		diagonal_stop_threshold = minusInfinity;
		assert(maxLevel_graph == -1);
	}

	class mScore {
	public:
		double D;
		double GraphGap;
		double SequenceGap;
	};

	class mScore_backtrace {
	public:
		backtraceStep_affine D;
		backtraceStep_affine GraphGap;
		backtraceStep_affine SequenceGap;
	};

	class mScore_alternatives {
	public:
		std::vector<double> D;
		std::vector<double> GraphGap;
		std::vector<double> SequenceGap;
	};

	class mScore_backtrace_alternatives {
	public:
		std::vector<backtraceStep_affine> D;
		std::vector<backtraceStep_affine> GraphGap;
		std::vector<backtraceStep_affine> SequenceGap;
	};

	std::map<int, std::map<int, std::map<int, mScore>> > scores;
	std::map<int, std::map<int, std::map<int, mScore_backtrace>> > scores_backtrace;

	std::vector<std::vector<int> > coordinates_for_backtracking;
	std::vector<std::vector<int> > m1_diagonal;
	std::vector<std::vector<int> > m2_diagonal;

	bool verbose = false;

	if(verbose)
	{
		std::cout << "fullNeedleman_affine_diagonal_extension(..) called, direction " << directionPositive << ".\n" << std::flush;
		std::cout << "\t" << "start_sequence" << ": " << start_sequence << "\n";
		std::cout << "\t" << "maxPosition_sequence" << ": " << maxPosition_sequence << "\n";
		std::cout << "\t" << "startLevel_graph" << ": " << startLevel_graph << "\n";
		std::cout << "\t" << "maxLevel_graph" << ": " << maxLevel_graph << "\n";
		std::cout << "\t" << "directionPositive" << ": " << directionPositive << "\n";
		std::cout << "\t" << "returnGlobalScore" << ": " << returnGlobalScore << "\n";
		std::cout << "\t" << "preferSequenceCompleAlignments" << ": " << preferSequenceCompleAlignments << "\n";

		if(directionPositive)
		{
			std::string coveredSequence = sequence.substr(start_sequence, maxPosition_sequence - start_sequence);
			std::cout << "\tSequeence: " << coveredSequence << "\n" << std::flush;
		}
		else
		{
			std::string coveredSequence = sequence.substr(maxPosition_sequence, start_sequence -  maxPosition_sequence);
			std::cout << "\tSequeence: " << coveredSequence << "\n" << std::flush;

		}

		std::cout << std::flush;
	}

	unsigned int levels = g->NodesPerLevel.size();
	unsigned int sequenceLength = sequence.length();
	int diagonals = sequenceLength + levels - 1;
	int max_levelI = levels - 1;
	int max_seqI = sequenceLength;
	int min_levelI = 0;
	int min_seqI = 0;

	if(maxLevel_graph != -1)
	{
		if(directionPositive)
		{
			max_levelI = maxLevel_graph;
		}
		else
		{
			min_levelI = maxLevel_graph;
		}
	}
	if(maxPosition_sequence != -1)
	{
		if(directionPositive)
		{
			max_seqI = maxPosition_sequence;
		}
		else
		{
			min_seqI = maxPosition_sequence;
		}
	}

	assert(startLevel_graph >= min_levelI);
	assert(start_sequence >= min_seqI);
	assert(startLevel_graph <= max_levelI);
	assert(start_sequence <= max_seqI);


	assert(min_levelI >= 0);
	assert(max_levelI <= (int)(levels - 1));
	assert(max_levelI > min_levelI);

	assert(min_seqI >= 0);
	assert(max_seqI <= (int)sequenceLength);
	assert(max_seqI > min_seqI);

	double currentMaximum = 0;
	std::vector<std::vector<int> > currentMaxima_coordinates;

	// init first cell
	unsigned int statesPerLevel0 = g->NodesPerLevel.at(startLevel_graph).size();
	assert((startZ_graph >= 0) && (startZ_graph < (int)statesPerLevel0));


	// parameters
	// threshold_for_filtering: will remove all cells from NW table in a given diagonal which have value > 15 difference from maximum
	int threshold_for_filtering = 15;
	int maximum_steps_nonIncrease = 40;


	std::set<std::string> achieved_complete_sequence_alignments;

	for(unsigned int stateI = 0; stateI < statesPerLevel0; stateI++)
	{
		if((int)stateI == startZ_graph)
		{
			scores[startLevel_graph][start_sequence][stateI].D = 0;
			scores_backtrace[startLevel_graph][start_sequence][stateI] = mScore_backtrace();
		}
		else
		{
			scores[startLevel_graph][start_sequence][stateI].D = minusInfinity;
		}

		scores[startLevel_graph][start_sequence][stateI].GraphGap = minusInfinity;
		scores[startLevel_graph][start_sequence][stateI].SequenceGap = minusInfinity;

		if((int)stateI == startZ_graph)
		{
			std::vector<int> existingCoordinates;
			existingCoordinates.push_back(startLevel_graph);
			existingCoordinates.push_back(start_sequence);
			existingCoordinates.push_back(stateI);
			m1_diagonal.push_back(existingCoordinates);
			currentMaxima_coordinates.push_back(existingCoordinates);
		}
	}

	std::map<NWPath*, std::pair<double, std::vector<int> > > hit_NW_paths;

//	std::cerr << "fullNeedleman_diagonal_extension:\n";
//	std::cerr << "\tmin_seqI: " << min_seqI << "\n";
//	std::cerr << "\tmin_seqI: " << min_seqI << "\n";
//	std::cerr << "\tdirectionPositive: " << directionPositive << "\n" << std::flush;

	int lastMaximumIncrease_at_diagonalI = 0;


	for(int diagonalI = 1; diagonalI <= diagonals; diagonalI++)
	{

		if(verbose)
		{
			std::cout << "\t diagonalI " << diagonalI << "/" << diagonals << ".\n" << std::flush;
		}


//		int scores_size = 0;
//		for(std::map<int, std::map<int, std::map<int, mScore>> >::iterator it1 = scores.begin(); it1 != scores.end(); it1++)
//		{
//			for(std::map<int, std::map<int, mScore>>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++)
//			{
//				for(std::map<int, mScore>::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++)
//				{
//					scores_size++;
//				}
//			}
//		}
		//std::cerr << "\t\tdiagonalI = " << diagonalI << " => scores_size: " << scores_size << "\n" << std::flush;

		if((diagonalI - lastMaximumIncrease_at_diagonalI) > maximum_steps_nonIncrease)
		{
			break;
		}

		std::map<int, std::map<int, std::map<int, mScore_alternatives > > >  thisDiagonal;
		std::map<int, std::map<int, std::map<int, mScore_backtrace_alternatives > > > thisDiagonal_backtrace;

		if(verbose)
			std::cout << "\t\tfrom m-2 diagonal" << "\n" << std::flush;

		// extend from m-2 diagonal
		for(int m2I = 0; m2I < (int)m2_diagonal.size(); m2I++)
		{
			std::vector<int>& previous_coordinates = m2_diagonal.at(m2I);

			int previous_levelI = previous_coordinates.at(0);
			int previous_seqI = previous_coordinates.at(1);
			int previous_stateI = previous_coordinates.at(2);

			int next_levelI = previous_coordinates.at(0) + (directionPositive ? 1 : -1);
			int next_seqI = previous_coordinates.at(1) + (directionPositive ? 1 : -1);

			if((next_levelI > max_levelI) || (next_seqI > max_seqI))
				continue;

			if((next_levelI < min_levelI) || (next_seqI < min_seqI))
				continue;

			std::string sequenceEmission = (directionPositive ? sequence.substr(previous_seqI, 1) : sequence.substr(previous_seqI-1, 1));

			std::vector<std::pair<int, Edge*> > nextZs = (directionPositive ? _graph_get_next_z_values_and_edges(previous_levelI, previous_stateI) : _graph_get_previous_z_values_and_edges(previous_levelI, previous_stateI));
			assert(nextZs.size() > 0);

			for(unsigned int zI = 0; zI < nextZs.size(); zI++)
			{
				std::pair<int, Edge*>& thisZjump = nextZs.at(zI);
				int next_stateI = thisZjump.first;

				std::string edgeEmission = thisZjump.second->getEmission();
				assert(edgeEmission.size() == 1);

				double score_MatchMismatch = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + ((edgeEmission == sequenceEmission) ? S_match : S_mismatch);

				backtraceStep_affine backtrack_MatchMismatch;
				backtrack_MatchMismatch.x = previous_levelI;
				backtrack_MatchMismatch.y = previous_seqI;
				backtrack_MatchMismatch.z = previous_stateI;
				backtrack_MatchMismatch.usedEdge = thisZjump.second;
				backtrack_MatchMismatch.sourceMatrix = 0;

				thisDiagonal[next_levelI][next_seqI][next_stateI].D.push_back(score_MatchMismatch);
				thisDiagonal_backtrace[next_levelI][next_seqI][next_stateI].D.push_back(backtrack_MatchMismatch);
			}
		}

		if(verbose)
			std::cout << "\t\tfrom m-1 diagonal" << "\n" << std::flush;

		// extend from m-1 diagonal
		for(int m1I = 0; m1I < (int)m1_diagonal.size(); m1I++)
		{
			std::vector<int>& previous_coordinates = m1_diagonal.at(m1I);

			int previous_levelI = previous_coordinates.at(0);
			int previous_seqI = previous_coordinates.at(1);
			int previous_stateI = previous_coordinates.at(2);

			// gap in graph
			int gapInGraph_next_levelI = previous_coordinates.at(0);
			int gapInGraph_next_seqI = previous_coordinates.at(1) + (directionPositive ? 1 : -1);
			if	(
					(directionPositive && (gapInGraph_next_levelI <= max_levelI) && (gapInGraph_next_seqI <= max_seqI)) ||
					((! directionPositive) && (gapInGraph_next_levelI >= min_levelI) && (gapInGraph_next_seqI >= min_seqI))
				)
			{
				double score_gapInGraph_open = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_openGap + S_extendGap;
				if((! directionPositive) && (gapInGraph_next_levelI == 0))
				{
					// score_gapInGraph_open = minusInfinity;
				}
				int gapInGraph_next_stateI = previous_stateI;

				backtraceStep_affine backtrack_gapInGraph_open;
				backtrack_gapInGraph_open.x = previous_levelI;
				backtrack_gapInGraph_open.y = previous_seqI;
				backtrack_gapInGraph_open.z = previous_stateI;
				backtrack_gapInGraph_open.usedEdge = 0;
				backtrack_gapInGraph_open.sourceMatrix = 0;

				thisDiagonal[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(score_gapInGraph_open);
				thisDiagonal_backtrace[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(backtrack_gapInGraph_open);

				double score_gapInGraph_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).GraphGap + S_extendGap;
				if((! directionPositive) && (gapInGraph_next_levelI == 0))
				{
					// score_gapInGraph_extend = minusInfinity;
				}

				backtraceStep_affine backtrack_gapInGraph_extend;
				backtrack_gapInGraph_extend.x = previous_levelI;
				backtrack_gapInGraph_extend.y = previous_seqI;
				backtrack_gapInGraph_extend.z = previous_stateI;
				backtrack_gapInGraph_extend.usedEdge = 0;
				backtrack_gapInGraph_extend.sourceMatrix = 1;

				thisDiagonal[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(score_gapInGraph_extend);
				thisDiagonal_backtrace[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(backtrack_gapInGraph_extend);
			}

			// gap in sequence
			int gapInSequence_next_levelI = previous_coordinates.at(0) + (directionPositive ? 1 : -1);
			int gapInSequence_next_seqI = previous_coordinates.at(1);

			if	(
					(directionPositive && (gapInSequence_next_levelI <= max_levelI) && (gapInSequence_next_seqI <= max_seqI)) ||
					((! directionPositive) && (gapInSequence_next_levelI >= min_levelI) && (gapInSequence_next_seqI >= min_seqI))
				)
			{
				std::vector<std::pair<int, Edge*> > nextZs = (directionPositive ? _graph_get_next_z_values_and_edges(previous_levelI, previous_stateI) : _graph_get_previous_z_values_and_edges(previous_levelI, previous_stateI));
				assert(nextZs.size() > 0);

				for(unsigned int zI = 0; zI < nextZs.size(); zI++)
				{
					std::pair<int, Edge*>& thisZjump = nextZs.at(zI);
					int next_stateI = thisZjump.first;

					std::string edgeEmission = thisZjump.second->getEmission();
					assert(edgeEmission.size() == 1);

					// open sequence gap

					double score_gapInSequence_open = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_openGap + S_extendGap;
					if(edgeEmission == "_")
					{
						score_gapInSequence_open = minusInfinity;
					}
					if((! directionPositive) && (gapInSequence_next_seqI == 0))
					{
						// score_gapInSequence_open = minusInfinity;
					}

					backtraceStep_affine backtrack_gapInSequence_open;
					backtrack_gapInSequence_open.x = previous_levelI;
					backtrack_gapInSequence_open.y = previous_seqI;
					backtrack_gapInSequence_open.z = previous_stateI;
					backtrack_gapInSequence_open.usedEdge = thisZjump.second;
					backtrack_gapInSequence_open.sourceMatrix = 0;

					thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(score_gapInSequence_open);
					thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(backtrack_gapInSequence_open);

					// extend sequence gap

					double score_gapInSequence_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap + S_extendGap;
//					if(edgeEmission == "_")
//					{
//						score_gapInSequence_extend = minusInfinity;
//					}
					if(edgeEmission == "_")
					{
						if(scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap == minusInfinity)
						{
							score_gapInSequence_extend = minusInfinity;
						}
						else
						{
							score_gapInSequence_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap +  S_graphGap;
						}
					}
					if((! directionPositive) && (gapInSequence_next_seqI == 0))
					{
						// score_gapInSequence_extend = minusInfinity;
					}

					backtraceStep_affine backtrack_gapInSequence_extend;
					backtrack_gapInSequence_extend.x = previous_levelI;
					backtrack_gapInSequence_extend.y = previous_seqI;
					backtrack_gapInSequence_extend.z = previous_stateI;
					backtrack_gapInSequence_extend.usedEdge = thisZjump.second;
					backtrack_gapInSequence_extend.sourceMatrix = 2;

					thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(score_gapInSequence_extend);
					thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(backtrack_gapInSequence_extend);

					// non-affine sequence gap
					if(edgeEmission == "_")
					{
						double score_gapInSequence_nonAffine = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_graphGap;

						backtraceStep_affine backtrack_gapInSequence_nonAffine;
						backtrack_gapInSequence_nonAffine.x = previous_levelI;
						backtrack_gapInSequence_nonAffine.y = previous_seqI;
						backtrack_gapInSequence_nonAffine.z = previous_stateI;
						backtrack_gapInSequence_nonAffine.usedEdge = thisZjump.second;
						backtrack_gapInSequence_nonAffine.sourceMatrix = 0;

						thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].D.push_back(score_gapInSequence_nonAffine);
						thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].D.push_back(backtrack_gapInSequence_nonAffine);
					}
				}
			}

			// gap in sequence, jump
			std::vector<std::pair<int, std::pair<int, Edge*>>> gapJumps = (directionPositive ? _graph_get_jumpNext_x_and_z_values_and_edges(previous_levelI, previous_stateI) : _graph_get_jumpPrevious_x_and_z_values_and_edges(previous_levelI, previous_stateI));

			for(unsigned int jumpI = 0; jumpI < gapJumps.size(); jumpI++)
			{
				const std::pair<int, std::pair<int, Edge*>>& thisJump = gapJumps.at(jumpI);
				int jump_levelI = thisJump.first;
				int jump_seqI = previous_seqI;
				int jump_stateI = thisJump.second.first;
				Edge* jump_edge = thisJump.second.second;

				if	(
						(directionPositive && (jump_levelI <= max_levelI) && (jump_seqI <= max_seqI)) ||
						((! directionPositive) && (jump_levelI >= min_levelI) && (jump_seqI >= min_seqI))
					)
				{
					int jump_length = g->completedGapEdgePaths.at(g->pseudoEdges_correspondingToGapEdgePaths.at(jump_edge)).size();
					double score_gapInSequence_nonAffine = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + (jump_length * S_graphGap);

					backtraceStep_affine backtrack_gapInSequence_nonAffine;
					backtrack_gapInSequence_nonAffine.x = previous_levelI;
					backtrack_gapInSequence_nonAffine.y = previous_seqI;
					backtrack_gapInSequence_nonAffine.z = previous_stateI;
					backtrack_gapInSequence_nonAffine.usedEdge = jump_edge;
					backtrack_gapInSequence_nonAffine.sourceMatrix = 0;

					thisDiagonal[jump_levelI][jump_seqI][jump_stateI].D.push_back(score_gapInSequence_nonAffine);
					thisDiagonal_backtrace[jump_levelI][jump_seqI][jump_stateI].D.push_back(backtrack_gapInSequence_nonAffine);

				}
			}
		}

		if(verbose)
			std::cout << "\t\tmaximal" << "\n" << std::flush;

		// call maxima for this diagonal
		std::vector<std::vector<int> > m_thisDiagonal;
		for(std::map<int, std::map<int, std::map<int, mScore_alternatives > > >::iterator diagIt = thisDiagonal.begin(); diagIt != thisDiagonal.end(); diagIt++)
		{
			int levelI = diagIt->first;
			for(std::map<int, std::map<int, mScore_alternatives > >::iterator diagIt2 = diagIt->second.begin(); diagIt2 != diagIt->second.end(); diagIt2++)
			{
				int seqI = diagIt2->first;
				for(std::map<int, mScore_alternatives >::iterator diagIt3 = diagIt2->second.begin(); diagIt3 != diagIt2->second.end(); diagIt3++)
				{
					int stateI = diagIt3->first;

					// call maxima for GraphGap
					std::vector<double>& scores_GraphGap = diagIt2->second.at(stateI).GraphGap;
					std::vector<backtraceStep_affine>& scores_GraphGap_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).GraphGap;
					double selectedScore_GraphGap;
					backtraceStep_affine selectedStep_GraphGap;

					if(scores_GraphGap.size() > 0)
					{
						std::pair<double, unsigned int> maxScore_GraphGap = Utilities::findVectorMaxP_nonCritical(scores_GraphGap, &(rng_seeds.at(omp_get_thread_num())));
						selectedScore_GraphGap = maxScore_GraphGap.first;
						selectedStep_GraphGap = scores_GraphGap_bt.at(maxScore_GraphGap.second);
					}
					else
					{
						selectedScore_GraphGap = minusInfinity;
					}

					// call maxima for SequenceGap
					std::vector<double>& scores_SequenceGap = diagIt2->second.at(stateI).SequenceGap;
					std::vector<backtraceStep_affine>& scores_SequenceGap_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).SequenceGap;

					double selectedScore_SequenceGap;
					backtraceStep_affine selectedStep_SequenceGap;
					if(scores_SequenceGap.size() > 0)
					{
						std::pair<double, unsigned int> maxScore_SequenceGap = Utilities::findVectorMaxP_nonCritical(scores_SequenceGap, &(rng_seeds.at(omp_get_thread_num())));
						selectedScore_SequenceGap = maxScore_SequenceGap.first;
						selectedStep_SequenceGap = scores_SequenceGap_bt.at(maxScore_SequenceGap.second);
					}
					else
					{
						selectedScore_SequenceGap = minusInfinity;
					}

					// two additional steps for D, jumping from the two gap matrices

					double score_intoD_fromGraphGap = selectedScore_GraphGap;
					backtraceStep_affine step_intoD_fromGraphGap;
					step_intoD_fromGraphGap.x = levelI;
					step_intoD_fromGraphGap.y = seqI;
					step_intoD_fromGraphGap.z = stateI;
					step_intoD_fromGraphGap.sourceMatrix = 1;
					step_intoD_fromGraphGap.usedEdge = 0;
					thisDiagonal.at(levelI).at(seqI).at(stateI).D.push_back(score_intoD_fromGraphGap);
					thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D.push_back(step_intoD_fromGraphGap);

					double score_intoD_fromSequenceGap = selectedScore_SequenceGap;
					backtraceStep_affine step_intoD_fromSequenceGap;
					step_intoD_fromSequenceGap.x = levelI;
					step_intoD_fromSequenceGap.y = seqI;
					step_intoD_fromSequenceGap.z = stateI;
					step_intoD_fromSequenceGap.sourceMatrix = 2;
					step_intoD_fromSequenceGap.usedEdge = 0;
					thisDiagonal.at(levelI).at(seqI).at(stateI).D.push_back(score_intoD_fromSequenceGap);
					thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D.push_back(step_intoD_fromSequenceGap);

					// find final maximum for D
					std::vector<double>& scores_D = diagIt2->second.at(stateI).D;
					std::vector<backtraceStep_affine>& scores_D_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D;
					assert(scores_D.size() == scores_D_bt.size());

					std::pair<double, unsigned int> maxScore_D = Utilities::findVectorMaxP_nonCritical(scores_D, &(rng_seeds.at(omp_get_thread_num())));

					bool blockOutCell = false;
					if(directionPositive)
					{
						std::set<NWEdge*> emanatingNWEdges = blockedPathsTable->getEdgesEmanatingFrom(levelI, seqI, stateI);
						if((emanatingNWEdges.size() > 0) && (! returnGlobalScore))
						{
							for(std::set<NWEdge*>::iterator edgeIt = emanatingNWEdges.begin(); edgeIt != emanatingNWEdges.end(); edgeIt++)
							{
								NWEdge* e = *edgeIt;
								NWPath* correspondingPath = e->path;
								if(hit_NW_paths.count(correspondingPath) == 0)
								{
									std::pair<double, std::vector<int> > pathEntry;
									pathEntry.first = maxScore_D.first;
									pathEntry.second.push_back(levelI);
									pathEntry.second.push_back(seqI);
									pathEntry.second.push_back(stateI);
									hit_NW_paths[correspondingPath] = pathEntry;
								}
							}
							blockOutCell = true;
						}
					}
					else
					{
						std::set<NWEdge*> incomingNWEdges = blockedPathsTable->getEdgesGoingInto(levelI, seqI, stateI);
						if((incomingNWEdges.size() > 0) && (! returnGlobalScore))
						{
							for(std::set<NWEdge*>::iterator edgeIt = incomingNWEdges.begin(); edgeIt != incomingNWEdges.end(); edgeIt++)
							{
								NWEdge* e = *edgeIt;
								NWPath* correspondingPath = e->path;
								if(hit_NW_paths.count(correspondingPath) == 0)
								{
									std::pair<double, std::vector<int> > pathEntry;
									pathEntry.first = maxScore_D.first;
									pathEntry.second.push_back(levelI);
									pathEntry.second.push_back(seqI);
									pathEntry.second.push_back(stateI);
									hit_NW_paths[correspondingPath] = pathEntry;
								}
							}
							blockOutCell = true;
						}
					}

					if(blockOutCell)
					{
						double scoreThatWillBeDeleted = maxScore_D.first;
						std::vector<int> thisCoordinates;
						thisCoordinates.push_back(levelI);
						thisCoordinates.push_back(seqI);
						thisCoordinates.push_back(stateI);
						if(scoreThatWillBeDeleted == currentMaximum)
						{
							currentMaxima_coordinates.push_back(thisCoordinates);
							lastMaximumIncrease_at_diagonalI = diagonalI;
						}
						else if(scoreThatWillBeDeleted > currentMaximum)
						{
							currentMaximum = scoreThatWillBeDeleted;
							currentMaxima_coordinates.clear();
							currentMaxima_coordinates.push_back(thisCoordinates);
							lastMaximumIncrease_at_diagonalI = diagonalI;
						}

						if((scores.count(levelI) == 0) || (scores.at(levelI).count(seqI) == 0) || (scores.at(levelI).at(seqI).count(stateI) == 0))
						{
							scores[levelI][seqI][stateI].D = minusInfinity;
							scores_backtrace[levelI][seqI][stateI].D = scores_D_bt.at(maxScore_D.second);

							assert(scores_backtrace[levelI][seqI][stateI].D.x != -1);

							scores[levelI][seqI][stateI].GraphGap = minusInfinity;
							scores_backtrace[levelI][seqI][stateI].GraphGap = selectedStep_GraphGap;

							scores[levelI][seqI][stateI].SequenceGap = minusInfinity;
							scores_backtrace[levelI][seqI][stateI].SequenceGap = selectedStep_SequenceGap;
						}
					}
					else
					{
						if(maxScore_D.first >= diagonal_stop_threshold)
						{
							bool newEntry = ( (scores.count(levelI) == 0) || (scores.at(levelI).count(seqI) == 0) || (scores.at(levelI).at(seqI).count(stateI) == 0) );
							bool overwrittenEntry = false;
							if(		newEntry ||
									(scores.at(levelI).at(seqI).at(stateI).D < maxScore_D.first)
							)
							{
								overwrittenEntry = ! newEntry;
								scores[levelI][seqI][stateI].D = maxScore_D.first;
								scores_backtrace[levelI][seqI][stateI].D = scores_D_bt.at(maxScore_D.second);
								assert(scores_backtrace[levelI][seqI][stateI].D.x != -1);
							}

							if(		newEntry ||
									(scores.at(levelI).at(seqI).at(stateI).GraphGap < selectedScore_GraphGap)
							)
							{
								overwrittenEntry = ! newEntry;
								scores[levelI][seqI][stateI].GraphGap = selectedScore_GraphGap;
								scores_backtrace[levelI][seqI][stateI].GraphGap = selectedStep_GraphGap;
							}

							if(		newEntry ||
									(scores.at(levelI).at(seqI).at(stateI).SequenceGap < selectedScore_SequenceGap)
							)
							{
								overwrittenEntry = ! newEntry;
								scores[levelI][seqI][stateI].SequenceGap = selectedScore_SequenceGap;
								scores_backtrace[levelI][seqI][stateI].SequenceGap = selectedStep_SequenceGap;
							}


							if(preferSequenceCompleAlignments)
							{
								std::string id_levelI_stateI = Utilities::ItoStr(levelI) + "/" + Utilities::ItoStr(stateI);
								if(directionPositive)
								{
									if(seqI == max_seqI)
									{
										achieved_complete_sequence_alignments.insert(id_levelI_stateI);
									}
								}
								else
								{
									if(seqI == min_seqI)
									{
										achieved_complete_sequence_alignments.insert(id_levelI_stateI);
									}
								}
							}

							std::vector<int> thisCoordinates;
							thisCoordinates.push_back(levelI);
							thisCoordinates.push_back(seqI);
							thisCoordinates.push_back(stateI);
							m_thisDiagonal.push_back(thisCoordinates);

							backtraceStep_affine oneRealStepBackwards = scores_backtrace[levelI][seqI][stateI].D;
							while((oneRealStepBackwards.x == levelI) && (oneRealStepBackwards.y == seqI))
							{
								assert(oneRealStepBackwards.sourceMatrix != 0);
								if(oneRealStepBackwards.sourceMatrix == 1)
								{
									oneRealStepBackwards = scores_backtrace[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].GraphGap;
								}
								else if(oneRealStepBackwards.sourceMatrix == 2)
								{
									oneRealStepBackwards = scores_backtrace[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].SequenceGap;
								}
								else
								{
									assert(1 == 0);
								}
							}
							int previousScore;
							if(oneRealStepBackwards.sourceMatrix == 0)
							{
								previousScore = scores[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].D;
							}
							else if(oneRealStepBackwards.sourceMatrix == 1)
							{
								previousScore = scores[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].GraphGap;
							}
							else if(oneRealStepBackwards.sourceMatrix == 2)
							{
								previousScore = scores[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].SequenceGap;
							}
							else
							{
								assert( 1 == 0 );
							}
							int scoreDifference = maxScore_D.first - previousScore;

							if(maxScore_D.first == currentMaximum)
							{
								if(scoreDifference != 0)
								{
									currentMaxima_coordinates.push_back(thisCoordinates);
									lastMaximumIncrease_at_diagonalI = diagonalI;
								}
							}
							else if(maxScore_D.first > currentMaximum)
							{
								currentMaximum = maxScore_D.first;
								currentMaxima_coordinates.clear();
								currentMaxima_coordinates.push_back(thisCoordinates);
								lastMaximumIncrease_at_diagonalI = diagonalI;
							}

							if(overwrittenEntry)
							{
								lastMaximumIncrease_at_diagonalI = diagonalI;
							}
						}
						else
						{
							// assert(1 == 0);
						}
					}
				}
			}
		}

		if(verbose)
			std::cout << "\t\tfiltering" << "\n" << std::flush;

		std::vector<std::vector<int> > m_thisDiagonal_filtered;
		if(m_thisDiagonal.size() > 0)
		{
			double max;
			for(unsigned int i = 0; i < m_thisDiagonal.size(); i++)
			{
				std::vector<int> coordinates = m_thisDiagonal.at(i);
				double S = scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D;
				if((i == 0) || (max < S))
				{
					max = S;
				}
			}

			for(unsigned int i = 0; i < m_thisDiagonal.size(); i++)
			{
				std::vector<int> coordinates = m_thisDiagonal.at(i);
				double S = scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D;
				assert(S <= max);
				if((max - S) <= threshold_for_filtering)
				{
					m_thisDiagonal_filtered.push_back(coordinates);
				}
			}

			m_thisDiagonal = m_thisDiagonal_filtered;
		}

		m2_diagonal = m1_diagonal;
		m1_diagonal = m_thisDiagonal;
	}

	std::vector<localExtension_pathDescription> forReturn;
	auto backtraceFrom = [&](int start_x, int start_y, int start_z, double StartScore) {
		int backtrace_x = start_x;
		int backtrace_y = start_y;
		int backtrace_z = start_z;
		int backtrace_matrix = 0;

		std::string reconstructed_graph;
		std::string reconstructed_sequence;
		std::vector<int> reconstructed_graph_levels;
		std::vector<std::vector<int>> edge_coordinates;
		std::vector<Edge*> utilizedEdges;

		std::vector<int> startCoordinates;
		startCoordinates.push_back(start_x);
		startCoordinates.push_back(start_y);
		startCoordinates.push_back(start_z);
		edge_coordinates.push_back(startCoordinates);

		if(verbose)
		{
			std::cout << "\tbacktraceFrom() called.\n";
			std::cout << "\t\t" << "start_x" << ": " << start_x << "\n";
			std::cout << "\t\t" << "start_y" << ": " << start_y << "\n";
			std::cout << "\t\t" << "start_z" << ": " << start_z << "\n";
			std::cout << std::flush;
		}
		while((backtrace_x != startLevel_graph) || (backtrace_y != start_sequence))
		{
			if(verbose)
			{
				std::cout << "\t\tbacktrace_x: " << backtrace_x << "\n";
				std::cout << "\t\tbacktrace_y: " << backtrace_y << "\n";
				std::cout << "\t\tbacktrace_z: " << backtrace_z << "\n";
				std::cout << "\t\tselected matrix: " << backtrace_matrix << "\n" << std::flush;
			}

			if(directionPositive)
			{
				assert((backtrace_x >= startLevel_graph) && (backtrace_x <= ((int)levels - 1)));
				assert((backtrace_y >= start_sequence) && (backtrace_y <= (int)sequenceLength));
				assert((backtrace_z >= 0) && (backtrace_z <= (int)g->NodesPerLevel.at(backtrace_x).size()));
				assert((backtrace_matrix >= 0) && (backtrace_matrix <= 2));
			}
			else
			{
				assert((backtrace_x <= startLevel_graph) && (backtrace_x >= 0));
				assert((backtrace_y <= start_sequence) && (backtrace_y >= 0));
				assert((backtrace_z >= 0) && (backtrace_z <= (int)g->NodesPerLevel.at(backtrace_x).size()));
				assert((backtrace_matrix >= 0) && (backtrace_matrix <= 2));
			}

			double Score;
			backtraceStep_affine step;
			if(backtrace_matrix == 0)
			{
				Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).D;
				step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).D;
			}
			else if(backtrace_matrix == 1)
			{
				Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).GraphGap;
				step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).GraphGap;
			}
			else if(backtrace_matrix == 2)
			{
				Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).SequenceGap;
				step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).SequenceGap;
			}

			if((backtrace_x == start_x) && (backtrace_y == start_y) && (backtrace_z == start_z))
			{
				Score = StartScore;
			}

			if(verbose)
			{
				std::cout << "\t\tScore: " << Score << "\n\n" << std::flush;
			}

			std::string edgeEmission;
			if((step.usedEdge != 0) && (g->pseudoEdges_correspondingToGapEdgePaths.count(step.usedEdge) == 0))
			{
				edgeEmission = step.usedEdge->getEmission();
				assert(edgeEmission.size() == 1);
			}
			std::string sequenceEmission;
			if((backtrace_y >= 1) && (directionPositive))
			{
				sequenceEmission = sequence.substr(backtrace_y - 1, 1);
			}
			if((backtrace_y < max_seqI) && (! directionPositive))
			{
				sequenceEmission = sequence.substr(backtrace_y, 1);
			}

			int next_x = step.x;
			int next_y = step.y;
			int next_z = step.z;
			int next_matrix = step.sourceMatrix;

			bool dontAddCoordinates = false; // if we jump from one matrix to the other without changing coordinates, we don't store the coordinates...
			if(g->pseudoEdges_correspondingToGapEdgePaths.count(step.usedEdge) == 0)
			{
				if(directionPositive)
				{
					if((next_x == (backtrace_x - 1)) && (next_y == (backtrace_y - 1)))
					{
						// match or mismatch
						assert(step.usedEdge != 0);
						reconstructed_graph.append(edgeEmission);
						reconstructed_graph_levels.push_back(backtrace_x - 1);
						reconstructed_sequence.append(sequenceEmission);
						utilizedEdges.push_back(step.usedEdge);
					}
					else if((next_x == backtrace_x) && (next_y == (backtrace_y - 1)))
					{
						// gap in graph
						reconstructed_graph.append("_");
						reconstructed_graph_levels.push_back(-1);
						reconstructed_sequence.append(sequenceEmission);
						utilizedEdges.push_back(0);

					}
					else if((next_x == (backtrace_x - 1)) && (next_y == backtrace_y))
					{
						// gap in sequence
						assert(step.usedEdge != 0);
						reconstructed_graph.append(edgeEmission);
						reconstructed_graph_levels.push_back(backtrace_x - 1);
						reconstructed_sequence.append("_");
						utilizedEdges.push_back(step.usedEdge);
					}
					else
					{
						dontAddCoordinates = true;
						assert((backtrace_x == next_x) && (backtrace_y == next_y) && (backtrace_z == next_z) && (next_matrix != backtrace_matrix));
					}
				}
				else
				{
					if((next_x == (backtrace_x + 1)) && (next_y == (backtrace_y + 1)))
					{
						// match or mismatch
						assert(step.usedEdge != 0);
						reconstructed_graph.append(edgeEmission);
						reconstructed_graph_levels.push_back(backtrace_x);
						reconstructed_sequence.append(sequenceEmission);
						utilizedEdges.push_back(step.usedEdge);
					}
					else if((next_x == backtrace_x) && (next_y == (backtrace_y + 1)))
					{
						// gap in graph
						reconstructed_graph.append("_");
						reconstructed_graph_levels.push_back(-1);
						reconstructed_sequence.append(sequenceEmission);
						utilizedEdges.push_back(0);
					}
					else if((next_x == (backtrace_x + 1)) && (next_y == backtrace_y))
					{
						// gap in sequence
						assert(step.usedEdge != 0);
						reconstructed_graph.append(edgeEmission);
						reconstructed_graph_levels.push_back(backtrace_x);
						reconstructed_sequence.append("_");
						utilizedEdges.push_back(step.usedEdge);
					}
					else
					{
						dontAddCoordinates = true;
						assert((backtrace_x == next_x) && (backtrace_y == next_y) && (backtrace_z == next_z) && (next_matrix != backtrace_matrix));
					}
				}
			}
			else
			{
				std::vector<Edge*> edgePath = g->completedGapEdgePaths.at(g->pseudoEdges_correspondingToGapEdgePaths.at(step.usedEdge));
				int L = edgePath.size();
				std::string gapString;
				gapString.resize(L, '_');
				std::vector<int> graph_levels;
				for(unsigned int eI = 0; eI < edgePath.size(); eI++)
				{
					Edge* e = edgePath.at(eI);
					assert(e->getEmission() == "_");
					int level = e->From->level;
					graph_levels.push_back(level);
				}

				if(directionPositive)
				{
					std::reverse(graph_levels.begin(), graph_levels.end());
					std::reverse(edgePath.begin(), edgePath.end());
				}

				reconstructed_graph_levels.insert(reconstructed_graph_levels.end(), graph_levels.begin(), graph_levels.end());
				reconstructed_graph.append(gapString);
				reconstructed_sequence.append(gapString);
				utilizedEdges.insert(utilizedEdges.end(), edgePath.begin(), edgePath.end());
			}

			if(verbose)
			{
				std::cout << "\t\t\t" << "Emission: " << edgeEmission << "\n";
				std::cout << "\t\t\t" << "next_x: " << next_x << "\n";
				std::cout << "\t\t\t" << "next_y: " << next_y << "\n";
				std::cout << "\t\t\t" << "next_z: " << next_z << "\n";
				std::cout << "\t\t\t" << "next_matrix: " << next_matrix << "\n\n" << std::flush;
			}

			backtrace_x = next_x;
			backtrace_y = next_y;
			backtrace_z = next_z;
			backtrace_matrix = next_matrix;

			if( ! dontAddCoordinates)
			{
				std::vector<int> nextCoordinates;
				nextCoordinates.push_back(next_x);
				nextCoordinates.push_back(next_y);
				nextCoordinates.push_back(next_z);
				edge_coordinates.push_back(nextCoordinates);
			}
		}

		if(directionPositive)
		{
			std::reverse(reconstructed_graph.begin(), reconstructed_graph.end());
			std::reverse(reconstructed_graph_levels.begin(), reconstructed_graph_levels.end());
			std::reverse(reconstructed_sequence.begin(), reconstructed_sequence.end());
			std::reverse(utilizedEdges.begin(), utilizedEdges.end());
			std::reverse(edge_coordinates.begin(), edge_coordinates.end());
		}

		assert(reconstructed_graph.size() == reconstructed_sequence.size());
		assert(reconstructed_graph_levels.size() == reconstructed_graph.size());

		localExtension_pathDescription pathReturn;
		pathReturn.Score = StartScore;
		pathReturn.usedEdges = utilizedEdges;
		pathReturn.coordinates = edge_coordinates;
		pathReturn.alignedSequence = reconstructed_sequence;
		pathReturn.alignedGraph = reconstructed_graph;
		pathReturn.alignedGraph_levels = reconstructed_graph_levels;

		forReturn.push_back(pathReturn);
	};

	if(verbose)
	{
			/*
		for(unsigned int x = 231412; x <= 231418; x++)
		{
			for(unsigned int y = 0; y <= 5; y++)
			{
				std::cerr << "x = " << x << ", y = " << y << ": ";
				if((scores.count(x) == 0) || (scores.at(x).count(y) == 0) || (scores.at(x).at(y).size() == 0))
				{
					std::cerr << "Undefined.";
				}
				else
				{
					for(unsigned int z = 0; z < scores.at(x).at(y).size(); z++)
					{
						std::cerr << "\n\t z = " << z << ": " << scores.at(x).at(y).at(z).D;
					}
				}
				std::cerr << "\n\n" << std::flush;
			}
		}
		*/
	}

	if(! returnGlobalScore)
	{
		bool have_backtrace_sequenceComplete = false;

		int sequenceCompleteBacktrace_levelI;
		int sequenceCompleteBacktrace_seqI;
		int sequenceCompleteBacktrace_stateI;
		double sequenceCompleteBacktrace_score;


		if(preferSequenceCompleAlignments)
		{

			int coordinate_seqI = (directionPositive ? max_seqI : min_seqI);

			int _dbg_local_chainI = 0;
			if(0 && directionPositive && (_dbg_local_chainI == 0))
			{
				std::cout << "fullNeedleman_affine_diagonal_extension(..) called, direction " << directionPositive << ".\n" << std::flush;
				std::cout << "\t" << "_dbg_local_chainI" << ": " << _dbg_local_chainI << "\n";
				std::cout << "\t" << "start_sequence" << ": " << start_sequence << "\n";
				std::cout << "\t" << "maxPosition_sequence" << ": " << maxPosition_sequence << "\n";
				std::cout << "\t" << "startLevel_graph" << ": " << startLevel_graph << "\n";
				std::cout << "\t" << "maxLevel_graph" << ": " << maxLevel_graph << "\n";
				std::cout << "\t" << "achieved_complete_sequence_alignments.size()" << ": " << achieved_complete_sequence_alignments.size() << "\n";

				for(std::set<std::string>::iterator coordinateIt = achieved_complete_sequence_alignments.begin(); coordinateIt != achieved_complete_sequence_alignments.end(); coordinateIt++)
				{
					std::string coordinates = *coordinateIt;
					std::vector<std::string> coordinate_components = Utilities::split(coordinates, "/");
					assert(coordinate_components.size() == 2);
					int coordinate_levelI = Utilities::StrtoI(coordinate_components.at(0));
					int coordinate_stateI = Utilities::StrtoI(coordinate_components.at(1));

					double S = scores.at(coordinate_levelI).at(coordinate_seqI).at(coordinate_stateI).D;

					std::cout << "\t\t\tAlternative:\n";
					std::cout << "\t\t\t\t" << "coordinate_seqI" << ": " << coordinate_seqI << "\n";
					std::cout << "\t\t\t\t" << "coordinate_levelI" << ": " << coordinate_levelI << "\n";
					std::cout << "\t\t\t\t" << "coordinate_stateI" << ": " << coordinate_stateI << "\n";
					std::cout << "\t\t\t\t" << "S" << ": " << S << "\n";
				}

				std::cout << std::flush;
			}

			std::vector<std::string> sequenceComplete_backtraceCoordinates_maxScores;
			double maxScore;


			for(std::set<std::string>::iterator coordinateIt = achieved_complete_sequence_alignments.begin(); coordinateIt != achieved_complete_sequence_alignments.end(); coordinateIt++)
			{
				std::string coordinates = *coordinateIt;
				std::vector<std::string> coordinate_components = Utilities::split(coordinates, "/");
				assert(coordinate_components.size() == 2);
				int coordinate_levelI = Utilities::StrtoI(coordinate_components.at(0));
				int coordinate_stateI = Utilities::StrtoI(coordinate_components.at(1));

				double S = scores.at(coordinate_levelI).at(coordinate_seqI).at(coordinate_stateI).D;

				if((coordinateIt == achieved_complete_sequence_alignments.begin()) || (S > maxScore))
				{
					sequenceComplete_backtraceCoordinates_maxScores.clear();
					sequenceComplete_backtraceCoordinates_maxScores.push_back(coordinates);
					maxScore = S;
				}
				else if(S == maxScore)
				{
					sequenceComplete_backtraceCoordinates_maxScores.push_back(coordinates);
				}
			}

			// std::cout << "\t\t" << "sequenceComplete_backtraceCoordinates_maxScores.size()" << ": " << sequenceComplete_backtraceCoordinates_maxScores.size() << "\n" << std::flush;

			if(sequenceComplete_backtraceCoordinates_maxScores.size() > 0)
			{
				have_backtrace_sequenceComplete = true;

				int selectedIndex = Utilities::randomNumber_nonCritical(sequenceComplete_backtraceCoordinates_maxScores.size() - 1, &(rng_seeds.at(omp_get_thread_num())));
				assert(selectedIndex >= 0);
				assert(selectedIndex < (int)sequenceComplete_backtraceCoordinates_maxScores.size());
				std::string coordinates = sequenceComplete_backtraceCoordinates_maxScores.at(selectedIndex);

				std::vector<std::string> coordinate_components = Utilities::split(coordinates, "/");
				assert(coordinate_components.size() == 2);

				sequenceCompleteBacktrace_levelI = Utilities::StrtoI(coordinate_components.at(0));
				sequenceCompleteBacktrace_seqI = coordinate_seqI;
				sequenceCompleteBacktrace_stateI = Utilities::StrtoI(coordinate_components.at(1));
				sequenceCompleteBacktrace_score = maxScore;

			}

		}
		if(preferSequenceCompleAlignments && have_backtrace_sequenceComplete)
		{
			backtraceFrom(sequenceCompleteBacktrace_levelI, sequenceCompleteBacktrace_seqI, sequenceCompleteBacktrace_stateI, sequenceCompleteBacktrace_score);
		}
		else
		{
			if(currentMaximum > 0)
			{
				if(verbose)
				{
					std::cout << "\tMaximum " << currentMaximum << ", achieved at " << currentMaxima_coordinates.size() << " positions!\n" << std::flush;
				}

				for(int maximumI = 0; maximumI < (int)currentMaxima_coordinates.size(); maximumI++)
				{
					std::vector<int> coordinates = currentMaxima_coordinates.at(maximumI);
					if(scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D != minusInfinity)
					{
						if(verbose)
							std::cout << " - Start maximum backtrace from " <<coordinates.at(0) << ", " <<  coordinates.at(1) << ", " << coordinates.at(2) << "\n" << std::flush;

						backtraceFrom(coordinates.at(0), coordinates.at(1), coordinates.at(2), scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D);
					}
				}
			}

			if(verbose)
				std::cout << "Have hit " << hit_NW_paths.size() << " NW paths, backtrace independent of achieved value!\n" << std::flush;

			for(std::map<NWPath*, std::pair<double, std::vector<int> > >::iterator hitPathsIt = hit_NW_paths.begin(); hitPathsIt != hit_NW_paths.end(); hitPathsIt++)
			{
				// NWPath* path = hitPathsIt->first;
				double score = hitPathsIt->second.first;
				std::vector<int> coordinates = hitPathsIt->second.second;
				assert(scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D == minusInfinity);

				if(verbose)
					std::cout << " - Start NW path backtrace from " << coordinates.at(0) << ", " <<  coordinates.at(1) << ", " << coordinates.at(2) << "\n" << std::flush;

				backtraceFrom(coordinates.at(0), coordinates.at(1), coordinates.at(2), score);
			}
		}
	}
	else
	{

		if(directionPositive)
		{
			std::map<int, double> finalState_scores;
			for(std::map<int, mScore>::iterator altIt = scores.at(levels-1).at(sequenceLength).begin(); altIt != scores.at(levels-1).at(sequenceLength).end(); altIt++)
			{
				finalState_scores[altIt->first] = altIt->second.D;
			}
			std::pair<double, int> maxFinalState = Utilities::findIntMapMaxP_nonCritical(finalState_scores, &(rng_seeds.at(omp_get_thread_num())));

			backtraceFrom(levels - 1, sequenceLength, maxFinalState.second, maxFinalState.first);
		}
		else
		{
			std::map<int, double> finalState_scores;
			for(std::map<int, mScore>::iterator altIt = scores.at(0).at(0).begin(); altIt != scores.at(0).at(0).end(); altIt++)
			{
				finalState_scores[altIt->first] = altIt->second.D;

				if(verbose)
					std::cout << "Final state z-value " << altIt->first << ", D value " << scores.at(0).at(0).at(altIt->first).D << ", GraphGap: " <<  scores.at(0).at(0).at(altIt->first).GraphGap << ", SequenceGap: " <<  scores.at(0).at(0).at(altIt->first).SequenceGap << "\n" << std::flush;
			}
			std::pair<double, int> maxFinalState = Utilities::findIntMapMaxP_nonCritical(finalState_scores, &(rng_seeds.at(omp_get_thread_num())));

			backtraceFrom(0, 0, maxFinalState.second, maxFinalState.first);
		}
	}

	std::vector<mapper::reads::verboseSeedChain> vC;
	for(unsigned int i = 0; i < forReturn.size(); i++)
	{
		vC.push_back(forReturn.at(0).toVerboseSeedChain());
	}

	return vC;
	// return forReturn;
}


std::vector<mapper::reads::verboseSeedChain> extensionAligner::fullNeedleman_diagonal_extension(const std::string& sequence, int start_sequence, int startLevel_graph, int startZ_graph, int maxLevel_graph, int maxPosition_sequence, int diagonal_stop_threshold, VirtualNWTable_Unique* blockedPathsTable, bool directionPositive, bool returnGlobalScore, bool preferSequenceCompleAlignments) const
{
	assert(omp_get_num_threads() == (int)threads_initialized);

	assert(! ((returnGlobalScore) && (preferSequenceCompleAlignments)));

	assert(g != 0);
	if(directionPositive)
	{
		if(!((maxLevel_graph == -1) || (maxLevel_graph > startLevel_graph)))
		{
			std::cerr << "maxLevel_graph: " << maxLevel_graph << "\n";
			std::cerr << "startLevel_graph: " << startLevel_graph << "\n" << std::flush;
		}
		if(!((maxPosition_sequence == -1) || (maxPosition_sequence > start_sequence)))
		{
			std::cerr << "maxPosition_sequence: " << maxPosition_sequence << "\n";
			std::cerr << "start_sequence: " << start_sequence << "\n" << std::flush;
		}
		assert((maxLevel_graph == -1) || (maxLevel_graph > startLevel_graph));
		assert((maxPosition_sequence == -1) || (maxPosition_sequence > start_sequence));
	}
	else
	{
		assert((maxLevel_graph == -1) || (maxLevel_graph < startLevel_graph));
		assert((maxPosition_sequence == -1) || (maxPosition_sequence < start_sequence));
	}

	double minusInfinity = -1 * numeric_limits<double>::max();

	if(returnGlobalScore)
	{
		diagonal_stop_threshold = minusInfinity;
		assert(maxLevel_graph == -1);
	}

	class mScore {
	public:
		double D;
		double GraphGap;
		double SequenceGap;
	};

	class mScore_backtrace {
	public:
		backtraceStep_affine D;
		backtraceStep_affine GraphGap;
		backtraceStep_affine SequenceGap;
	};

	class mScore_alternatives {
	public:
		std::vector<double> D;
		std::vector<double> GraphGap;
		std::vector<double> SequenceGap;
	};

	class mScore_backtrace_alternatives {
	public:
		std::vector<backtraceStep_affine> D;
		std::vector<backtraceStep_affine> GraphGap;
		std::vector<backtraceStep_affine> SequenceGap;
	};

	std::map<int, std::map<int, std::map<int, mScore>> > scores;
	std::map<int, std::map<int, std::map<int, mScore_backtrace>> > scores_backtrace;

	std::vector<std::vector<int> > coordinates_for_backtracking;
	std::vector<std::vector<int> > m1_diagonal;
	std::vector<std::vector<int> > m2_diagonal;

	bool verbose = false;

	if(verbose)
	{
		std::cout << "fullNeedleman_affine_diagonal_extension(..) called, direction " << directionPositive << ".\n" << std::flush;
		std::cout << "\t" << "start_sequence" << ": " << start_sequence << "\n";
		std::cout << "\t" << "maxPosition_sequence" << ": " << maxPosition_sequence << "\n";
		std::cout << "\t" << "startLevel_graph" << ": " << startLevel_graph << "\n";
		std::cout << "\t" << "maxLevel_graph" << ": " << maxLevel_graph << "\n";
		std::cout << "\t" << "directionPositive" << ": " << directionPositive << "\n";
		std::cout << "\t" << "returnGlobalScore" << ": " << returnGlobalScore << "\n";
		std::cout << "\t" << "preferSequenceCompleAlignments" << ": " << preferSequenceCompleAlignments << "\n";

		if(directionPositive)
		{
			std::string coveredSequence = sequence.substr(start_sequence, maxPosition_sequence - start_sequence);
			std::cout << "\tSequeence: " << coveredSequence << "\n" << std::flush;
		}
		else
		{
			std::string coveredSequence = sequence.substr(maxPosition_sequence, start_sequence -  maxPosition_sequence);
			std::cout << "\tSequeence: " << coveredSequence << "\n" << std::flush;

		}

		std::cout << std::flush;
	}

	unsigned int levels = g->NodesPerLevel.size();
	unsigned int sequenceLength = sequence.length();
	int diagonals = sequenceLength + levels - 1;
	int max_levelI = levels - 1;
	int max_seqI = sequenceLength;
	int min_levelI = 0;
	int min_seqI = 0;

	if(maxLevel_graph != -1)
	{
		if(directionPositive)
		{
			max_levelI = maxLevel_graph;
		}
		else
		{
			min_levelI = maxLevel_graph;
		}
	}
	if(maxPosition_sequence != -1)
	{
		if(directionPositive)
		{
			max_seqI = maxPosition_sequence;
		}
		else
		{
			min_seqI = maxPosition_sequence;
		}
	}

	assert(startLevel_graph >= min_levelI);
	assert(start_sequence >= min_seqI);
	assert(startLevel_graph <= max_levelI);
	assert(start_sequence <= max_seqI);


	assert(min_levelI >= 0);
	assert(max_levelI <= (int)(levels - 1));
	assert(max_levelI > min_levelI);

	assert(min_seqI >= 0);
	assert(max_seqI <= (int)sequenceLength);
	assert(max_seqI > min_seqI);

	double currentMaximum = 0;
	std::vector<std::vector<int> > currentMaxima_coordinates;

	// init first cell
	unsigned int statesPerLevel0 = g->NodesPerLevel.at(startLevel_graph).size();
	assert((startZ_graph >= 0) && (startZ_graph < (int)statesPerLevel0));


	// parameters
	// threshold_for_filtering: will remove all cells from NW table in a given diagonal which have value > 15 difference from maximum
	int threshold_for_filtering = 15;
	int maximum_steps_nonIncrease = 40;


	std::set<std::string> achieved_complete_sequence_alignments;

	for(unsigned int stateI = 0; stateI < statesPerLevel0; stateI++)
	{
		if((int)stateI == startZ_graph)
		{
			scores[startLevel_graph][start_sequence][stateI].D = 0;
			scores_backtrace[startLevel_graph][start_sequence][stateI] = mScore_backtrace();
		}
		else
		{
			scores[startLevel_graph][start_sequence][stateI].D = minusInfinity;
		}

		scores[startLevel_graph][start_sequence][stateI].GraphGap = minusInfinity;
		scores[startLevel_graph][start_sequence][stateI].SequenceGap = minusInfinity;

		if((int)stateI == startZ_graph)
		{
			std::vector<int> existingCoordinates;
			existingCoordinates.push_back(startLevel_graph);
			existingCoordinates.push_back(start_sequence);
			existingCoordinates.push_back(stateI);
			m1_diagonal.push_back(existingCoordinates);
			currentMaxima_coordinates.push_back(existingCoordinates);
		}
	}

	std::map<NWPath*, std::pair<double, std::vector<int> > > hit_NW_paths;

//	std::cerr << "fullNeedleman_diagonal_extension:\n";
//	std::cerr << "\tmin_seqI: " << min_seqI << "\n";
//	std::cerr << "\tmin_seqI: " << min_seqI << "\n";
//	std::cerr << "\tdirectionPositive: " << directionPositive << "\n" << std::flush;

	int lastMaximumIncrease_at_diagonalI = 0;


	for(int diagonalI = 1; diagonalI <= diagonals; diagonalI++)
	{

		if(verbose)
		{
			std::cout << "\t diagonalI " << diagonalI << "/" << diagonals << ".\n" << std::flush;
		}


//		int scores_size = 0;
//		for(std::map<int, std::map<int, std::map<int, mScore>> >::iterator it1 = scores.begin(); it1 != scores.end(); it1++)
//		{
//			for(std::map<int, std::map<int, mScore>>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++)
//			{
//				for(std::map<int, mScore>::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++)
//				{
//					scores_size++;
//				}
//			}
//		}
		//std::cerr << "\t\tdiagonalI = " << diagonalI << " => scores_size: " << scores_size << "\n" << std::flush;

		if((diagonalI - lastMaximumIncrease_at_diagonalI) > maximum_steps_nonIncrease)
		{
			break;
		}

		std::map<int, std::map<int, std::map<int, mScore_alternatives > > >  thisDiagonal;
		std::map<int, std::map<int, std::map<int, mScore_backtrace_alternatives > > > thisDiagonal_backtrace;

		if(verbose)
			std::cout << "\t\tfrom m-2 diagonal" << "\n" << std::flush;

		// extend from m-2 diagonal
		for(int m2I = 0; m2I < (int)m2_diagonal.size(); m2I++)
		{
			std::vector<int>& previous_coordinates = m2_diagonal.at(m2I);

			int previous_levelI = previous_coordinates.at(0);
			int previous_seqI = previous_coordinates.at(1);
			int previous_stateI = previous_coordinates.at(2);

			int next_levelI = previous_coordinates.at(0) + (directionPositive ? 1 : -1);
			int next_seqI = previous_coordinates.at(1) + (directionPositive ? 1 : -1);

			if((next_levelI > max_levelI) || (next_seqI > max_seqI))
				continue;

			if((next_levelI < min_levelI) || (next_seqI < min_seqI))
				continue;

			std::string sequenceEmission = (directionPositive ? sequence.substr(previous_seqI, 1) : sequence.substr(previous_seqI-1, 1));

			std::vector<std::pair<int, Edge*> > nextZs = (directionPositive ? _graph_get_next_z_values_and_edges(previous_levelI, previous_stateI) : _graph_get_previous_z_values_and_edges(previous_levelI, previous_stateI));
			assert(nextZs.size() > 0);

			for(unsigned int zI = 0; zI < nextZs.size(); zI++)
			{
				std::pair<int, Edge*>& thisZjump = nextZs.at(zI);
				int next_stateI = thisZjump.first;

				std::string edgeEmission = thisZjump.second->getEmission();
				assert(edgeEmission.size() == 1);

				double score_MatchMismatch = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + ((edgeEmission == sequenceEmission) ? S_match : S_mismatch);

				backtraceStep_affine backtrack_MatchMismatch;
				backtrack_MatchMismatch.x = previous_levelI;
				backtrack_MatchMismatch.y = previous_seqI;
				backtrack_MatchMismatch.z = previous_stateI;
				backtrack_MatchMismatch.usedEdge = thisZjump.second;
				backtrack_MatchMismatch.sourceMatrix = 0;

				thisDiagonal[next_levelI][next_seqI][next_stateI].D.push_back(score_MatchMismatch);
				thisDiagonal_backtrace[next_levelI][next_seqI][next_stateI].D.push_back(backtrack_MatchMismatch);
			}
		}

		if(verbose)
			std::cout << "\t\tfrom m-1 diagonal" << "\n" << std::flush;

		// extend from m-1 diagonal
		for(int m1I = 0; m1I < (int)m1_diagonal.size(); m1I++)
		{
			std::vector<int>& previous_coordinates = m1_diagonal.at(m1I);

			int previous_levelI = previous_coordinates.at(0);
			int previous_seqI = previous_coordinates.at(1);
			int previous_stateI = previous_coordinates.at(2);

			// gap in graph
			int gapInGraph_next_levelI = previous_coordinates.at(0);
			int gapInGraph_next_seqI = previous_coordinates.at(1) + (directionPositive ? 1 : -1);
			if	(
					(directionPositive && (gapInGraph_next_levelI <= max_levelI) && (gapInGraph_next_seqI <= max_seqI)) ||
					((! directionPositive) && (gapInGraph_next_levelI >= min_levelI) && (gapInGraph_next_seqI >= min_seqI))
				)
			{
				double score_gapInGraph_open = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_openGap + S_extendGap;
				if((! directionPositive) && (gapInGraph_next_levelI == 0))
				{
					// score_gapInGraph_open = minusInfinity;
				}
				int gapInGraph_next_stateI = previous_stateI;

				backtraceStep_affine backtrack_gapInGraph_open;
				backtrack_gapInGraph_open.x = previous_levelI;
				backtrack_gapInGraph_open.y = previous_seqI;
				backtrack_gapInGraph_open.z = previous_stateI;
				backtrack_gapInGraph_open.usedEdge = 0;
				backtrack_gapInGraph_open.sourceMatrix = 0;

				thisDiagonal[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(score_gapInGraph_open);
				thisDiagonal_backtrace[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(backtrack_gapInGraph_open);

				double score_gapInGraph_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).GraphGap + S_extendGap;
				if((! directionPositive) && (gapInGraph_next_levelI == 0))
				{
					// score_gapInGraph_extend = minusInfinity;
				}

				backtraceStep_affine backtrack_gapInGraph_extend;
				backtrack_gapInGraph_extend.x = previous_levelI;
				backtrack_gapInGraph_extend.y = previous_seqI;
				backtrack_gapInGraph_extend.z = previous_stateI;
				backtrack_gapInGraph_extend.usedEdge = 0;
				backtrack_gapInGraph_extend.sourceMatrix = 1;

				thisDiagonal[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(score_gapInGraph_extend);
				thisDiagonal_backtrace[gapInGraph_next_levelI][gapInGraph_next_seqI][gapInGraph_next_stateI].GraphGap.push_back(backtrack_gapInGraph_extend);
			}

			// gap in sequence
			int gapInSequence_next_levelI = previous_coordinates.at(0) + (directionPositive ? 1 : -1);
			int gapInSequence_next_seqI = previous_coordinates.at(1);

			if	(
					(directionPositive && (gapInSequence_next_levelI <= max_levelI) && (gapInSequence_next_seqI <= max_seqI)) ||
					((! directionPositive) && (gapInSequence_next_levelI >= min_levelI) && (gapInSequence_next_seqI >= min_seqI))
				)
			{
				std::vector<std::pair<int, Edge*> > nextZs = (directionPositive ? _graph_get_next_z_values_and_edges(previous_levelI, previous_stateI) : _graph_get_previous_z_values_and_edges(previous_levelI, previous_stateI));
				assert(nextZs.size() > 0);

				for(unsigned int zI = 0; zI < nextZs.size(); zI++)
				{
					std::pair<int, Edge*>& thisZjump = nextZs.at(zI);
					int next_stateI = thisZjump.first;

					std::string edgeEmission = thisZjump.second->getEmission();
					assert(edgeEmission.size() == 1);

					// open sequence gap

					double score_gapInSequence_open = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_openGap + S_extendGap;
					if(edgeEmission == "_")
					{
						score_gapInSequence_open = minusInfinity;
					}
					if((! directionPositive) && (gapInSequence_next_seqI == 0))
					{
						// score_gapInSequence_open = minusInfinity;
					}

					backtraceStep_affine backtrack_gapInSequence_open;
					backtrack_gapInSequence_open.x = previous_levelI;
					backtrack_gapInSequence_open.y = previous_seqI;
					backtrack_gapInSequence_open.z = previous_stateI;
					backtrack_gapInSequence_open.usedEdge = thisZjump.second;
					backtrack_gapInSequence_open.sourceMatrix = 0;

					thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(score_gapInSequence_open);
					thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(backtrack_gapInSequence_open);

					// extend sequence gap

					double score_gapInSequence_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap + S_extendGap;
//					if(edgeEmission == "_")
//					{
//						score_gapInSequence_extend = minusInfinity;
//					}
					if(edgeEmission == "_")
					{
						if(scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap == minusInfinity)
						{
							score_gapInSequence_extend = minusInfinity;
						}
						else
						{
							score_gapInSequence_extend = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).SequenceGap +  S_graphGap;
						}
					}
					if((! directionPositive) && (gapInSequence_next_seqI == 0))
					{
						// score_gapInSequence_extend = minusInfinity;
					}

					backtraceStep_affine backtrack_gapInSequence_extend;
					backtrack_gapInSequence_extend.x = previous_levelI;
					backtrack_gapInSequence_extend.y = previous_seqI;
					backtrack_gapInSequence_extend.z = previous_stateI;
					backtrack_gapInSequence_extend.usedEdge = thisZjump.second;
					backtrack_gapInSequence_extend.sourceMatrix = 2;

					thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(score_gapInSequence_extend);
					thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].SequenceGap.push_back(backtrack_gapInSequence_extend);

					// non-affine sequence gap
					if(edgeEmission == "_")
					{
						double score_gapInSequence_nonAffine = scores.at(previous_levelI).at(previous_seqI).at(previous_stateI).D + S_graphGap;

						backtraceStep_affine backtrack_gapInSequence_nonAffine;
						backtrack_gapInSequence_nonAffine.x = previous_levelI;
						backtrack_gapInSequence_nonAffine.y = previous_seqI;
						backtrack_gapInSequence_nonAffine.z = previous_stateI;
						backtrack_gapInSequence_nonAffine.usedEdge = thisZjump.second;
						backtrack_gapInSequence_nonAffine.sourceMatrix = 0;

						thisDiagonal[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].D.push_back(score_gapInSequence_nonAffine);
						thisDiagonal_backtrace[gapInSequence_next_levelI][gapInSequence_next_seqI][next_stateI].D.push_back(backtrack_gapInSequence_nonAffine);
					}
				}
			}
		}

		if(verbose)
			std::cout << "\t\tmaximal" << "\n" << std::flush;

		// call maxima for this diagonal
		std::vector<std::vector<int> > m_thisDiagonal;
		for(std::map<int, std::map<int, std::map<int, mScore_alternatives > > >::iterator diagIt = thisDiagonal.begin(); diagIt != thisDiagonal.end(); diagIt++)
		{
			int levelI = diagIt->first;
			for(std::map<int, std::map<int, mScore_alternatives > >::iterator diagIt2 = diagIt->second.begin(); diagIt2 != diagIt->second.end(); diagIt2++)
			{
				int seqI = diagIt2->first;
				for(std::map<int, mScore_alternatives >::iterator diagIt3 = diagIt2->second.begin(); diagIt3 != diagIt2->second.end(); diagIt3++)
				{
					int stateI = diagIt3->first;

					// call maxima for GraphGap
					std::vector<double>& scores_GraphGap = diagIt2->second.at(stateI).GraphGap;
					std::vector<backtraceStep_affine>& scores_GraphGap_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).GraphGap;
					double selectedScore_GraphGap;
					backtraceStep_affine selectedStep_GraphGap;

					if(scores_GraphGap.size() > 0)
					{
						std::pair<double, unsigned int> maxScore_GraphGap = Utilities::findVectorMaxP_nonCritical(scores_GraphGap, &(rng_seeds.at(omp_get_thread_num())));
						selectedScore_GraphGap = maxScore_GraphGap.first;
						selectedStep_GraphGap = scores_GraphGap_bt.at(maxScore_GraphGap.second);
					}
					else
					{
						selectedScore_GraphGap = minusInfinity;
					}

					// call maxima for SequenceGap
					std::vector<double>& scores_SequenceGap = diagIt2->second.at(stateI).SequenceGap;
					std::vector<backtraceStep_affine>& scores_SequenceGap_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).SequenceGap;

					double selectedScore_SequenceGap;
					backtraceStep_affine selectedStep_SequenceGap;
					if(scores_SequenceGap.size() > 0)
					{
						std::pair<double, unsigned int> maxScore_SequenceGap = Utilities::findVectorMaxP_nonCritical(scores_SequenceGap, &(rng_seeds.at(omp_get_thread_num())));
						selectedScore_SequenceGap = maxScore_SequenceGap.first;
						selectedStep_SequenceGap = scores_SequenceGap_bt.at(maxScore_SequenceGap.second);
					}
					else
					{
						selectedScore_SequenceGap = minusInfinity;
					}

					// two additional steps for D, jumping from the two gap matrices

					double score_intoD_fromGraphGap = selectedScore_GraphGap;
					backtraceStep_affine step_intoD_fromGraphGap;
					step_intoD_fromGraphGap.x = levelI;
					step_intoD_fromGraphGap.y = seqI;
					step_intoD_fromGraphGap.z = stateI;
					step_intoD_fromGraphGap.sourceMatrix = 1;
					step_intoD_fromGraphGap.usedEdge = 0;
					thisDiagonal.at(levelI).at(seqI).at(stateI).D.push_back(score_intoD_fromGraphGap);
					thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D.push_back(step_intoD_fromGraphGap);

					double score_intoD_fromSequenceGap = selectedScore_SequenceGap;
					backtraceStep_affine step_intoD_fromSequenceGap;
					step_intoD_fromSequenceGap.x = levelI;
					step_intoD_fromSequenceGap.y = seqI;
					step_intoD_fromSequenceGap.z = stateI;
					step_intoD_fromSequenceGap.sourceMatrix = 2;
					step_intoD_fromSequenceGap.usedEdge = 0;
					thisDiagonal.at(levelI).at(seqI).at(stateI).D.push_back(score_intoD_fromSequenceGap);
					thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D.push_back(step_intoD_fromSequenceGap);

					// find final maximum for D
					std::vector<double>& scores_D = diagIt2->second.at(stateI).D;
					std::vector<backtraceStep_affine>& scores_D_bt = thisDiagonal_backtrace.at(levelI).at(seqI).at(stateI).D;
					assert(scores_D.size() == scores_D_bt.size());

					std::pair<double, unsigned int> maxScore_D = Utilities::findVectorMaxP_nonCritical(scores_D, &(rng_seeds.at(omp_get_thread_num())));

					bool blockOutCell = false;
					if(directionPositive)
					{
						std::set<NWEdge*> emanatingNWEdges = blockedPathsTable->getEdgesEmanatingFrom(levelI, seqI, stateI);
						if((emanatingNWEdges.size() > 0) && (! returnGlobalScore))
						{
							for(std::set<NWEdge*>::iterator edgeIt = emanatingNWEdges.begin(); edgeIt != emanatingNWEdges.end(); edgeIt++)
							{
								NWEdge* e = *edgeIt;
								NWPath* correspondingPath = e->path;
								if(hit_NW_paths.count(correspondingPath) == 0)
								{
									std::pair<double, std::vector<int> > pathEntry;
									pathEntry.first = maxScore_D.first;
									pathEntry.second.push_back(levelI);
									pathEntry.second.push_back(seqI);
									pathEntry.second.push_back(stateI);
									hit_NW_paths[correspondingPath] = pathEntry;
								}
							}
							blockOutCell = true;
						}
					}
					else
					{
						std::set<NWEdge*> incomingNWEdges = blockedPathsTable->getEdgesGoingInto(levelI, seqI, stateI);
						if((incomingNWEdges.size() > 0) && (! returnGlobalScore))
						{
							for(std::set<NWEdge*>::iterator edgeIt = incomingNWEdges.begin(); edgeIt != incomingNWEdges.end(); edgeIt++)
							{
								NWEdge* e = *edgeIt;
								NWPath* correspondingPath = e->path;
								if(hit_NW_paths.count(correspondingPath) == 0)
								{
									std::pair<double, std::vector<int> > pathEntry;
									pathEntry.first = maxScore_D.first;
									pathEntry.second.push_back(levelI);
									pathEntry.second.push_back(seqI);
									pathEntry.second.push_back(stateI);
									hit_NW_paths[correspondingPath] = pathEntry;
								}
							}
							blockOutCell = true;
						}
					}

					if(blockOutCell)
					{
						double scoreThatWillBeDeleted = maxScore_D.first;
						std::vector<int> thisCoordinates;
						thisCoordinates.push_back(levelI);
						thisCoordinates.push_back(seqI);
						thisCoordinates.push_back(stateI);
						if(scoreThatWillBeDeleted == currentMaximum)
						{
							currentMaxima_coordinates.push_back(thisCoordinates);
							lastMaximumIncrease_at_diagonalI = diagonalI;
						}
						else if(scoreThatWillBeDeleted > currentMaximum)
						{
							currentMaximum = scoreThatWillBeDeleted;
							currentMaxima_coordinates.clear();
							currentMaxima_coordinates.push_back(thisCoordinates);
							lastMaximumIncrease_at_diagonalI = diagonalI;
						}

						scores[levelI][seqI][stateI].D = minusInfinity;
						scores_backtrace[levelI][seqI][stateI].D = scores_D_bt.at(maxScore_D.second);

						assert(scores_backtrace[levelI][seqI][stateI].D.x != -1);

						scores[levelI][seqI][stateI].GraphGap = minusInfinity;
						scores_backtrace[levelI][seqI][stateI].GraphGap = selectedStep_GraphGap;

						scores[levelI][seqI][stateI].SequenceGap = minusInfinity;
						scores_backtrace[levelI][seqI][stateI].SequenceGap = selectedStep_SequenceGap;
					}
					else
					{
						if(maxScore_D.first >= diagonal_stop_threshold)
						{
							scores[levelI][seqI][stateI].D = maxScore_D.first;
							scores_backtrace[levelI][seqI][stateI].D = scores_D_bt.at(maxScore_D.second);

							assert(scores_backtrace[levelI][seqI][stateI].D.x != -1);

							scores[levelI][seqI][stateI].GraphGap = selectedScore_GraphGap;
							scores_backtrace[levelI][seqI][stateI].GraphGap = selectedStep_GraphGap;

							scores[levelI][seqI][stateI].SequenceGap = selectedScore_SequenceGap;
							scores_backtrace[levelI][seqI][stateI].SequenceGap = selectedStep_SequenceGap;

							if(preferSequenceCompleAlignments)
							{
								std::string id_levelI_stateI = Utilities::ItoStr(levelI) + "/" + Utilities::ItoStr(stateI);
								if(directionPositive)
								{
									if(seqI == max_seqI)
									{
										achieved_complete_sequence_alignments.insert(id_levelI_stateI);
									}
								}
								else
								{
									if(seqI == min_seqI)
									{
										achieved_complete_sequence_alignments.insert(id_levelI_stateI);
									}
								}
							}

							std::vector<int> thisCoordinates;
							thisCoordinates.push_back(levelI);
							thisCoordinates.push_back(seqI);
							thisCoordinates.push_back(stateI);
							m_thisDiagonal.push_back(thisCoordinates);

							backtraceStep_affine oneRealStepBackwards = scores_backtrace[levelI][seqI][stateI].D;
							while((oneRealStepBackwards.x == levelI) && (oneRealStepBackwards.y == seqI))
							{
								assert(oneRealStepBackwards.sourceMatrix != 0);
								if(oneRealStepBackwards.sourceMatrix == 1)
								{
									oneRealStepBackwards = scores_backtrace[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].GraphGap;
								}
								else if(oneRealStepBackwards.sourceMatrix == 2)
								{
									oneRealStepBackwards = scores_backtrace[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].SequenceGap;
								}
								else
								{
									assert(1 == 0);
								}
							}
							int previousScore;
							if(oneRealStepBackwards.sourceMatrix == 0)
							{
								previousScore = scores[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].D;
							}
							else if(oneRealStepBackwards.sourceMatrix == 1)
							{
								previousScore = scores[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].GraphGap;
							}
							else if(oneRealStepBackwards.sourceMatrix == 2)
							{
								previousScore = scores[oneRealStepBackwards.x][oneRealStepBackwards.y][oneRealStepBackwards.z].SequenceGap;
							}
							else
							{
								assert( 1 == 0 );
							}
							int scoreDifference = maxScore_D.first - previousScore;

							if(maxScore_D.first == currentMaximum)
							{
								if(scoreDifference != 0)
								{
									currentMaxima_coordinates.push_back(thisCoordinates);
									lastMaximumIncrease_at_diagonalI = diagonalI;
								}
							}
							else if(maxScore_D.first > currentMaximum)
							{
								currentMaximum = maxScore_D.first;
								currentMaxima_coordinates.clear();
								currentMaxima_coordinates.push_back(thisCoordinates);
								lastMaximumIncrease_at_diagonalI = diagonalI;
							}
						}
						else
						{
							// assert(1 == 0);
						}
					}
				}
			}
		}

		if(verbose)
			std::cout << "\t\tfiltering" << "\n" << std::flush;

		std::vector<std::vector<int> > m_thisDiagonal_filtered;
		if(m_thisDiagonal.size() > 0)
		{
			double max;
			for(unsigned int i = 0; i < m_thisDiagonal.size(); i++)
			{
				std::vector<int> coordinates = m_thisDiagonal.at(i);
				double S = scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D;
				if((i == 0) || (max < S))
				{
					max = S;
				}
			}

			for(unsigned int i = 0; i < m_thisDiagonal.size(); i++)
			{
				std::vector<int> coordinates = m_thisDiagonal.at(i);
				double S = scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D;
				assert(S <= max);
				if((max - S) <= threshold_for_filtering)
				{
					m_thisDiagonal_filtered.push_back(coordinates);
				}
			}

			m_thisDiagonal = m_thisDiagonal_filtered;
		}

		m2_diagonal = m1_diagonal;
		m1_diagonal = m_thisDiagonal;
	}

	std::vector<localExtension_pathDescription> forReturn;
	auto backtraceFrom = [&](int start_x, int start_y, int start_z, double StartScore) {
		int backtrace_x = start_x;
		int backtrace_y = start_y;
		int backtrace_z = start_z;
		int backtrace_matrix = 0;

		std::string reconstructed_graph;
		std::string reconstructed_sequence;
		std::vector<int> reconstructed_graph_levels;
		std::vector<std::vector<int>> edge_coordinates;
		std::vector<Edge*> utilizedEdges;

		std::vector<int> startCoordinates;
		startCoordinates.push_back(start_x);
		startCoordinates.push_back(start_y);
		startCoordinates.push_back(start_z);
		edge_coordinates.push_back(startCoordinates);

		if(verbose)
		{
			std::cout << "\tbacktraceFrom() called.\n";
			std::cout << "\t\t" << "start_x" << ": " << start_x << "\n";
			std::cout << "\t\t" << "start_y" << ": " << start_y << "\n";
			std::cout << "\t\t" << "start_z" << ": " << start_z << "\n";
			std::cout << std::flush;
		}
		while((backtrace_x != startLevel_graph) || (backtrace_y != start_sequence))
		{
			if(verbose)
			{
				std::cout << "\t\tbacktrace_x: " << backtrace_x << "\n";
				std::cout << "\t\tbacktrace_y: " << backtrace_y << "\n";
				std::cout << "\t\tbacktrace_z: " << backtrace_z << "\n";
				std::cout << "\t\tselected matrix: " << backtrace_matrix << "\n" << std::flush;
			}

			if(directionPositive)
			{
				assert((backtrace_x >= startLevel_graph) && (backtrace_x <= ((int)levels - 1)));
				assert((backtrace_y >= start_sequence) && (backtrace_y <= (int)sequenceLength));
				assert((backtrace_z >= 0) && (backtrace_z <= (int)g->NodesPerLevel.at(backtrace_x).size()));
				assert((backtrace_matrix >= 0) && (backtrace_matrix <= 2));
			}
			else
			{
				assert((backtrace_x <= startLevel_graph) && (backtrace_x >= 0));
				assert((backtrace_y <= start_sequence) && (backtrace_y >= 0));
				assert((backtrace_z >= 0) && (backtrace_z <= (int)g->NodesPerLevel.at(backtrace_x).size()));
				assert((backtrace_matrix >= 0) && (backtrace_matrix <= 2));
			}

			double Score;
			backtraceStep_affine step;
			if(backtrace_matrix == 0)
			{
				Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).D;
				step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).D;
			}
			else if(backtrace_matrix == 1)
			{
				Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).GraphGap;
				step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).GraphGap;
			}
			else if(backtrace_matrix == 2)
			{
				Score = scores.at(backtrace_x).at(backtrace_y).at(backtrace_z).SequenceGap;
				step = scores_backtrace.at(backtrace_x).at(backtrace_y).at(backtrace_z).SequenceGap;
			}

			if((backtrace_x == start_x) && (backtrace_y == start_y) && (backtrace_z == start_z))
			{
				Score = StartScore;
			}

			if(verbose)
			{
				std::cout << "\t\tScore: " << Score << "\n\n" << std::flush;
			}

			std::string edgeEmission;
			if(step.usedEdge != 0)
			{
				edgeEmission = step.usedEdge->getEmission();
				assert(edgeEmission.size() == 1);
			}
			std::string sequenceEmission;
			if((backtrace_y >= 1) && (directionPositive))
			{
				sequenceEmission = sequence.substr(backtrace_y - 1, 1);
			}
			if((backtrace_y < max_seqI) && (! directionPositive))
			{
				sequenceEmission = sequence.substr(backtrace_y, 1);
			}

			int next_x = step.x;
			int next_y = step.y;
			int next_z = step.z;
			int next_matrix = step.sourceMatrix;

			bool dontAddCoordinates = false; // if we jump from one matrix to the other without changing coordinates, we don't store the coordinates...
			if(directionPositive)
			{
				if((next_x == (backtrace_x - 1)) && (next_y == (backtrace_y - 1)))
				{
					// match or mismatch
					assert(step.usedEdge != 0);
					reconstructed_graph.append(edgeEmission);
					reconstructed_graph_levels.push_back(backtrace_x - 1);
					reconstructed_sequence.append(sequenceEmission);
					utilizedEdges.push_back(step.usedEdge);
				}
				else if((next_x == backtrace_x) && (next_y == (backtrace_y - 1)))
				{
					// gap in graph
					reconstructed_graph.append("_");
					reconstructed_graph_levels.push_back(-1);
					reconstructed_sequence.append(sequenceEmission);
					utilizedEdges.push_back(0);

				}
				else if((next_x == (backtrace_x - 1)) && (next_y == backtrace_y))
				{
					// gap in sequence
					assert(step.usedEdge != 0);
					reconstructed_graph.append(edgeEmission);
					reconstructed_graph_levels.push_back(backtrace_x - 1);
					reconstructed_sequence.append("_");
					utilizedEdges.push_back(step.usedEdge);
				}
				else
				{
					dontAddCoordinates = true;
					assert((backtrace_x == next_x) && (backtrace_y == next_y) && (backtrace_z == next_z) && (next_matrix != backtrace_matrix));
				}
			}
			else
			{
				if((next_x == (backtrace_x + 1)) && (next_y == (backtrace_y + 1)))
				{
					// match or mismatch
					assert(step.usedEdge != 0);
					reconstructed_graph.append(edgeEmission);
					reconstructed_graph_levels.push_back(backtrace_x);
					reconstructed_sequence.append(sequenceEmission);
					utilizedEdges.push_back(step.usedEdge);
				}
				else if((next_x == backtrace_x) && (next_y == (backtrace_y + 1)))
				{
					// gap in graph
					reconstructed_graph.append("_");
					reconstructed_graph_levels.push_back(-1);
					reconstructed_sequence.append(sequenceEmission);
					utilizedEdges.push_back(0);
				}
				else if((next_x == (backtrace_x + 1)) && (next_y == backtrace_y))
				{
					// gap in sequence
					assert(step.usedEdge != 0);
					reconstructed_graph.append(edgeEmission);
					reconstructed_graph_levels.push_back(backtrace_x);
					reconstructed_sequence.append("_");
					utilizedEdges.push_back(step.usedEdge);
				}
				else
				{
					dontAddCoordinates = true;
					assert((backtrace_x == next_x) && (backtrace_y == next_y) && (backtrace_z == next_z) && (next_matrix != backtrace_matrix));
				}
			}

			if(verbose)
			{
				std::cout << "\t\t\t" << "Emission: " << edgeEmission << "\n";
				std::cout << "\t\t\t" << "next_x: " << next_x << "\n";
				std::cout << "\t\t\t" << "next_y: " << next_y << "\n";
				std::cout << "\t\t\t" << "next_z: " << next_z << "\n";
				std::cout << "\t\t\t" << "next_matrix: " << next_matrix << "\n\n" << std::flush;
			}

			backtrace_x = next_x;
			backtrace_y = next_y;
			backtrace_z = next_z;
			backtrace_matrix = next_matrix;

			if( ! dontAddCoordinates)
			{
				std::vector<int> nextCoordinates;
				nextCoordinates.push_back(next_x);
				nextCoordinates.push_back(next_y);
				nextCoordinates.push_back(next_z);
				edge_coordinates.push_back(nextCoordinates);
			}
		}

		if(directionPositive)
		{
			std::reverse(reconstructed_graph.begin(), reconstructed_graph.end());
			std::reverse(reconstructed_graph_levels.begin(), reconstructed_graph_levels.end());
			std::reverse(reconstructed_sequence.begin(), reconstructed_sequence.end());
			std::reverse(utilizedEdges.begin(), utilizedEdges.end());
			std::reverse(edge_coordinates.begin(), edge_coordinates.end());
		}

		assert(reconstructed_graph.size() == reconstructed_sequence.size());
		assert(reconstructed_graph_levels.size() == reconstructed_graph.size());

		localExtension_pathDescription pathReturn;
		pathReturn.Score = StartScore;
		pathReturn.usedEdges = utilizedEdges;
		pathReturn.coordinates = edge_coordinates;
		pathReturn.alignedSequence = reconstructed_sequence;
		pathReturn.alignedGraph = reconstructed_graph;
		pathReturn.alignedGraph_levels = reconstructed_graph_levels;

		forReturn.push_back(pathReturn);
	};

	if(verbose)
	{
			/*
		for(unsigned int x = 231412; x <= 231418; x++)
		{
			for(unsigned int y = 0; y <= 5; y++)
			{
				std::cerr << "x = " << x << ", y = " << y << ": ";
				if((scores.count(x) == 0) || (scores.at(x).count(y) == 0) || (scores.at(x).at(y).size() == 0))
				{
					std::cerr << "Undefined.";
				}
				else
				{
					for(unsigned int z = 0; z < scores.at(x).at(y).size(); z++)
					{
						std::cerr << "\n\t z = " << z << ": " << scores.at(x).at(y).at(z).D;
					}
				}
				std::cerr << "\n\n" << std::flush;
			}
		}
		*/
	}

	if(! returnGlobalScore)
	{
		bool have_backtrace_sequenceComplete = false;

		int sequenceCompleteBacktrace_levelI;
		int sequenceCompleteBacktrace_seqI;
		int sequenceCompleteBacktrace_stateI;
		double sequenceCompleteBacktrace_score;


		if(preferSequenceCompleAlignments)
		{

			int coordinate_seqI = (directionPositive ? max_seqI : min_seqI);

			int _dbg_local_chainI = 0;
			if(0 && directionPositive && (_dbg_local_chainI == 0))
			{
				std::cout << "fullNeedleman_affine_diagonal_extension(..) called, direction " << directionPositive << ".\n" << std::flush;
				std::cout << "\t" << "_dbg_local_chainI" << ": " << _dbg_local_chainI << "\n";
				std::cout << "\t" << "start_sequence" << ": " << start_sequence << "\n";
				std::cout << "\t" << "maxPosition_sequence" << ": " << maxPosition_sequence << "\n";
				std::cout << "\t" << "startLevel_graph" << ": " << startLevel_graph << "\n";
				std::cout << "\t" << "maxLevel_graph" << ": " << maxLevel_graph << "\n";
				std::cout << "\t" << "achieved_complete_sequence_alignments.size()" << ": " << achieved_complete_sequence_alignments.size() << "\n";

				for(std::set<std::string>::iterator coordinateIt = achieved_complete_sequence_alignments.begin(); coordinateIt != achieved_complete_sequence_alignments.end(); coordinateIt++)
				{
					std::string coordinates = *coordinateIt;
					std::vector<std::string> coordinate_components = Utilities::split(coordinates, "/");
					assert(coordinate_components.size() == 2);
					int coordinate_levelI = Utilities::StrtoI(coordinate_components.at(0));
					int coordinate_stateI = Utilities::StrtoI(coordinate_components.at(1));

					double S = scores.at(coordinate_levelI).at(coordinate_seqI).at(coordinate_stateI).D;

					std::cout << "\t\t\tAlternative:\n";
					std::cout << "\t\t\t\t" << "coordinate_seqI" << ": " << coordinate_seqI << "\n";
					std::cout << "\t\t\t\t" << "coordinate_levelI" << ": " << coordinate_levelI << "\n";
					std::cout << "\t\t\t\t" << "coordinate_stateI" << ": " << coordinate_stateI << "\n";
					std::cout << "\t\t\t\t" << "S" << ": " << S << "\n";
				}

				std::cout << std::flush;
			}

			std::vector<std::string> sequenceComplete_backtraceCoordinates_maxScores;
			double maxScore;


			for(std::set<std::string>::iterator coordinateIt = achieved_complete_sequence_alignments.begin(); coordinateIt != achieved_complete_sequence_alignments.end(); coordinateIt++)
			{
				std::string coordinates = *coordinateIt;
				std::vector<std::string> coordinate_components = Utilities::split(coordinates, "/");
				assert(coordinate_components.size() == 2);
				int coordinate_levelI = Utilities::StrtoI(coordinate_components.at(0));
				int coordinate_stateI = Utilities::StrtoI(coordinate_components.at(1));

				double S = scores.at(coordinate_levelI).at(coordinate_seqI).at(coordinate_stateI).D;

				if((coordinateIt == achieved_complete_sequence_alignments.begin()) || (S > maxScore))
				{
					sequenceComplete_backtraceCoordinates_maxScores.clear();
					sequenceComplete_backtraceCoordinates_maxScores.push_back(coordinates);
					maxScore = S;
				}
				else if(S == maxScore)
				{
					sequenceComplete_backtraceCoordinates_maxScores.push_back(coordinates);
				}
			}

			// std::cout << "\t\t" << "sequenceComplete_backtraceCoordinates_maxScores.size()" << ": " << sequenceComplete_backtraceCoordinates_maxScores.size() << "\n" << std::flush;

			if(sequenceComplete_backtraceCoordinates_maxScores.size() > 0)
			{
				have_backtrace_sequenceComplete = true;

				int selectedIndex = Utilities::randomNumber_nonCritical(sequenceComplete_backtraceCoordinates_maxScores.size() - 1, &(rng_seeds.at(omp_get_thread_num())));
				assert(selectedIndex >= 0);
				assert(selectedIndex < (int)sequenceComplete_backtraceCoordinates_maxScores.size());
				std::string coordinates = sequenceComplete_backtraceCoordinates_maxScores.at(selectedIndex);

				std::vector<std::string> coordinate_components = Utilities::split(coordinates, "/");
				assert(coordinate_components.size() == 2);

				sequenceCompleteBacktrace_levelI = Utilities::StrtoI(coordinate_components.at(0));
				sequenceCompleteBacktrace_seqI = coordinate_seqI;
				sequenceCompleteBacktrace_stateI = Utilities::StrtoI(coordinate_components.at(1));
				sequenceCompleteBacktrace_score = maxScore;

			}

		}
		if(preferSequenceCompleAlignments && have_backtrace_sequenceComplete)
		{
			backtraceFrom(sequenceCompleteBacktrace_levelI, sequenceCompleteBacktrace_seqI, sequenceCompleteBacktrace_stateI, sequenceCompleteBacktrace_score);
		}
		else
		{
			if(currentMaximum > 0)
			{
				if(verbose)
				{
					std::cout << "\tMaximum " << currentMaximum << ", achieved at " << currentMaxima_coordinates.size() << " positions!\n" << std::flush;
				}

				for(int maximumI = 0; maximumI < (int)currentMaxima_coordinates.size(); maximumI++)
				{
					std::vector<int> coordinates = currentMaxima_coordinates.at(maximumI);
					if(scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D != minusInfinity)
					{
						if(verbose)
							std::cout << " - Start maximum backtrace from " <<coordinates.at(0) << ", " <<  coordinates.at(1) << ", " << coordinates.at(2) << "\n" << std::flush;

						backtraceFrom(coordinates.at(0), coordinates.at(1), coordinates.at(2), scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D);
					}
				}
			}

			if(verbose)
				std::cout << "Have hit " << hit_NW_paths.size() << " NW paths, backtrace independent of achieved value!\n" << std::flush;

			for(std::map<NWPath*, std::pair<double, std::vector<int> > >::iterator hitPathsIt = hit_NW_paths.begin(); hitPathsIt != hit_NW_paths.end(); hitPathsIt++)
			{
				// NWPath* path = hitPathsIt->first;
				double score = hitPathsIt->second.first;
				std::vector<int> coordinates = hitPathsIt->second.second;
				assert(scores.at(coordinates.at(0)).at(coordinates.at(1)).at(coordinates.at(2)).D == minusInfinity);

				if(verbose)
					std::cout << " - Start NW path backtrace from " << coordinates.at(0) << ", " <<  coordinates.at(1) << ", " << coordinates.at(2) << "\n" << std::flush;

				backtraceFrom(coordinates.at(0), coordinates.at(1), coordinates.at(2), score);
			}
		}
	}
	else
	{

		if(directionPositive)
		{
			std::map<int, double> finalState_scores;
			for(std::map<int, mScore>::iterator altIt = scores.at(levels-1).at(sequenceLength).begin(); altIt != scores.at(levels-1).at(sequenceLength).end(); altIt++)
			{
				finalState_scores[altIt->first] = altIt->second.D;
			}
			std::pair<double, int> maxFinalState = Utilities::findIntMapMaxP_nonCritical(finalState_scores, &(rng_seeds.at(omp_get_thread_num())));

			backtraceFrom(levels - 1, sequenceLength, maxFinalState.second, maxFinalState.first);
		}
		else
		{
			std::map<int, double> finalState_scores;
			for(std::map<int, mScore>::iterator altIt = scores.at(0).at(0).begin(); altIt != scores.at(0).at(0).end(); altIt++)
			{
				finalState_scores[altIt->first] = altIt->second.D;

				if(verbose)
					std::cout << "Final state z-value " << altIt->first << ", D value " << scores.at(0).at(0).at(altIt->first).D << ", GraphGap: " <<  scores.at(0).at(0).at(altIt->first).GraphGap << ", SequenceGap: " <<  scores.at(0).at(0).at(altIt->first).SequenceGap << "\n" << std::flush;
			}
			std::pair<double, int> maxFinalState = Utilities::findIntMapMaxP_nonCritical(finalState_scores, &(rng_seeds.at(omp_get_thread_num())));

			backtraceFrom(0, 0, maxFinalState.second, maxFinalState.first);
		}
	}

	std::vector<mapper::reads::verboseSeedChain> vC;
	for(unsigned int i = 0; i < forReturn.size(); i++)
	{
		vC.push_back(forReturn.at(0).toVerboseSeedChain());
	}

	return vC;
	// return forReturn;
}


std::vector<std::pair<int, std::pair<int, Edge*>> > extensionAligner::_graph_get_jumpNext_x_and_z_values_and_edges(int x, int z) const
{
	assert(g != 0);
	std::vector< std::pair<int, std::pair<int, Edge*> > > forReturn;
	unsigned int uZ = (unsigned int)z;
	assert(x >= 0);
	assert(x < (int)nodesPerLevel_ordered.size());
	assert(z >= 0);
	assert(z < (int)nodesPerLevel_ordered.at(x).size());
	Node* thisZ = nodesPerLevel_ordered.at(x).at(uZ);
	if(g->gapEdgePaths_connectedNodes_forwards.count(thisZ))
	{
		for(std::map<Node*, Edge*>::iterator targetNodeIt = g->gapEdgePaths_connectedNodes_forwards.at(thisZ).begin(); targetNodeIt != g->gapEdgePaths_connectedNodes_forwards.at(thisZ).end(); targetNodeIt++)
		{
			Node* toNode = targetNodeIt->first;
			Edge* pseudoEdge = targetNodeIt->second;
			int toNode_level = toNode->level;
			int z_for_toNode = nodesPerLevel_ordered_rev.at(toNode_level).at(toNode);
			forReturn.push_back(make_pair(toNode_level, make_pair(z_for_toNode, pseudoEdge)));
		}
	}

	return forReturn;
}

std::vector<std::pair<int, std::pair<int, Edge*>>> extensionAligner::_graph_get_jumpPrevious_x_and_z_values_and_edges(int x, int z) const
{
	assert(g != 0);
	std::vector<std::pair<int, std::pair<int, Edge*> >> forReturn;
	unsigned int uZ = (unsigned int)z;
	assert(x >= 0);
	assert(x < (int)nodesPerLevel_ordered.size());
	assert(z >= 0);
	assert(z < (int)nodesPerLevel_ordered.at(x).size());
	Node* thisZ = nodesPerLevel_ordered.at(x).at(uZ);
	if(g->gapEdgePaths_connectedNodes_backwards.count(thisZ))
	{
		for(std::map<Node*, Edge*>::iterator targetNodeIt = g->gapEdgePaths_connectedNodes_backwards.at(thisZ).begin(); targetNodeIt != g->gapEdgePaths_connectedNodes_backwards.at(thisZ).end(); targetNodeIt++)
		{
			Node* toNode = targetNodeIt->first;
			Edge* pseudoEdge = targetNodeIt->second;
			int toNode_level = toNode->level;
			int z_for_toNode = nodesPerLevel_ordered_rev.at(toNode_level).at(toNode);
			forReturn.push_back(make_pair(toNode_level, make_pair(z_for_toNode, pseudoEdge)));
		}
	}

	return forReturn;
}

std::vector<std::pair<int, int>> extensionAligner::_graph_get_jumpPrevious_x_and_z_values(int x, int z) const
{
	assert(g != 0);

	std::vector<std::pair<int, std::pair<int, Edge*>>> previous_X_and_Zs = _graph_get_jumpPrevious_x_and_z_values_and_edges(x, z);

	std::vector<std::pair<int, int>> forReturn;
	for(unsigned int i = 0; i < previous_X_and_Zs.size(); i++)
	{
		std::pair<int, int> previous_X_and_Z = make_pair(previous_X_and_Zs.at(i).first, previous_X_and_Zs.at(i).second.first);
		if(std::find(forReturn.begin(), forReturn.end(), previous_X_and_Z) == forReturn.end())
		{
			forReturn.push_back(previous_X_and_Z);
		}
	}

	return forReturn;
}


} /* namespace aligner */
} /* namespace mapper */
