/*
 * VirtualNW.cpp
 *
 *  Created on: 30.06.2013
 *      Author: AlexanderDilthey
 */

#include "VirtualNWUnique.h"

#include <assert.h>
#include <string>
#include <iostream>
#include <algorithm>

#include "../../Utilities.h"

namespace mapper {
namespace aligner {

mapper::reads::verboseSeedChain localExtension_pathDescription::toVerboseSeedChain()
{
	mapper::reads::verboseSeedChain forReturn;

	forReturn.graph_aligned_edges = usedEdges;
	forReturn.graph_aligned_levels = alignedGraph_levels;
	forReturn.graph_aligned = alignedGraph;
	forReturn.sequence_aligned = alignedSequence;
	forReturn.sequence_begin = coordinates.front().at(1);
	forReturn.sequence_end = coordinates.back().at(1) - 1;
	if(!(forReturn.sequence_begin <= forReturn.sequence_end))
	{
		_printExtension();
		forReturn.print();
	}
	assert(forReturn.sequence_begin <= forReturn.sequence_end);
	forReturn.reverse = false;

	return forReturn;
}

NWEdge::NWEdge()
{
	path = 0;
	usedGraphEdge = 0;
	scoreComputed = false;
	endsFree_previousEdgeAffineSequenceGap = false;
}

void VirtualNWTable_Unique::_enumeratePathsAndEdges()
{
	for(std::set<NWPath*>::iterator pathIt = paths.begin(); pathIt != paths.end(); pathIt++)
	{
		NWPath* p = *pathIt;
		if(_ordered_paths.count(p) == 0)
		{
			_ordered_paths[p] = _ordered_paths.size() + 1;
		}

		std::set<NWEdge*> path_edges = p->all_edges;
		std::map<NWEdge*, int> ordered_edges;
		for(std::set<NWEdge*>::iterator edgeIt = path_edges.begin(); edgeIt != path_edges.end(); edgeIt++)
		{
			NWEdge* e = *edgeIt;
			if(ordered_edges.count(e) == 0)
			{
				ordered_edges[e] = ordered_edges.size() + 1;
			}
		}

		_ordered_paths_edges[p] = ordered_edges;
	}
}

void VirtualNWTable_Unique::printSequenceCoverageStats()
{
	vector<int> coveragePerLevel;
	int uncovered_positions = 0;
	assert(Context->sequence != 0);
	for(int y = 0; y <= (int)Context->sequence->length(); y++)
	{
		std::set<NWEdge*> edges_starting_thisY = (index_edges_start_y.count(y) ? index_edges_start_y[y] : std::set<NWEdge*>());;
		int relevantSequenceEdges = 0;
		for(std::set<NWEdge*>::iterator eIt = edges_starting_thisY.begin(); eIt != edges_starting_thisY.end(); eIt++)
		{
			NWEdge* e = *eIt;
			assert(e != 0);
			if(! e->isSequenceGap())
			{
				relevantSequenceEdges++;
			}
		}
		coveragePerLevel.push_back(relevantSequenceEdges);
		if(relevantSequenceEdges == 0)
		{
			uncovered_positions++;
		}
	}

	std::sort(coveragePerLevel.begin(), coveragePerLevel.end());

	int coveragePerLevel_sum = 0;
	std::for_each(coveragePerLevel.begin(), coveragePerLevel.end(), [&](int i){coveragePerLevel_sum += i;});

	assert(coveragePerLevel.size() > 0);
	int coveragePerLevel_median = coveragePerLevel.at(coveragePerLevel.size() / 2);

	std::cout << "printSequenceCoverageStats():\n";
	std::cout << "\tSequence levels: " << coveragePerLevel.size() << "\n";
	std::cout << "\tMean coveragePerLevel: " << ( (double)coveragePerLevel_sum/(double)coveragePerLevel.size() ) << "\n";
	std::cout << "\tMedian: " << coveragePerLevel_median << "\n";
	std::cout << "\tSequence positions uncovered: " << uncovered_positions << "\n" << std::flush;
}

void VirtualNWTable_Unique::print()
{
	int maxX = Context->aligner->g->NodesPerLevel.size() - 1;
	int maxY = Context->sequence->length();

	std::vector<std::vector<std::vector<std::string>> > printMatrix;
	printMatrix.resize(maxX + 1);
	for(int levelI = 0; levelI <= maxX; levelI++)
	{
		int zMax = Context->aligner->g->NodesPerLevel.at(levelI).size();
		printMatrix.at(levelI).resize(maxY+1);
		for(int seqI = 0; seqI <= maxY; seqI++)
		{
			printMatrix.at(levelI).at(seqI).resize(zMax);
		}
	}



	for(std::set<NWPath*>::iterator pathIt = paths.begin(); pathIt != paths.end(); pathIt++)
	{
		NWPath* p = *pathIt;
		int path_index = _ordered_paths.at(p);

		std::set<NWEdge*> path_edges = p->all_edges;
		for(std::set<NWEdge*>::iterator edgeIt = path_edges.begin(); edgeIt != path_edges.end(); edgeIt++)
		{
			NWEdge* e = *edgeIt;
			int edge_index = _ordered_paths_edges.at(p).at(e);
			std::string edgeName = Utilities::ItoStr(path_index)+"|"+Utilities::ItoStr(edge_index);

			std::string edgeBacktrackID;
			std::string backtrackPathID;
			if(e->scoreBacktrack == 0)
			{
				edgeBacktrackID = "ORIGIN";
				backtrackPathID = "ORIGIN";
			}
			else
			{
				edgeBacktrackID = Utilities::ItoStr(_ordered_paths.at(e->scoreBacktrack->path))+"|"+Utilities::ItoStr(_ordered_paths_edges.at(e->scoreBacktrack->path).at(e->scoreBacktrack));
				backtrackPathID = Utilities::ItoStr(_ordered_paths.at(e->scoreBacktrack->path));
			}
//
//			printMatrix.at(e->from_x).at(e->from_y).at(e->from_z) += ( edgeName + "|from:" + edgeBacktrackID + ", ");
//			assert(e->scoreComputed);
//			printMatrix.at(e->to_x).at(e->to_y).at(e->to_z) += ( edgeName + "|score:" + Utilities::DtoStr(e->scoreAfterEdge) + ", ");

			if(p->entry_edges.count(e))
			{
				printMatrix.at(e->from_x).at(e->from_y).at(e->from_z) += ( "ENTRY " + Utilities::ItoStr(path_index) + "; FROM "+backtrackPathID + ". ");
			}
			std::string edgeLabel;
			if(e->usedGraphEdge != 0)
			{
				edgeLabel = e->usedGraphEdge->getEmission();
			}

			std::string exitLabel;
			if(p->exit_edges.count(e))
			{
				exitLabel = ":ExitEdge";
			}
			printMatrix.at(e->to_x).at(e->to_y).at(e->to_z) += (  Utilities::ItoStr(path_index) + ":" + edgeLabel+ ":"+ Utilities::DtoStr(e->scoreAfterEdge) + exitLabel + " ");


		}
	}

	std::cout << "\t";
	for(int levelI = 0; levelI <= maxX; levelI++)
	{
		std::cout << "" << levelI << "\t";
	}
	std::cout << "\n";

	for(int seqI = 0; seqI <= maxY; seqI++)
	{
		std::cout << "seqI=" << seqI << "\t";
		for(int levelI = 0; levelI <= maxX; levelI++)
		{
			bool printCell = false;
			for(int z = 0; z < (int)printMatrix.at(levelI).at(seqI).size(); z++)
			{
				if(printMatrix.at(levelI).at(seqI).at(z).length() > 0)
				{
					printCell = true;
					break;
				}
			}

			if(printCell)
			{
				if(printMatrix.at(levelI).at(seqI).size() > 1)
				{
					for(int z = 0; z < (int)printMatrix.at(levelI).at(seqI).size(); z++)
					{
						//std::cout << "{" << z << ": " << printMatrix.at(levelI).at(seqI).at(z) << "} ";
						std::cout << "{ " << printMatrix.at(levelI).at(seqI).at(z) << " } ";
					}
				}
				else
				{
					std::cout << printMatrix.at(levelI).at(seqI).at(0);
				}
			}
			std::cout << "\t";
		}
		std::cout << "\n" << std::flush;
	}

	std::cout << std::flush;
}

void VirtualNWTable_Unique::addPath(NWPath* p)
{
	assert(paths.count(p) == 0);

	p->table = this;
	for(std::set<NWEdge*>::iterator edgeIt = p->all_edges.begin(); edgeIt != p->all_edges.end(); edgeIt++)
	{
		NWEdge* edge = *edgeIt;
		index_edges_start_x[edge->from_x].insert(edge);
		index_edges_stop_x[edge->to_x].insert(edge);
		index_edges_start_y[edge->from_y].insert(edge);
		index_edges_stop_y[edge->to_y].insert(edge);

		index_edges_from[edge->from_x][edge->from_y][edge->from_z].insert(edge);
		index_edges_to[edge->to_x][edge->to_y][edge->to_z].insert(edge);

		if(edge->usedGraphEdge != 0)
		{
			index_graphEdges_2_NWEdges[edge->usedGraphEdge].insert(edge);
		}
	}

	p->checkConsistency();
	paths.insert(p);
}


void VirtualNWTable_Unique::removePath(NWPath* p)
{
	assert(p->table == this);
	assert(paths.count(p) > 0);

	for(std::set<NWEdge*>::iterator edgeIt = p->all_edges.begin(); edgeIt != p->all_edges.end(); edgeIt++)
	{
		NWEdge* edge = *edgeIt;
		index_edges_start_x[edge->from_x].erase(edge);
		index_edges_stop_x[edge->to_x].erase(edge);
		index_edges_start_y[edge->from_y].erase(edge);
		index_edges_stop_y[edge->to_y].erase(edge);

		index_edges_from[edge->from_x][edge->from_y][edge->from_z].erase(edge);
		index_edges_to[edge->to_x][edge->to_y][edge->to_z].erase(edge);


		if(edge->usedGraphEdge != 0)
		{
			index_graphEdges_2_NWEdges[edge->usedGraphEdge].erase(edge);
		}
	}



	paths.erase(p);
	p->table = 0;
}

bool VirtualNWTable_Unique::hasEdgeEmanatingFrom(int x, int y, int z)
{
	if((index_edges_from.count(x) == 0) || (index_edges_from.at(x).count(y) == 0) || (index_edges_from.at(x).at(y).count(z) == 0))
	{
		return false;
	}
	else
	{
		return true;
	}

//	std::set<NWEdge*>& candidateEdges = index_edges_from.at(x).at(y).at(z);
//	for(std::set<NWEdge*>::iterator edgeIt = candidateEdges.begin(); edgeIt != candidateEdges.end(); edgeIt++)
//	{
//		NWEdge* e = *edgeIt;
//		assert(e->from_x == x);
//		if((e->from_y == y) && (e->from_z == z))
//		{
//			return true;
//		}
//	}
//	return false;
}

bool VirtualNWTable_Unique::hasEdgeGoingInto(int x, int y, int z)
{
	if((index_edges_to.count(x) == 0) || (index_edges_to.at(x).count(y) == 0) || (index_edges_to.at(x).at(y).count(z) == 0))
	{
		return false;
	}
	else
	{
		return true;
	}


//	std::set<NWEdge*>& candidateEdges = index_edges_to.at(x).at(y).at(z);
//	for(std::set<NWEdge*>::iterator edgeIt = candidateEdges.begin(); edgeIt != candidateEdges.end(); edgeIt++)
//	{
//		NWEdge* e = *edgeIt;
//		assert(e->to == x);
//		if((e->to == y) && (e->from_z == z))
//		{
//			return true;
//		}
//	}
//	return false;
}

std::set<NWEdge*> VirtualNWTable_Unique::getEdgesEmanatingFrom(int x, int y, int z)
{
	if((index_edges_from.count(x) == 0) || (index_edges_from.at(x).count(y) == 0) || (index_edges_from.at(x).at(y).count(z) == 0))
	{
		return std::set<NWEdge*>();
	}
	else
	{
		return index_edges_from.at(x).at(y).at(z);
	}
}

std::set<NWEdge*> VirtualNWTable_Unique::getEdgesGoingInto(int x, int y, int z)
{
	if((index_edges_to.count(x) == 0) || (index_edges_to.at(x).count(y) == 0) || (index_edges_to.at(x).at(y).count(z) == 0))
	{
		return std::set<NWEdge*>();
	}
	else
	{
		return index_edges_to.at(x).at(y).at(z);
	}

}

void VirtualNWTable_Unique::_testTracePath(NWPath* p, std::string& reconstructed_sequence, std::string& reconstructed_graph, std::vector<int>& reconstructed_graph_levels)
{
	assert(paths.count(p) > 0);
	p->checkConsistency();

	reconstructed_sequence = "";
	reconstructed_graph = "";
	reconstructed_graph_levels.clear();

	NWEdge* thisStepEdge;
	NWEdge* nextStepEdge = *(p->first_edges.begin());

	do {
		thisStepEdge = nextStepEdge;

		int this_x = thisStepEdge->from_x;
		int this_y = thisStepEdge->from_y;

		int next_x = thisStepEdge->to_x;
		int next_y = thisStepEdge->to_y;

		Edge* usedGraphEdge = thisStepEdge->usedGraphEdge;

		std::string edgeEmission;
		if(usedGraphEdge != 0)
		{
			edgeEmission = usedGraphEdge->getEmission();
			assert(edgeEmission.size() == 1);
		}

		std::string sequenceEmission;
		if(next_y >= 1)
		{
			sequenceEmission = Context->sequence->substr(next_y - 1, 1);
		}

		if((next_x == (this_x + 1)) && (next_y == (this_y + 1)))
		{
			// match or mismatch
			assert(usedGraphEdge != 0);
			reconstructed_graph.append(edgeEmission);
			reconstructed_graph_levels.push_back(this_x);
			reconstructed_sequence.append(sequenceEmission);
		}
		else if((next_x == this_x) && (next_y == (this_y + 1)))
		{
			// gap in graph
			reconstructed_graph.append("_");
			reconstructed_graph_levels.push_back(-1);
			reconstructed_sequence.append(sequenceEmission);
		}
		else if((next_x == (this_x + 1)) && (next_y == this_y))
		{
			// gap in sequence
			assert(usedGraphEdge != 0);
			reconstructed_graph.append(edgeEmission);
			reconstructed_graph_levels.push_back(this_x);
			reconstructed_sequence.append("_");
		}
		else
		{
			assert(1 == 0);
		}

		std::set<NWEdge*> nextPossibleEgdes = p->forwardEdges(thisStepEdge);
		if(nextPossibleEgdes.size() > 0)
		{
			nextStepEdge = *(nextPossibleEgdes.begin());
		}
		else
		{
			assert(p->last_edges.count(thisStepEdge) > 0);
		}

	} while(p->last_edges.count(thisStepEdge) == 0);
}

void VirtualNWTable_Unique::checkConsistency()
{
	for(std::set<NWPath*>::iterator pathIt = paths.begin(); pathIt != paths.end(); pathIt++)
	{
		NWPath *p = *pathIt;
		assert(p->table == this);
		p->checkConsistency();
	}
}


unsigned int VirtualNWTable_Unique::getNumEdges()
{
	unsigned int forReturn = 0;
	for(std::set<NWPath*>::iterator pIt = paths.begin(); pIt != paths.end(); pIt++)
	{
		NWPath* p = *pIt;
		forReturn += p->getNumEdges();
	}
	return forReturn;
}

unsigned int VirtualNWTable_Unique::getNumEntryEgdes()
{
	unsigned int forReturn = 0;
	for(std::set<NWPath*>::iterator pIt = paths.begin(); pIt != paths.end(); pIt++)
	{
		NWPath* p = *pIt;
		forReturn += p->entry_edges.size();
	}
	return forReturn;
}

std::vector<NWEdge*> VirtualNWTable_Unique::getEntryEdges()
{
	std::set<NWEdge*> forReturn_set;
	for(std::set<NWPath*>::iterator pIt = paths.begin(); pIt != paths.end(); pIt++)
	{
		NWPath* p = *pIt;
		forReturn_set.insert(p->entry_edges.begin(), p->entry_edges.end());
	}
	std::vector<NWEdge*> forReturn(forReturn_set.begin(), forReturn_set.end());
	return forReturn;
}

std::vector<NWEdge*> VirtualNWTable_Unique::getExitEdges()
{
	std::set<NWEdge*> forReturn_set;
	for(std::set<NWPath*>::iterator pIt = paths.begin(); pIt != paths.end(); pIt++)
	{
		NWPath* p = *pIt;
		forReturn_set.insert(p->exit_edges.begin(), p->exit_edges.end());
	}
	std::vector<NWEdge*> forReturn(forReturn_set.begin(), forReturn_set.end());
	return forReturn;
}

std::set<NWEdge*> VirtualNWTable_Unique::getAllEdges()
{
	std::set<NWEdge*> forReturn_set;
	for(std::set<NWPath*>::iterator pIt = paths.begin(); pIt != paths.end(); pIt++)
	{
		NWPath* p = *pIt;
		forReturn_set.insert(p->all_edges.begin(), p->all_edges.end());
	}
	return forReturn_set;
}


unsigned int VirtualNWTable_Unique::getNumExitEdges()
{
	unsigned int forReturn = 0;
	for(std::set<NWPath*>::iterator pIt = paths.begin(); pIt != paths.end(); pIt++)
	{
		NWPath* p = *pIt;
		forReturn += p->exit_edges.size();
	}
	return forReturn;
}



NWPath::NWPath()
{
	table = 0;
}

std::vector<Edge*> NWPath::graphEdgesPath()
{
	assert(first_edges.size() == 1);
	NWEdge* currentEdge = *(first_edges.begin());

	std::vector<Edge*> forReturn;
	if(last_edges.count(currentEdge) > 0)
	{
		forReturn.push_back(currentEdge->usedGraphEdge);
		return forReturn;
	}

	do
	{
		forReturn.push_back(currentEdge->usedGraphEdge);
		assert(edges_from.at(currentEdge->to_x).at(currentEdge->to_y).at(currentEdge->to_z).size() == 1);
		currentEdge = *(edges_from.at(currentEdge->to_x).at(currentEdge->to_y).at(currentEdge->to_z).begin());

	} while(last_edges.count(currentEdge) == 0);

	assert(last_edges.count(currentEdge) > 0);

	forReturn.push_back(currentEdge->usedGraphEdge);

	return forReturn;
}

void NWEdge::checkConsistency()
{
	assert(to_x >= 0);
	assert(to_y >= 0);
	assert(to_z >= 0);

	assert(to_x >= from_x);
	assert(to_y >= from_y);
	assert(from_z >= 0);

	assert((to_x - from_x) <= 1);
	assert((to_y - from_y) <= 1);

	if(path != 0)
	{
		NWEdge* selfPointer = this;
		assert(path->edges_from.at(from_x).at(from_y).at(from_z).count(selfPointer) > 0);
		assert(path->edges_to.at(to_x).at(to_y).at(to_z).count(selfPointer) > 0);

		if(path->table != 0)
		{
			if(! (to_y <= (int)path->table->Context->sequence->length()))
			{
				std::cerr << "to_y: "<< to_y << "\n";
				std::cerr << "path->table->Context->sequence->length(): "<< path->table->Context->sequence->length() << "\n" << std::flush;
			}
			assert(to_y <= (int)path->table->Context->sequence->length());
			assert(to_x < (int)path->table->Context->aligner->nodesPerLevel_ordered.size());

			assert(from_z < (int)path->table->Context->aligner->nodesPerLevel_ordered.at(from_x).size());
			assert(to_z < (int)path->table->Context->aligner->nodesPerLevel_ordered.at(to_x).size());


			assert(path->table->index_edges_start_x.at(from_x).count(selfPointer) > 0);
			assert(path->table->index_edges_stop_x.at(to_x).count(selfPointer) > 0);
			assert(path->table->index_edges_start_y.at(from_y).count(selfPointer) > 0);
			assert(path->table->index_edges_stop_y.at(to_y).count(selfPointer) > 0);

			assert(path->table->index_edges_from.at(from_x).at(from_y).at(from_z).count(selfPointer) > 0);
			assert(path->table->index_edges_to.at(to_x).at(to_y).at(to_z).count(selfPointer) > 0);

			if(usedGraphEdge != 0)
			{
				assert(path->table->Context->aligner->g->Edges.count(usedGraphEdge) > 0);
				assert(path->table->Context->aligner->nodesPerLevel_ordered.at(from_x).at(from_z)->Outgoing_Edges.count(usedGraphEdge) > 0);
				if(! (path->table->Context->aligner->nodesPerLevel_ordered.at(to_x).at(to_z)->Incoming_Edges.count(usedGraphEdge) > 0))
				{
					std::cerr << "Problem with path: specified edge not present in graph! to_x = " << to_x << ", to_z = " << to_z << "\n" << std::flush;
				}
				assert(path->table->Context->aligner->nodesPerLevel_ordered.at(to_x).at(to_z)->Incoming_Edges.count(usedGraphEdge) > 0);

				assert(path->table->index_graphEdges_2_NWEdges[usedGraphEdge].count(selfPointer) > 0);
			}
		}
	}
}

bool NWPath::edgeExists(int from_x, int from_y, int from_z, int to_x, int to_y, int to_z, Edge* graphEdge)
{
	if(
			(edges_from.count(from_x) && edges_from.at(from_x).count(from_y) && edges_from.at(from_x).at(from_y).count(from_z)) &&
			(edges_to.count(to_x) && edges_to.at(to_x).count(to_y) && edges_to.at(to_x).at(to_y).count(to_z))
	)
	{
		std::set<NWEdge*>& candidates_from = edges_from.at(from_x).at(from_y).at(from_z);
		std::set<NWEdge*>& candidates_to =  edges_to.at(to_x).at(to_y).at(to_z);
		for(std::set<NWEdge*>::iterator eIt = candidates_from.begin(); eIt != candidates_from.end(); eIt++)
		{
			NWEdge* e = *eIt;
			if(candidates_to.count(e) > 0)
			{
				if(e->usedGraphEdge == graphEdge)
				{
					return true;
				}
			}
		}
		return false;
	}
	else
	{
		return false;
	}
}

NWEdge* NWPath::retrieveEdge(int from_x, int from_y, int from_z, int to_x, int to_y, int to_z, Edge* graphEdge)
{
	std::set<NWEdge*>& candidates_from = edges_from.at(from_x).at(from_y).at(from_z);
	std::set<NWEdge*>& candidates_to =  edges_to.at(to_x).at(to_y).at(to_z);
	for(std::set<NWEdge*>::iterator eIt = candidates_from.begin(); eIt != candidates_from.end(); eIt++)
	{
		NWEdge* e = *eIt;
		if(candidates_to.count(e) > 0)
		{
			if(e->usedGraphEdge == graphEdge)
			{
				return e;
			}
		}
	}

	std::cerr << "Error: could not retrieve edge!\n";

	assert( 1 == 0 ); // should not happen!
	return 0;
}

void NWPath::extendToZero(int whichEnd, int lastYCoordinate)
{
	if(whichEnd == -1)
	{
		std::set<NWEdge*> entryEdges_copy = entry_edges;
		for(std::set<NWEdge*>::iterator entryEdgeIt = entryEdges_copy.begin(); entryEdgeIt != entryEdges_copy.end(); entryEdgeIt++)
		{
			NWEdge* entryEdge = *entryEdgeIt;
			if(entryEdge->from_y != 0)
			{
				entry_edges.erase(entryEdge);
				// int missingYs = entryEdge->from_y;
				for(int y = 1; y <= entryEdge->from_y; y++)
				{
					int entryExitStatus = (y == 1) ? -1 : 0;
					createAndAddEdge(entryEdge->from_x, y - 1, entryEdge->from_z, entryEdge->from_x, y, entryEdge->from_z, 0, 0, entryExitStatus);
					// std::cerr << "CreateAdd edge (" << entryEdge->from_x << ", " << (y - 1) << ", " << entryEdge->from_z << ") -> " <<
				}
			}
		}
	}
	else
	{
		assert(whichEnd == 1);
		assert(lastYCoordinate != 0);

		std::set<NWEdge*> exitEdges_copy = exit_edges;

		for(std::set<NWEdge*>::iterator exitEdgeIt = exitEdges_copy.begin(); exitEdgeIt != exitEdges_copy.end(); exitEdgeIt++)
		{
			NWEdge* exitEdge = *exitEdgeIt;
			assert(exitEdge != 0);

			if(exitEdge->to_y != lastYCoordinate)
			{
				exit_edges.erase(exitEdge);

				for(int y = exitEdge->to_y + 1; y <= lastYCoordinate; y++)
				{
					int entryExitStatus = (y == lastYCoordinate) ? 1 : 0;
					createAndAddEdge(exitEdge->to_x, y - 1, exitEdge->to_z, exitEdge->to_x, y, exitEdge->to_z, 0, 0, entryExitStatus);
				}
			}
		}
	}
}

void NWPath::takeInExtensionPath(localExtension_pathDescription* pD, int entryExit, int completeToEnd)
{
	assert((entryExit == -1) || (entryExit == 1));
	assert((completeToEnd == 0) || (completeToEnd == 1) || (completeToEnd == -1));

	assert(pD->usedEdges.size() > 0);
	assert(pD->coordinates.size() == (pD->usedEdges.size() + 1));
	if(completeToEnd == 0)
	{
		for(unsigned int coordinateI = 0; coordinateI < (pD->coordinates.size() - 1); coordinateI++)
		{
			std::vector<int>& thisCoordinates = pD->coordinates.at(coordinateI);
			std::vector<int>& nextCoordinates = pD->coordinates.at(coordinateI+1);
			Edge* usedGraphEdge = pD->usedEdges.at(coordinateI);
			if(! edgeExists(thisCoordinates.at(0), thisCoordinates.at(1), thisCoordinates.at(2), nextCoordinates.at(0), nextCoordinates.at(1), nextCoordinates.at(2), usedGraphEdge))
			{
				int entryExitStatus = 0;
				if((entryExit == -1) && (coordinateI == 0))
				{
					entryExitStatus = -1;
				}
				if((entryExit == 1) && (coordinateI == (pD->coordinates.size() - 2)))
				{
					entryExitStatus = 1;
				}

				createAndAddEdge(thisCoordinates.at(0), thisCoordinates.at(1), thisCoordinates.at(2), nextCoordinates.at(0), nextCoordinates.at(1), nextCoordinates.at(2), usedGraphEdge, 0, entryExitStatus);
			}
			else
			{
				int entryExitStatus = 0;
				if((entryExit == -1) && (coordinateI == 0))
				{
					entryExitStatus = -1;
				}
				if((entryExit == 1) && (coordinateI == (pD->coordinates.size() - 2)))
				{
					entryExitStatus = 1;
				}
				if(entryExitStatus != 0)
				{
					addStatusEdge(thisCoordinates.at(0), thisCoordinates.at(1), thisCoordinates.at(2), nextCoordinates.at(0), nextCoordinates.at(1), nextCoordinates.at(2), usedGraphEdge, 0, entryExitStatus);
				}

			}
		}
	}
	else
	{
		for(unsigned int coordinateI = 0; coordinateI < (pD->coordinates.size() - 1); coordinateI++)
		{
			std::vector<int>& thisCoordinates = pD->coordinates.at(coordinateI);
			std::vector<int>& nextCoordinates = pD->coordinates.at(coordinateI+1);
			Edge* usedGraphEdge = pD->usedEdges.at(coordinateI);
			if(! edgeExists(thisCoordinates.at(0), thisCoordinates.at(1), thisCoordinates.at(2), nextCoordinates.at(0), nextCoordinates.at(1), nextCoordinates.at(2), usedGraphEdge))
			{
				int entryExitStatus = 0;
				createAndAddEdge(thisCoordinates.at(0), thisCoordinates.at(1), thisCoordinates.at(2), nextCoordinates.at(0), nextCoordinates.at(1), nextCoordinates.at(2), usedGraphEdge, 0, entryExitStatus);
			}
			else
			{
				// nothing
			}
		}

		if(completeToEnd == 1)
		{
			std::vector<int> lastCoordinates = pD->coordinates.at(pD->coordinates.size() - 1);

			int lastYCoordinate = table->Context->sequence->length();
			int missingYs = table->Context->sequence->length() - lastCoordinates.at(1);
			assert(missingYs >= 0);

			for(int y = lastCoordinates.at(1) + 1; y <= lastYCoordinate; y++)
			{
				int entryExitStatus = 0;
				std::cerr << "[completeToEnd == 1] Add edge " << lastCoordinates.at(0) << ", " << (y-1) << ", " << lastCoordinates.at(2) << " to " << lastCoordinates.at(0) << ", " << y << ", " << lastCoordinates.at(2) << "\n" << std::flush;
				createAndAddEdge(lastCoordinates.at(0), y - 1, lastCoordinates.at(2), lastCoordinates.at(0), y, lastCoordinates.at(2), 0, 0, entryExitStatus);
			}

			if(entryExit != 0)
			{
				if(missingYs > 0)
				{
					addStatusEdge(lastCoordinates.at(0), lastYCoordinate - 1, lastCoordinates.at(2), lastCoordinates.at(0), lastYCoordinate, lastCoordinates.at(2), 0, 0, entryExit);
				}
				else
				{
					Edge* usedGraphEdge = pD->usedEdges.at(pD->usedEdges.size() - 1);
					std::vector<int> secondLastCoordinates = pD->coordinates.at(pD->coordinates.size() - 1);
					addStatusEdge(secondLastCoordinates.at(0), secondLastCoordinates.at(1), secondLastCoordinates.at(2), lastCoordinates.at(0), lastCoordinates.at(1), lastCoordinates.at(2), usedGraphEdge, 0, entryExit);
				}
			}
		}
		else
		{
			assert(completeToEnd == -1);
			std::vector<int> firstCoordinates = pD->coordinates.at(0);

			int missingYs = firstCoordinates.at(1);
			assert(missingYs >= 0);


			for(int y = 1; y <= firstCoordinates.at(1); y++)
			{
				int entryExitStatus = 0;
				std::cerr << "[completeToEnd == -1] Add edge " << firstCoordinates.at(0) << ", " << (y-1) << ", " << firstCoordinates.at(2) << " to " << firstCoordinates.at(0) << ", " << y << ", " << firstCoordinates.at(2) << "\n" << std::flush;
				createAndAddEdge(firstCoordinates.at(0), y - 1, firstCoordinates.at(2), firstCoordinates.at(0), y, firstCoordinates.at(2), 0, 0, entryExitStatus);
			}

			if(entryExit != 0)
			{
				if(missingYs > 0)
				{
					addStatusEdge(firstCoordinates.at(0), 0, firstCoordinates.at(2), firstCoordinates.at(0), 1, firstCoordinates.at(2), 0, 0, entryExit);
				}
				else
				{
					// Edge* usedGraphEdge = pD->usedEdges.at(0);

					std::vector<int> secondCoordinates = pD->coordinates.at(1);
					addStatusEdge(firstCoordinates.at(0), firstCoordinates.at(1), firstCoordinates.at(2), secondCoordinates.at(0), secondCoordinates.at(1), secondCoordinates.at(2), 0, 0, entryExit);
				}
			}
		}

	}
//	std::cerr << "takeInExtensionPath: entryExit is " << entryExit << " and # entry edges = " << entry_edges.size() << " / # exit edges: " << exit_edges.size() << "\n" << std::flush;
}

void NWPath::eraseAllEdgeStatus()
{
	first_edges.clear();
	last_edges.clear();
	entry_edges.clear();
	exit_edges.clear();
}

NWPath* NWPath::clonePathWithoutTable()
{
	NWPath* newPath = new NWPath();

	for(std::set<NWEdge*>::iterator eIt = all_edges.begin(); eIt != all_edges.end(); eIt++)
	{
		NWEdge* existingE = *eIt;
		NWEdge* newEdge = new NWEdge(*existingE);
		assert(newEdge->path == this);

		newEdge->path = newPath;

		newPath->all_edges.insert(newEdge);
		if(first_edges.count(existingE))
		{
			newPath->first_edges.insert(newEdge);
		}
		if(last_edges.count(existingE))
		{
			newPath->last_edges.insert(newEdge);
		}
		if(entry_edges.count(existingE))
		{
			newPath->entry_edges.insert(newEdge);
		}
		if(exit_edges.count(existingE))
		{
			newPath->exit_edges.insert(newEdge);
		}
		newPath->edges_from[existingE->from_x][existingE->from_y][existingE->from_z].insert(newEdge);
		newPath->edges_to[existingE->to_x][existingE->to_y][existingE->to_z].insert(newEdge);
	}

	newPath->checkConsistency();

	return newPath;
}

void NWPath::eraseStatusEdge(int from_x, int from_y, int from_z, int to_x, int to_y, int to_z, Edge* graphEdge)
{
	NWEdge* e = retrieveEdge(from_x, from_y, from_z, to_x, to_y, to_z, graphEdge);
	first_edges.erase(e);
	last_edges.erase(e);
	entry_edges.erase(e);
	exit_edges.erase(e);
}


void NWPath::addStatusEdge(int from_x, int from_y, int from_z, int to_x, int to_y, int to_z, Edge* graphEdge, int firstOrLast, int entryExitStatus)
{
	NWEdge* e = retrieveEdge(from_x, from_y, from_z, to_x, to_y, to_z, graphEdge);

	if(firstOrLast == -1)
	{
		first_edges.insert(e);
	}
	if(firstOrLast == 1)
	{
		last_edges.insert(e);
	}
	if(entryExitStatus == -1)
	{
		entry_edges.insert(e);
	}
	if(entryExitStatus == 1)
	{
		exit_edges.insert(e);
	}
}

NWEdge* NWPath::createAndAddEdge(int from_x, int from_y, int from_z, int to_x, int to_y, int to_z, Edge* graphEdge, int firstOrLast, int entryExitStatus)
{
	NWEdge* edge = new NWEdge();
	edge->from_x = from_x;
	edge->from_y = from_y;
	edge->from_z = from_z;
	edge->to_x = to_x;
	edge->to_y = to_y;
	edge->to_z = to_z;
	edge->usedGraphEdge = graphEdge;
	addEdge(edge);

	assert((firstOrLast == -1) || (firstOrLast == 0) || (firstOrLast == 1));
	assert((entryExitStatus == -1) || (entryExitStatus == 0) || (entryExitStatus == 1));

	if(firstOrLast == -1)
	{
		first_edges.insert(edge);
	}
	if(firstOrLast == 1)
	{
		last_edges.insert(edge);
	}
	if(entryExitStatus == -1)
	{
		entry_edges.insert(edge);
	}
	if(entryExitStatus == 1)
	{
		exit_edges.insert(edge);
	}

	return edge;
}


void NWPath::addEdge(NWEdge* e)
{
	assert(table == 0);
	assert(all_edges.count(e) == 0);

	e->path = this;
	edges_from[e->from_x][e->from_y][e->from_z].insert(e);
	edges_to[e->to_x][e->to_y][e->to_z].insert(e);
	all_edges.insert(e);
}


std::set<NWEdge*> NWPath::forwardEdges(int x, int y, int z)
{
	if(edges_from.count(x) && edges_from.at(x).count(y) && edges_from.at(x).at(y).count(z))
	{
		return edges_from.at(x).at(y).at(z);
	}
	else
	{
		return std::set<NWEdge*>();
	}
}

std::set<NWEdge*> NWPath::forwardEdges(NWEdge* e)
{
	assert(e->path == this);
	return forwardEdges(e->to_x, e->to_y, e->to_z);
}

std::set<NWEdge*> NWPath::backwardEdges(int x, int y, int z)
{
	if(edges_to.count(x) && edges_to.at(x).count(y) && edges_to.at(x).at(y).count(z))
	{
		return edges_to.at(x).at(y).at(z);
	}
	else
	{
		return std::set<NWEdge*>();
	}
}

std::set<NWEdge*> NWPath::backwardEdges(NWEdge* e)
{
	assert(e->path == this);
	return backwardEdges(e->from_x, e->from_y, e->from_z);
}

void NWPath::recalculateFirstLast()
{
	first_edges.clear();
	last_edges.clear();
	for(std::set<NWEdge*>::iterator edgeIt = all_edges.begin(); edgeIt != all_edges.end(); edgeIt++)
	{
		NWEdge* e = *edgeIt;
		if(forwardEdges(e).size() == 0)
		{
			last_edges.insert(e);
		}
		if(backwardEdges(e).size() == 0)
		{
			first_edges.insert(e);
		}
	}
}

void NWPath::verify_edges()
{
	std::set<NWEdge*> seen_edges_fromTo;

	for(std::map<int, std::map<int, std::map<int, std::set<NWEdge*> > > >::iterator xIt = edges_from.begin(); xIt != edges_from.end(); xIt++)
	{
		for(std::map<int, std::map<int, std::set<NWEdge*> > >::iterator yIt = xIt->second.begin(); yIt != xIt->second.end(); yIt++)
		{
			for(std::map<int, std::set<NWEdge*>>::iterator zIt = yIt->second.begin(); zIt != yIt->second.end(); zIt++)
			{
				for(std::set<NWEdge*>::iterator edgeIt = zIt->second.begin(); edgeIt != zIt->second.end(); edgeIt++)
				{
					NWEdge* e = *edgeIt;
					assert(all_edges.count(e) > 0);
					seen_edges_fromTo.insert(e);
				}
			}
		}
	}

	for(std::map<int, std::map<int, std::map<int, std::set<NWEdge*> > > >::iterator xIt = edges_to.begin(); xIt != edges_to.end(); xIt++)
	{
		for(std::map<int, std::map<int, std::set<NWEdge*> > >::iterator yIt = xIt->second.begin(); yIt != xIt->second.end(); yIt++)
		{
			for(std::map<int, std::set<NWEdge*>>::iterator zIt = yIt->second.begin(); zIt != yIt->second.end(); zIt++)
			{
				for(std::set<NWEdge*>::iterator edgeIt = zIt->second.begin(); edgeIt != zIt->second.end(); edgeIt++)
				{
					NWEdge* e = *edgeIt;
					assert(all_edges.count(e) > 0);
					seen_edges_fromTo.insert(e);
				}
			}
		}
	}

	for(std::set<NWEdge*>::iterator edgeIt = all_edges.begin(); edgeIt != all_edges.end(); edgeIt++)
	{
		NWEdge* e = *edgeIt;
		assert(seen_edges_fromTo.count(e) > 0);
		assert(e->path == this);

		e->checkConsistency();
		if(forwardEdges(e).size() == 0)
		{
			if(!(last_edges.count(e) > 0))
			{
				std::cerr << "!(last_edges.count(e) > 0)" << "\n";
				std::cerr << "e: " << e << "\n";
				_printPath();
				std::cerr << std::flush;
			}
			assert(last_edges.count(e) > 0);
		}
		else
		{
			if(last_edges.count(e) != 0)
			{
				std::cerr << "Forward edges, but a last edge?\n";
				std::cerr << "last_edges.size(): " << last_edges.size() << "\n";
				std::cerr << "forwardEdges(e).size(): " << forwardEdges(e).size() << "\n";

				e->_printEdge();

				_printPath();
			}
			assert(last_edges.count(e) == 0);
		}
		if(backwardEdges(e).size() == 0)
		{
			assert(first_edges.count(e) > 0);
		}
		else
		{
			assert(first_edges.count(e) == 0);
		}
	}

	for(std::set<NWEdge*>::iterator firstEdgeIt = first_edges.begin(); firstEdgeIt != first_edges.end(); firstEdgeIt++)
	{
		NWEdge* firstEdge = *firstEdgeIt;
		assert(all_edges.count(firstEdge) > 0);
		assert(entry_edges.count(firstEdge) > 0);
	}

	for(std::set<NWEdge*>::iterator lastEdgeIt = last_edges.begin(); lastEdgeIt != last_edges.end(); lastEdgeIt++)
	{
		NWEdge* lastEdge = *lastEdgeIt;
		assert(all_edges.count(lastEdge) > 0);
		assert(exit_edges.count(lastEdge) > 0);
	}

	for(std::set<NWEdge*>::iterator entryEdgeIt = entry_edges.begin(); entryEdgeIt != entry_edges.end(); entryEdgeIt++)
	{
		NWEdge* entryEdge = *entryEdgeIt;
		assert(all_edges.count(entryEdge) > 0);
	}

	for(std::set<NWEdge*>::iterator exitEdgeIt = exit_edges.begin(); exitEdgeIt != exit_edges.end(); exitEdgeIt++)
	{
		NWEdge* exitEdge = *exitEdgeIt;
		assert(all_edges.count(exitEdge) > 0);
	}

	assert(first_edges.size() > 0);
	assert(last_edges.size() > 0);
	assert(entry_edges.size() > 0);
	assert(exit_edges.size() > 0);

//	for(std::set<NWEdge*>::iterator exitEdgeIt = exit_edges.begin(); exitEdgeIt != exit_edges.end(); exitEdgeIt++)
//	{
//		NWEdge* exitEdge = *exitEdgeIt;
//		for(std::set<NWEdge*>::iterator entryEdgeIt = entry_edges.begin(); entryEdgeIt != entry_edges.end(); entryEdgeIt++)
//		{
//			NWEdge* entryEdge = *entryEdgeIt;
//			assert(exitEdge->from_x > entryEdge->to_x);
//			assert(exitEdge->from_y > entryEdge->to_y);
//		}
//	}

}

void NWEdge::makeExitEdge()
{
	assert(path != 0);
	path->exit_edges.insert(this);
}


void NWPath::_printPath()
{
	std::cout << "NWPath " << this << "\n\n";
	for(std::map<int, std::map<int, std::map<int, std::set<NWEdge*> > > >::iterator xIt = edges_from.begin(); xIt != edges_from.end(); xIt++)
	{
		for(std::map<int, std::map<int, std::set<NWEdge*> > >::iterator yIt = xIt->second.begin(); yIt != xIt->second.end(); yIt++)
		{
			for(std::map<int, std::set<NWEdge*>>::iterator zIt = yIt->second.begin(); zIt != yIt->second.end(); zIt++)
			{
				for(std::set<NWEdge*>::iterator edgeIt = zIt->second.begin(); edgeIt != zIt->second.end(); edgeIt++)
				{
					NWEdge* e = *edgeIt;
					std::string firstString = first_edges.count(e) ?  " FIRST " : "";
					std::string lastString = last_edges.count(e) ?  " LAST " : "";
					std::string entryString = entry_edges.count(e) ?  " ENTRY " : "";
					std::string exitString = exit_edges.count(e) ?  " EXIT " : "";
					std::cout << "\t" << e << ": (" << e->from_x << ", " << e->from_y << ", " << e->from_z << ") -> (" << e->to_x << ", " << e->to_y << ", " << e->to_z << ") " << firstString << lastString << entryString << exitString << "\n";
				}
			}
		}
	}
	std::cout << "\n" << std::flush;
}

void NWPath::checkConsistency()
{
	verify_edges();
	std::set<NWEdge*> seen_edges;
	std::set<NWEdge*> visit_edges;
	if(all_edges.size() > 0)
	{
		visit_edges.insert(*(all_edges.begin()));
	}

	while(visit_edges.size() > 0)
	{
		NWEdge* e = *(visit_edges.begin());
		visit_edges.erase(e);
		if(! seen_edges.count(e))
		{
			seen_edges.insert(e);
			std::set<NWEdge*> nextVisit = forwardEdges(e);
			std::set<NWEdge*> nextVisit_2 = backwardEdges(e);
			visit_edges.insert(nextVisit.begin(), nextVisit.end());
			visit_edges.insert(nextVisit_2.begin(), nextVisit_2.end());
		}
	}

	assert(all_edges.size() == seen_edges.size());
}

void VirtualNWTable_Unique::freeMemory()
{
	for(std::set<NWPath*>::iterator pathIt = paths.begin(); pathIt != paths.end(); pathIt++)
	{
		NWPath* path = *pathIt;
		path->freeMemory();
		delete(path);
	}
}

void NWPath::freeMemory()
{
	for(std::set<NWEdge*>::iterator eIt = all_edges.begin(); eIt != all_edges.end(); eIt++)
	{
		NWEdge* e = *eIt;
		delete(e);
	}
}

void NWEdge::_printEdge()
{
	std::cout << "NWEdge " << this << "\n";
	std::cout << "\t" << "from_x: " << from_x << "\n";
	std::cout << "\t" << "from_y: " << from_y << "\n";
	std::cout << "\t" << "from_z: " << from_z << "\n";

	std::cout << "\t" << "to_x: " << to_x << "\n";
	std::cout << "\t" << "to_y: " << to_y << "\n";
	std::cout << "\t" << "to_z: " << to_z << "\n";

	std::cout << "\t" << "usedGraphEdge: " << usedGraphEdge << "\n" << std::flush;

}


void NWEdge::takeScore(double Score, NWEdge* backtrackEdge, double S_match, double S_mismatch, double S_gap)
{
	assert( 1 == 0);
//	std::cout << "enter takeScore...\n" << std::flush;
	auto scoreEdge = [&](NWEdge* e) -> double {
		assert(e != 0);

//		std::cout << "scoreEdge " << "1" << "\n" << std::flush;

		int diff_x = e->to_x - e->from_x;
		int diff_y = e->to_y - e->from_y;
		assert((diff_x >= 0) && (diff_x <= 1));
		assert((diff_y >= 0) && (diff_y <= 1));

//		std::cout << "scoreEdge " << "2" << "\n" << std::flush;

		if((diff_x == 1) && (diff_y == 1))
		{

//			std::cout << "scoreEdge " << "3.1" << "\n" << std::flush;

			assert(e->usedGraphEdge != 0);
			assert(e->path != 0);
			assert(e->path->table != 0);
			assert(e->path->table->Context->sequence != 0);

//			std::cout << "scoreEdge " << "3.1.1" << "\n" << std::flush;

			std::string edge_emission = e->usedGraphEdge->getEmission();

//			std::cout << "scoreEdge " << "3.1.2" << "\n" << std::flush;
			assert(edge_emission.length() == 1);
			std::string sequenceEmission = e->path->table->Context->sequence->substr(e->from_y, 1);
//			std::cout << "scoreEdge " << "3.1.3" << "\n" << std::flush;

			if(edge_emission == sequenceEmission)
			{
				return S_match;
			}
			else
			{
//				std::cerr << "OK " << edge_emission << " != " << sequenceEmission << " "  << S_mismatch << "\n";
				return S_mismatch;
			}
		}
		else
		{

//			std::cout << "scoreEdge " << "3.2" << "\n" << std::flush;

			if((diff_x == 1) && (diff_y == 0))
			{
				// gap in sequence
				assert(e->usedGraphEdge != 0);
				std::string edge_emission = e->usedGraphEdge->getEmission();
				assert(edge_emission.length() == 1);
				if(edge_emission == "_")
				{
					return 0;
				}
				else
				{
					return S_gap;
				}
			}
			else
			{
				assert((diff_x == 0) && (diff_y == 1));
				return S_gap;
			}
		}
	};

	assert(path != 0);
	std::set<NWEdge*> backwardEdges_set = path->backwardEdges(this);
	std::vector<NWEdge*> backwardEdges(backwardEdges_set.begin(), backwardEdges_set.end());

	std::vector<double> ScoreAlternatives;
	std::vector<NWEdge*> ScoreAlternatives_bt;

	for(unsigned int bwI = 0; bwI < backwardEdges.size(); bwI++)
	{
		NWEdge* bwE = backwardEdges.at(bwI);
		assert(bwE->scoreComputed);
		ScoreAlternatives.push_back(bwE->scoreAfterEdge);
		ScoreAlternatives_bt.push_back(bwE);
	}

	ScoreAlternatives.push_back(Score);
	ScoreAlternatives_bt.push_back(backtrackEdge);

	std::pair<double, unsigned int> maxScore = Utilities::findVectorMaxP(ScoreAlternatives);

	double selectedScore = maxScore.first;
	NWEdge* selectedBacktrack = ScoreAlternatives_bt.at(maxScore.second);

	scoreAfterEdge = selectedScore + scoreEdge(this);
	scoreBacktrack = selectedBacktrack;
	scoreComputed = true;

//	std::cout << "\t halfway takeScore...\n" << std::flush;

	auto canComputeEdge = [&](NWEdge* e) -> bool {
		assert(e != 0);
		if(e->path->entry_edges.count(e))
		{
			return false;
		}
		else
		{
			std::set<NWEdge*> e_backward = e->path->backwardEdges(e);
			bool all_computed = true;
			for(std::set<NWEdge*>::iterator eIt = e_backward.begin(); eIt != e_backward.end(); eIt++)
			{
				NWEdge* previousEdge = *eIt;
				assert(previousEdge != 0);
				if(! previousEdge->scoreComputed)
				{
					all_computed = false;
					break;
				}
			}
			return all_computed;
		}
	};

//	std::cout << "\t B1...\n" << std::flush;
	std::set<NWEdge*> checkForComputability = path->forwardEdges(this);
//	std::cout << "\t B2...\n" << std::flush;

	while(checkForComputability.size() > 0)
	{

		NWEdge* thisE = *(checkForComputability.begin());
		assert(thisE != 0);
//		std::cout << "\t\tCheck edge " << thisE << "\n" << std::flush;

		checkForComputability.erase(thisE);
		if(canComputeEdge(thisE))
		{
//			std::cout << "\t\t\t" << "c1" << "\n" << std::flush;


			std::vector<double> ScoreAlternatives;
			std::vector<NWEdge*> ScoreAlternatives_bt;

			std::set<NWEdge*> thisE_backward_edges_set = thisE->path->backwardEdges(thisE);
			std::vector<NWEdge*> thisE_backward_edges(thisE_backward_edges_set.begin(), thisE_backward_edges_set.end());

//			std::cout << "\t\t\t" << "c2" << "\n" << std::flush;


			for(unsigned int bwI = 0; bwI < thisE_backward_edges.size(); bwI++)
			{
				NWEdge* bwE = thisE_backward_edges.at(bwI);
				assert(bwE != 0);
				assert(bwE->scoreComputed);
				ScoreAlternatives.push_back(bwE->scoreAfterEdge);
				ScoreAlternatives_bt.push_back(bwE);
			}

//			std::cout << "\t\t\t" << "c3" << "\n" << std::flush;

			assert(ScoreAlternatives.size() > 0);

			std::pair<double, unsigned int> maxScore = Utilities::findVectorMaxP(ScoreAlternatives);

			double selectedScore = maxScore.first;
			NWEdge* selectedBacktrack = ScoreAlternatives_bt.at(maxScore.second);

//			std::cout << "\t\t\t" << "c4" << "\n" << std::flush;


			thisE->scoreAfterEdge = selectedScore + scoreEdge(thisE);

//			std::cout << "\t\t\t" << "c4.0" << "\n" << std::flush;


			thisE->scoreBacktrack = selectedBacktrack;
			thisE->scoreComputed = true;

//			std::cout << "\t\t\t" << "c4.1" << "\n" << std::flush;


			assert(thisE->path != 0);
			std::set<NWEdge*> next_checkForComputability = thisE->path->forwardEdges(thisE);

//			std::cout << "\t\t\t" << "c4.2" << "\n" << std::flush;


			checkForComputability.insert(next_checkForComputability.begin(), next_checkForComputability.end());

//			std::cout << "\t\t\t" << "c5" << "\n" << std::flush;
		}
	}

//	std::cout << "\t\texit takeScore...\n" << std::flush;
}


void NWEdge::takeScore_affine(double Score, NWEdge* backtrackEdge, double S_match, double S_mismatch, double S_gapOpen, double S_gapExtend)
{
	assert( 1 == 0);
	assert(path != 0);
	std::set<NWEdge*> backwardEdges_set = path->backwardEdges(this);
	std::vector<NWEdge*> backwardEdges(backwardEdges_set.begin(), backwardEdges_set.end());

	std::vector<double> ScoreAlternatives;
	std::vector<NWEdge*> ScoreAlternatives_bt;

	for(unsigned int bwI = 0; bwI < backwardEdges.size(); bwI++)
	{
		NWEdge* bwE = backwardEdges.at(bwI);
		assert(bwE->scoreComputed);
		double thisEdgeScoreAfterThatEdge = calculateScore(bwE->isSequenceGap_affine(), bwE->isGraphGap(), S_match, S_mismatch, S_gapOpen, S_gapExtend);
		double edgesCombinedScore = thisEdgeScoreAfterThatEdge + bwE->scoreAfterEdge;
		ScoreAlternatives.push_back(edgesCombinedScore);
		ScoreAlternatives_bt.push_back(bwE);
	}

	ScoreAlternatives.push_back(Score);
	ScoreAlternatives_bt.push_back(backtrackEdge);

	std::pair<double, unsigned int> maxScore = Utilities::findVectorMaxP(ScoreAlternatives);

	double selectedScore = maxScore.first;
	NWEdge* selectedBacktrack = ScoreAlternatives_bt.at(maxScore.second);

	scoreAfterEdge = selectedScore;
	scoreBacktrack = selectedBacktrack;
	scoreComputed = true;

//	std::cout << "\t halfway takeScore...\n" << std::flush;

	auto canComputeEdge = [&](NWEdge* e) -> bool {
		assert(e != 0);
		if(e->path->entry_edges.count(e))
		{
			return false;
		}
		else
		{
			std::set<NWEdge*> e_backward = e->path->backwardEdges(e);
			bool all_computed = true;
			for(std::set<NWEdge*>::iterator eIt = e_backward.begin(); eIt != e_backward.end(); eIt++)
			{
				NWEdge* previousEdge = *eIt;
				assert(previousEdge != 0);
				if(! previousEdge->scoreComputed)
				{
					all_computed = false;
					break;
				}
			}
			return all_computed;
		}
	};

	std::set<NWEdge*> checkForComputability = path->forwardEdges(this);

	while(checkForComputability.size() > 0)
	{

		NWEdge* thisE = *(checkForComputability.begin());
		assert(thisE != 0);

		checkForComputability.erase(thisE);
		if(canComputeEdge(thisE))
		{
			std::vector<double> ScoreAlternatives;
			std::vector<NWEdge*> ScoreAlternatives_bt;

			std::set<NWEdge*> thisE_backward_edges_set = thisE->path->backwardEdges(thisE);
			std::vector<NWEdge*> thisE_backward_edges(thisE_backward_edges_set.begin(), thisE_backward_edges_set.end());

			for(unsigned int bwI = 0; bwI < thisE_backward_edges.size(); bwI++)
			{
				NWEdge* bwE = thisE_backward_edges.at(bwI);
				assert(bwE != 0);
				assert(bwE->scoreComputed);

				double thisEdgeScoreAfterThatEdge = thisE->calculateScore(bwE->isSequenceGap_affine(), bwE->isGraphGap(), S_match, S_mismatch, S_gapOpen, S_gapExtend);
				double edgesCombinedScore = thisEdgeScoreAfterThatEdge + bwE->scoreAfterEdge;
				ScoreAlternatives.push_back(edgesCombinedScore);
				ScoreAlternatives_bt.push_back(bwE);
			}

			assert(ScoreAlternatives.size() > 0);

			std::pair<double, unsigned int> maxScore = Utilities::findVectorMaxP(ScoreAlternatives);

			double selectedScore = maxScore.first;
			NWEdge* selectedBacktrack = ScoreAlternatives_bt.at(maxScore.second);


			thisE->scoreAfterEdge = selectedScore;
			thisE->scoreBacktrack = selectedBacktrack;
			thisE->scoreComputed = true;

			assert(thisE->path != 0);
			std::set<NWEdge*> next_checkForComputability = thisE->path->forwardEdges(thisE);

			checkForComputability.insert(next_checkForComputability.begin(), next_checkForComputability.end());
		}
	}
}


void NWEdge::takeScore_endsFree_nonCritical(double Score, bool previousEdgeAffineSequenceGap, NWEdge* backtrackEdge, double S_match, double S_mismatch, double S_gapOpen, double S_gapExtend, unsigned int* thisSeed)
{
	assert(path != 0);
	std::set<NWEdge*> backwardEdges_set = path->backwardEdges(this);
	std::vector<NWEdge*> backwardEdges(backwardEdges_set.begin(), backwardEdges_set.end());

	std::vector<double> ScoreAlternatives;
	std::vector<NWEdge*> ScoreAlternatives_bt;
	std::vector<bool> ScoreAlternatives_previousEdgeAffineSequenceGap;

	for(unsigned int bwI = 0; bwI < backwardEdges.size(); bwI++)
	{
		NWEdge* bwE = backwardEdges.at(bwI);
		assert(bwE->scoreComputed);
		double thisEdgeScoreAfterThatEdge = calculateScore_endsFree(bwE->isSequenceGap_endsFree(), bwE->isGraphGap(), S_match, S_mismatch, S_gapOpen, S_gapExtend);
		double edgesCombinedScore = thisEdgeScoreAfterThatEdge + bwE->scoreAfterEdge;
		ScoreAlternatives.push_back(edgesCombinedScore);
		ScoreAlternatives_bt.push_back(bwE);
		ScoreAlternatives_previousEdgeAffineSequenceGap.push_back(bwE->isSequenceGap_endsFree());
	}

	ScoreAlternatives.push_back(Score);
	ScoreAlternatives_bt.push_back(backtrackEdge);
	ScoreAlternatives_previousEdgeAffineSequenceGap.push_back(previousEdgeAffineSequenceGap);

	std::pair<double, unsigned int> maxScore = Utilities::findVectorMaxP_nonCritical(ScoreAlternatives, thisSeed);

	double selectedScore = maxScore.first;
	NWEdge* selectedBacktrack = ScoreAlternatives_bt.at(maxScore.second);

	scoreAfterEdge = selectedScore;
	scoreBacktrack = selectedBacktrack;
	scoreComputed = true;
	endsFree_previousEdgeAffineSequenceGap = ScoreAlternatives_previousEdgeAffineSequenceGap.at(maxScore.second);

//	std::cout << "\t halfway takeScore...\n" << std::flush;

	auto canComputeEdge = [&](NWEdge* e) -> bool {
		assert(e != 0);
		if(e->path->entry_edges.count(e))
		{
			return false;
		}
		else
		{
			std::set<NWEdge*> e_backward = e->path->backwardEdges(e);
			bool all_computed = true;
			for(std::set<NWEdge*>::iterator eIt = e_backward.begin(); eIt != e_backward.end(); eIt++)
			{
				NWEdge* previousEdge = *eIt;
				assert(previousEdge != 0);
				if(! previousEdge->scoreComputed)
				{
					all_computed = false;
					break;
				}
			}
			return all_computed;
		}
	};

	std::set<NWEdge*> checkForComputability = path->forwardEdges(this);

	while(checkForComputability.size() > 0)
	{

		NWEdge* thisE = *(checkForComputability.begin());
		assert(thisE != 0);

		checkForComputability.erase(thisE);
		if(canComputeEdge(thisE))
		{
			std::vector<double> ScoreAlternatives;
			std::vector<NWEdge*> ScoreAlternatives_bt;
			std::vector<bool> ScoreAlternatives_previousEdgeAffineSequenceGap;

			std::set<NWEdge*> thisE_backward_edges_set = thisE->path->backwardEdges(thisE);
			std::vector<NWEdge*> thisE_backward_edges(thisE_backward_edges_set.begin(), thisE_backward_edges_set.end());

			for(unsigned int bwI = 0; bwI < thisE_backward_edges.size(); bwI++)
			{
				NWEdge* bwE = thisE_backward_edges.at(bwI);
				assert(bwE != 0);
				assert(bwE->scoreComputed);

				double thisEdgeScoreAfterThatEdge = thisE->calculateScore_endsFree(bwE->isSequenceGap_endsFree(), bwE->isGraphGap(), S_match, S_mismatch, S_gapOpen, S_gapExtend);
				double edgesCombinedScore = thisEdgeScoreAfterThatEdge + bwE->scoreAfterEdge;
				ScoreAlternatives.push_back(edgesCombinedScore);
				ScoreAlternatives_bt.push_back(bwE);
				ScoreAlternatives_previousEdgeAffineSequenceGap.push_back(bwE->isSequenceGap_endsFree());
			}

			assert(ScoreAlternatives.size() > 0);

			std::pair<double, unsigned int> maxScore = Utilities::findVectorMaxP_nonCritical(ScoreAlternatives, thisSeed);

			double selectedScore = maxScore.first;
			NWEdge* selectedBacktrack = ScoreAlternatives_bt.at(maxScore.second);

			thisE->scoreAfterEdge = selectedScore;
			thisE->scoreBacktrack = selectedBacktrack;
			thisE->scoreComputed = true;
			thisE->endsFree_previousEdgeAffineSequenceGap = ScoreAlternatives_previousEdgeAffineSequenceGap.at(maxScore.second);

			assert(thisE->path != 0);
			std::set<NWEdge*> next_checkForComputability = thisE->path->forwardEdges(thisE);

			checkForComputability.insert(next_checkForComputability.begin(), next_checkForComputability.end());
		}
	}
}



void NWEdge::takeScore_endsFree(double Score, bool previousEdgeAffineSequenceGap, NWEdge* backtrackEdge, double S_match, double S_mismatch, double S_gapOpen, double S_gapExtend)
{
	assert(path != 0);
	std::set<NWEdge*> backwardEdges_set = path->backwardEdges(this);
	std::vector<NWEdge*> backwardEdges(backwardEdges_set.begin(), backwardEdges_set.end());

	std::vector<double> ScoreAlternatives;
	std::vector<NWEdge*> ScoreAlternatives_bt;
	std::vector<bool> ScoreAlternatives_previousEdgeAffineSequenceGap;

	for(unsigned int bwI = 0; bwI < backwardEdges.size(); bwI++)
	{
		NWEdge* bwE = backwardEdges.at(bwI);
		assert(bwE->scoreComputed);
		double thisEdgeScoreAfterThatEdge = calculateScore_endsFree(bwE->isSequenceGap_endsFree(), bwE->isGraphGap(), S_match, S_mismatch, S_gapOpen, S_gapExtend);
		double edgesCombinedScore = thisEdgeScoreAfterThatEdge + bwE->scoreAfterEdge;
		ScoreAlternatives.push_back(edgesCombinedScore);
		ScoreAlternatives_bt.push_back(bwE);
		ScoreAlternatives_previousEdgeAffineSequenceGap.push_back(bwE->isSequenceGap_endsFree());
	}

	ScoreAlternatives.push_back(Score);
	ScoreAlternatives_bt.push_back(backtrackEdge);
	ScoreAlternatives_previousEdgeAffineSequenceGap.push_back(previousEdgeAffineSequenceGap);

	std::pair<double, unsigned int> maxScore = Utilities::findVectorMaxP(ScoreAlternatives);

	double selectedScore = maxScore.first;
	NWEdge* selectedBacktrack = ScoreAlternatives_bt.at(maxScore.second);

	scoreAfterEdge = selectedScore;
	scoreBacktrack = selectedBacktrack;
	scoreComputed = true;
	endsFree_previousEdgeAffineSequenceGap = ScoreAlternatives_previousEdgeAffineSequenceGap.at(maxScore.second);

//	std::cout << "\t halfway takeScore...\n" << std::flush;

	auto canComputeEdge = [&](NWEdge* e) -> bool {
		assert(e != 0);
		if(e->path->entry_edges.count(e))
		{
			return false;
		}
		else
		{
			std::set<NWEdge*> e_backward = e->path->backwardEdges(e);
			bool all_computed = true;
			for(std::set<NWEdge*>::iterator eIt = e_backward.begin(); eIt != e_backward.end(); eIt++)
			{
				NWEdge* previousEdge = *eIt;
				assert(previousEdge != 0);
				if(! previousEdge->scoreComputed)
				{
					all_computed = false;
					break;
				}
			}
			return all_computed;
		}
	};

	std::set<NWEdge*> checkForComputability = path->forwardEdges(this);

	while(checkForComputability.size() > 0)
	{

		NWEdge* thisE = *(checkForComputability.begin());
		assert(thisE != 0);

		checkForComputability.erase(thisE);
		if(canComputeEdge(thisE))
		{
			std::vector<double> ScoreAlternatives;
			std::vector<NWEdge*> ScoreAlternatives_bt;
			std::vector<bool> ScoreAlternatives_previousEdgeAffineSequenceGap;

			std::set<NWEdge*> thisE_backward_edges_set = thisE->path->backwardEdges(thisE);
			std::vector<NWEdge*> thisE_backward_edges(thisE_backward_edges_set.begin(), thisE_backward_edges_set.end());

			for(unsigned int bwI = 0; bwI < thisE_backward_edges.size(); bwI++)
			{
				NWEdge* bwE = thisE_backward_edges.at(bwI);
				assert(bwE != 0);
				assert(bwE->scoreComputed);

				double thisEdgeScoreAfterThatEdge = thisE->calculateScore_endsFree(bwE->isSequenceGap_endsFree(), bwE->isGraphGap(), S_match, S_mismatch, S_gapOpen, S_gapExtend);
				double edgesCombinedScore = thisEdgeScoreAfterThatEdge + bwE->scoreAfterEdge;
				ScoreAlternatives.push_back(edgesCombinedScore);
				ScoreAlternatives_bt.push_back(bwE);
				ScoreAlternatives_previousEdgeAffineSequenceGap.push_back(bwE->isSequenceGap_endsFree());
			}

			assert(ScoreAlternatives.size() > 0);

			std::pair<double, unsigned int> maxScore = Utilities::findVectorMaxP(ScoreAlternatives);

			double selectedScore = maxScore.first;
			NWEdge* selectedBacktrack = ScoreAlternatives_bt.at(maxScore.second);

			thisE->scoreAfterEdge = selectedScore;
			thisE->scoreBacktrack = selectedBacktrack;
			thisE->scoreComputed = true;
			thisE->endsFree_previousEdgeAffineSequenceGap = ScoreAlternatives_previousEdgeAffineSequenceGap.at(maxScore.second);

			assert(thisE->path != 0);
			std::set<NWEdge*> next_checkForComputability = thisE->path->forwardEdges(thisE);

			checkForComputability.insert(next_checkForComputability.begin(), next_checkForComputability.end());
		}
	}
}


double NWEdge::calculateScore(bool comingFromAffineSequenceGap, bool comingFromAffineGraphGap, double S_match, double S_mismatch, double S_gapOpen, double S_gapExtend)
{
	assert(!(comingFromAffineGraphGap && comingFromAffineSequenceGap));
	assert(! (isGraphGap() && isSequenceGap()));
	if(isGraphGap())
	{
		if(comingFromAffineGraphGap)
		{
			return S_gapExtend;
		}
		else
		{
			return (S_gapOpen + S_gapExtend);
		}
	}
	else if(isSequenceGap())
	{
		if(isSequenceGap_affine())
		{
			if(comingFromAffineSequenceGap)
			{
				return S_gapExtend;
			}
			else
			{
				return (S_gapOpen + S_gapExtend);
			}
		}
		else
		{
			return 0;
		}
	}
	else
	{
		assert(((to_x - from_x) == 1) && ((to_y - from_y) == 1));
		assert(usedGraphEdge != 0);
		assert(path && path->table && path->table->Context && path->table->Context->sequence && path->table->Context->aligner);

		std::string edgeLabel = usedGraphEdge->getEmission();
		assert(edgeLabel.length() == 1);

		std::string sequenceEmission = path->table->Context->sequence->substr(from_y, 1);

		if(edgeLabel == sequenceEmission)
		{
			return S_match;
		}
		else
		{
			return S_mismatch;
		}
	}

}


double NWEdge::calculateScore_endsFree(bool comingFromAffineSequenceGap, bool comingFromAffineGraphGap, double S_match, double S_mismatch, double S_gapOpen, double S_gapExtend)
{
	assert(!(comingFromAffineGraphGap && comingFromAffineSequenceGap));
	assert(! (isGraphGap() && isSequenceGap()));

	assert(path && path->table && path->table->Context && path->table->Context->sequence && path->table->Context->aligner);
	if(((to_y == 0) && (from_y == 0)) || ((to_y == (int)path->table->Context->sequence->length()) && (from_y == (int)path->table->Context->sequence->length())))
	{
		// return 0;
	}

	if(isGraphGap())
	{
		if(comingFromAffineGraphGap)
		{
			return S_gapExtend;
		}
		else
		{
			return (S_gapOpen + S_gapExtend);
		}
	}
	else if(isSequenceGap())
	{
		assert(usedGraphEdge != 0);
		assert(path && path->table && path->table->Context && path->table->Context->sequence && path->table->Context->aligner);
		std::string edgeLabel = usedGraphEdge->getEmission();
		assert(edgeLabel.length() == 1);

		if(edgeLabel != "_")
		{
			if(comingFromAffineSequenceGap)
			{
				return S_gapExtend;
			}
			else
			{
				return (S_gapOpen + S_gapExtend);
			}
		}
		else
		{
			return 0;
		}
	}
	else
	{
		assert(((to_x - from_x) == 1) && ((to_y - from_y) == 1));
		assert(usedGraphEdge != 0);
		assert(path && path->table && path->table->Context && path->table->Context->sequence && path->table->Context->aligner);

		std::string edgeLabel = usedGraphEdge->getEmission();
		assert(edgeLabel.length() == 1);

		std::string sequenceEmission = path->table->Context->sequence->substr(from_y, 1);

		if(edgeLabel == sequenceEmission)
		{
			return S_match;
		}
		else
		{
			return S_mismatch;
		}
	}

}

bool NWEdge::isSequenceGap()
{
	return
			(((to_x - from_x) == 1) &&
			((to_y - from_y) == 0));
}

bool NWEdge::isSequenceGap_affine()
{
	if(! isSequenceGap())
	{
		return false;
	}

	assert(usedGraphEdge != 0);
	assert(path && path->table && path->table->Context && path->table->Context->sequence && path->table->Context->aligner);
	std::string edgeLabel = usedGraphEdge->getEmission();

	if(edgeLabel == "_")
	{
		return false;
	}
	else
	{
		return true;
	}

}

bool NWEdge::isSequenceGap_endsFree()
{
	if(! isSequenceGap())
	{
		return false;
	}

	assert(usedGraphEdge != 0);
	assert(path && path->table && path->table->Context && path->table->Context->sequence && path->table->Context->aligner);
	std::string edgeLabel = usedGraphEdge->getEmission();

	if(edgeLabel == "_")
	{
		return endsFree_previousEdgeAffineSequenceGap;
	}
	else
	{
		return true;
	}

}


bool NWEdge::isGraphGap()
{
	return
			(((to_x - from_x) == 0) &&
			((to_y - from_y) == 1));
}

}
}

