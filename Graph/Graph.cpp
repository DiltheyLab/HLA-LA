/*
 * Graph.cpp
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#include "Graph.h"
#include "../Utilities.h"

#include <assert.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <stdlib.h>
#include <time.h>
#include <exception>
#include <stdexcept>

#include "HaplotypePanel.h"
typedef std::basic_string <unsigned char> ustring;

std::string separatorForSerialization = "|||";

class positionsSorter {
public:
	map<string, int>* p;
	bool operator() (string i, string j) {
		if(p->count(i) == 0)
			throw std::runtime_error("No position information for locus: "+i);
		if(p->count(j) == 0)
			throw std::runtime_error("No position information for locus: "+j);
		return ((*p)[i] < (*p)[j]);
	}
};

void Graph::checkAlignmentBackbonePresences(std::vector<std::pair<std::string, std::vector<int>>> backbones)
{
	std::cout << "Graph::checkAlignmentBackbonePresences(..):\n";
	for(unsigned int i = 0; i < backbones.size(); i++)
	{
		std::cout << "\t" << i << "\t" << backbones.at(i).first << "\n";
		std::set<unsigned int> occurrences = checkPartialAlignmentBackbonePresence(backbones.at(i), false);
		if(occurrences.size() == 0)
		{
			std::cout << "\t\tNOT FOUND\n";
		}
		else
		{
			for(auto o : occurrences)
			{
				std::cout << "\t\t" << o << "\n";
			}
		}
		std::cout << std::flush;
	}
}

void Graph::checkAlignmentBackbonePresences_ignoreGraphGaps(std::vector<std::pair<std::string, std::vector<int>>> backbones)
{
	std::cout << "Graph::checkAlignmentBackbonePresences_ignoreGraphGaps(..):\n" << std::flush;
	for(unsigned int i = 0; i < backbones.size(); i++)
	{
		std::cout << "\t" << i << "\t" << backbones.at(i).first << "\n" << std::flush;
		std::set<unsigned int> occurrences = checkPartialAlignmentBackbonePresence_ignoreGraphGaps(backbones.at(i), false);
		if(occurrences.size() == 0)
		{
			std::cout << "\t\tNOT FOUND\n";
		}
		else
		{
			for(auto o : occurrences)
			{
				std::cout << "\t\t" << o << "\n";
			}
		}
		std::cout << std::flush;
	}
}

std::set<unsigned int> Graph::checkPartialAlignmentBackbonePresence(std::pair<std::string, std::vector<int>> alignmentBackbone, bool greedy, bool verbose)
{
	assert(alignmentBackbone.first.length() == alignmentBackbone.second.size());

	std::string alignmentSequence_noRealGaps;
	alignmentSequence_noRealGaps.reserve(alignmentBackbone.first.length());

	for(unsigned int i = 0; i < alignmentBackbone.first.length(); i++)
	{
		if(alignmentBackbone.second.at(i) != -1)
		{
			alignmentSequence_noRealGaps.push_back(alignmentBackbone.first.at(i));
		}
	}

	return checkPartialSequencePresence(alignmentSequence_noRealGaps, greedy, verbose);
}

std::set<unsigned int> Graph::checkPartialAlignmentBackbonePresence_ignoreGraphGaps(std::pair<std::string, std::vector<int>> alignmentBackbone, bool greedy, bool verbose)
{
	assert(alignmentBackbone.first.length() == alignmentBackbone.second.size());

	std::string alignmentSequence_noGaps;
	alignmentSequence_noGaps.reserve(alignmentBackbone.first.length());

	for(unsigned int i = 0; i < alignmentBackbone.first.length(); i++)
	{
		if(alignmentBackbone.first.at(i) != '_')
		{
			alignmentSequence_noGaps.push_back(alignmentBackbone.first.at(i));
		}
	}

	return checkPartialSequencePresence_ignoreGraphGaps(alignmentSequence_noGaps, greedy, verbose);
}

std::set<unsigned int> Graph::checkPartialSequencePresence(std::string S, bool greedy, bool verbose)
{
	std::set<unsigned int> found_positions;
	for(unsigned int lI = 0; lI < NodesPerLevel.size(); lI++)
	{
		bool found = checkSequencePresence(S, lI, false);
		if(found)
		{
			found_positions.insert(lI);
			if(greedy)
			{
				break;
			}
		}
	}

	return found_positions;
}

std::set<unsigned int> Graph::checkPartialSequencePresence_ignoreGraphGaps(std::string S, bool greedy, bool verbose)
{
	std::set<unsigned int> found_positions;
	for(unsigned int lI = 0; lI < NodesPerLevel.size(); lI++)
	{
		if((lI % 10000) == 0)
		{
			// std::cout << "\t\tScanning " << lI << "\n" << std::flush;
		}
		bool found = checkSequencePresence_ignoreGraphGaps(S, lI, false);
		if(found)
		{
			found_positions.insert(lI);
			if(greedy)
			{
				break;
			}
		}
	}

	return found_positions;
}

bool Graph::checkSequencePresence(std::string S, int startingAt, bool verbose)
{
	if((startingAt + S.length() - 1) > (NodesPerLevel.size() - 1))
	{
		return false;
	}

	std::set<Node*> runningNodes = NodesPerLevel.at(startingAt);
	for(unsigned int lI = 0; lI < S.length(); lI++)
	{
		if(verbose)
			std::cout << "Level " << lI << " runningNodes.size(): " << runningNodes.size() << "\n" << std::flush;

		if(runningNodes.size() == 0)
		{
			return false;
		}
		std::set<Node*> nodesNextLevel;
		std::string C = S.substr(lI, 1);

		if(verbose)
			std::cout << "\tLooking for " << C << "\n" << std::flush;

		for(std::set<Node*>::iterator nIt = runningNodes.begin(); nIt != runningNodes.end(); nIt++)
		{
			Node* n = *nIt;
			for(std::set<Edge*>::iterator eIt = n->Outgoing_Edges.begin(); eIt != n->Outgoing_Edges.end(); eIt++)
			{
				Edge* e = *eIt;

				if(verbose)
					std::cout << "\t\tEdge " << e << " emission " << e->getEmission() << "\n" << std::flush;
				if(e->getEmission() == C)
				{
					nodesNextLevel.insert(e->To);
				}
			}
		}
		runningNodes = nodesNextLevel;
	}

	if(runningNodes.size())
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool Graph::checkSequencePresence_ignoreGraphGaps(std::string S, int startingAt, bool verbose)
{
	if((startingAt + S.length() - 1) > (NodesPerLevel.size() - 1))
	{
		return false;
	}

	std::vector<std::pair<Node*, int>> runningPositions;
	for(auto initialNode : NodesPerLevel.at(startingAt))
	{
		runningPositions.push_back(make_pair(initialNode, 0));
	}

//	int graphLevel = startingAt;
	while(runningPositions.size())
	{
//		if(verbose)
//			std::cout << "Graph level " << graphLevel << " runningPositions.size(): " << runningPositions.size() << "\n" << std::flush;

		std::vector<std::pair<Node*, int>> runningPositionsNextLevel;

		for(std::vector<std::pair<Node*, int>>::iterator nIt = runningPositions.begin(); nIt != runningPositions.end(); nIt++)
		{
			Node* n = nIt->first;
			int positionInS = nIt->second;

//			assert((int)n->level == graphLevel);

			assert(positionInS < (int)S.length());
			std::string C = S.substr(positionInS, 1);
			assert(C != "_");

			if(verbose)
				std::cout << "\tLooking for " << C << "\n" << std::flush;

			std::set<Node*> nodes_oneCharacterConsumed;
			std::set<Node*> nodes_gaps;

			if(gapEdgePaths_connectedNodes_forwards.count(n))
			{
				for(auto afterGapNode : gapEdgePaths_connectedNodes_forwards.at(n))
				{
					nodes_gaps.insert(afterGapNode.first);
				}
			}

			for(std::set<Edge*>::iterator eIt = n->Outgoing_Edges.begin(); eIt != n->Outgoing_Edges.end(); eIt++)
			{
				Edge* e = *eIt;

				if(verbose)
					std::cout << "\t\tEdge " << e << " emission " << e->getEmission() << "\n" << std::flush;

				if(e->getEmission() == C)
				{
					Node* toNode = e->To;
					nodes_oneCharacterConsumed.insert(toNode);
				}
			}

			for(auto toNode : nodes_oneCharacterConsumed)
			{
				if((positionInS + 1) == (int)S.length())
				{
					return true;
				}
				runningPositionsNextLevel.push_back(make_pair(toNode, positionInS + 1));
			}

			for(auto toNode : nodes_gaps)
			{
				runningPositionsNextLevel.push_back(make_pair(toNode, positionInS));
			}
		}
		runningPositions = runningPositionsNextLevel;
//		graphLevel++;
	}

	return false;
}


bool Graph::checkSequencePresence(std::string S, bool verbose)
{
	if(S.length() > (NodesPerLevel.size() - 1))
	{
		return false;
	}
	else
	{
		std::set<Node*> runningNodes = NodesPerLevel.at(0);
		for(unsigned int lI = 0; lI < S.length(); lI++)
		{
			if(verbose)
				std::cout << "Level " << lI << " runningNodes.size(): " << runningNodes.size() << "\n" << std::flush;
			
			if(runningNodes.size() == 0)
			{
				return false;
			}
			std::set<Node*> nodesNextLevel;
			std::string C = S.substr(lI, 1);
			
			if(verbose)
				std::cout << "\tLooking for " << C << "\n" << std::flush;
			for(std::set<Node*>::iterator nIt = runningNodes.begin(); nIt != runningNodes.end(); nIt++)
			{
				Node* n = *nIt;
				for(std::set<Edge*>::iterator eIt = n->Outgoing_Edges.begin(); eIt != n->Outgoing_Edges.end(); eIt++)
				{
					Edge* e = *eIt;
					
					if(verbose)
						std::cout << "\t\tEdge " << e << " emission " << e->getEmission() << "\n" << std::flush;
					if(e->getEmission() == C)
					{
						nodesNextLevel.insert(e->To);
					}
				}
			}
			runningNodes = nodesNextLevel;
		}

		if(runningNodes.size())
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}

void Graph::computeGapEdgePaths()
{
	if(! haveComputedGapEdgePaths)
	{
		Node* firstPseudoNode = new Node();
		Node* lastPseudoNode = new Node();
		pseudoNodes.insert(firstPseudoNode);
		pseudoNodes.insert(lastPseudoNode);

		haveComputedGapEdgePaths = true;

		// find graph gap edge paths

		std::map<Node*, std::map<Node*, std::vector<Edge*>>> runningPaths;
		std::cerr << Utilities::timestamp() << "Start graph gap path analysis.\n" << std::flush;

		for(unsigned int lI = 0; lI < NodesPerLevel.size(); lI++)
		{
			if((lI % 1000) == 0)
			{
				std::cerr << "\t" << lI << "/" << NodesPerLevel.size() << "\t" << "runningPaths.size(): " << runningPaths.size() << "\t" << "completedGapEdgePaths.size(): " << completedGapEdgePaths.size() << "\n" << std::flush;
			}
			std::map<Node*, std::map<Node*, std::vector<Edge*>>> runningPaths_nextLevel;
			std::set<Edge*> seen_gap_edge;
			for(std::map<Node*, std::map<Node*, std::vector<Edge*>>>::iterator nodeIt = runningPaths.begin(); nodeIt != runningPaths.end(); nodeIt++)
			{
				Node* thisLevelNode = nodeIt->first;
				assert(thisLevelNode->level == lI);

				std::set<Edge*> egdges = thisLevelNode->Outgoing_Edges;
				int non_gap_edges = 0;
				for(std::set<Edge*>::iterator edgeIt = egdges.begin(); edgeIt != egdges.end(); edgeIt++)
				{
					Edge* e = *edgeIt;
					if(e->getEmission() == "_")
					{
						seen_gap_edge.insert(e);
						Node* targetNode = e->To;
						for(std::map<Node*, std::vector<Edge*>>::iterator fromNodeIt = runningPaths.at(thisLevelNode).begin(); fromNodeIt != runningPaths.at(thisLevelNode).end(); fromNodeIt++)
						{
							Node* fromNode = fromNodeIt->first;

							if( (runningPaths_nextLevel.count(targetNode) == 0) ||
								(runningPaths_nextLevel.at(targetNode).count(fromNode) == 0))
							{
								runningPaths_nextLevel[targetNode][fromNode] = fromNodeIt->second;
								runningPaths_nextLevel.at(targetNode).at(fromNode).push_back(e);
							}
						}
					}
					else
					{
						non_gap_edges++;
					}
				}

				if((non_gap_edges != 0) || (lI == (NodesPerLevel.size() - 1)))
				{
					for(std::map<Node*, std::vector<Edge*>>::iterator fromNodeIt = runningPaths.at(thisLevelNode).begin(); fromNodeIt != runningPaths.at(thisLevelNode).end(); fromNodeIt++)
					{
						Node* fromNode = fromNodeIt->first;
						Node* firstNode = fromNodeIt->second.front()->From;
						assert(firstNode == fromNode);

						Node* lastNode = fromNodeIt->second.back()->To;
						assert(lastNode == thisLevelNode);

						for(unsigned int i = 0; i < fromNodeIt->second.size(); i++)
						{
							Edge* e = fromNodeIt->second.at(i);
							assert(e->getEmission() == "_");
							if(i < (fromNodeIt->second.size()-1))
							{
								assert(fromNodeIt->second.at(i)->To == fromNodeIt->second.at(i+1)->From);
							}
	//						std::cerr << "\t" << e << "\n";
						}
						completedGapEdgePaths.push_back(fromNodeIt->second);
	//					std::cerr << "Completed path length " << fromNodeIt->second.size() << ", " << fromNode << ", level " << fromNode->level << " to " << thisLevelNode << ", level " << thisLevelNode->level << "\n" << std::flush;
					}
				}
			}

			std::set<Edge*> edges_thisLevel = getEdgesEmanatingFromLevel(lI);
			for(std::set<Edge*>::iterator edgeIt = edges_thisLevel.begin(); edgeIt != edges_thisLevel.end(); edgeIt++)
			{
				Edge* e = *edgeIt;
				if(e->getEmission() == "_")
				{
					if(seen_gap_edge.count(e) == 0)
					{
						Node* fromNode = e->From;
						assert(fromNode->level == lI);
						Node* targetNode = e->To;
						assert(targetNode->level == (lI + 1));

						if( (runningPaths_nextLevel.count(targetNode) == 0) ||
							(runningPaths_nextLevel.at(targetNode).count(fromNode) == 0))
						{
							std::vector<Edge*> newEdgePath;
							newEdgePath.push_back(e);
							runningPaths_nextLevel[targetNode][fromNode] = newEdgePath;
						}
					}
				}
			}

			runningPaths = runningPaths_nextLevel;
		}

		std::cerr << Utilities::timestamp() << "extensionAligner::extensionAligner(..): Indexing of edge paths; have " << completedGapEdgePaths.size() << " completed edge paths.\n" << std::flush;

		for(unsigned int pathI = 0; pathI < completedGapEdgePaths.size(); pathI++)
		{
			Edge* pseudoEdge = new Edge();
			pseudoEdges_correspondingToGapEdgePaths[pseudoEdge] = pathI;
			Node* firstNode = completedGapEdgePaths.at(pathI).front()->From;
			Node* lastNode = completedGapEdgePaths.at(pathI).back()->To;

			assert((gapEdgePaths_connectedNodes_forwards.count(firstNode) == 0) || (gapEdgePaths_connectedNodes_forwards.at(firstNode).count(lastNode) == 0));
			assert((gapEdgePaths_connectedNodes_backwards.count(lastNode) == 0) || (gapEdgePaths_connectedNodes_backwards.at(lastNode).count(firstNode) == 0));

			gapEdgePaths_connectedNodes_forwards[firstNode][lastNode] = pseudoEdge;
			gapEdgePaths_connectedNodes_backwards[lastNode][firstNode] = pseudoEdge;

			pseudoEdge->From = firstNode;
			pseudoEdge->To = lastNode;
		}
	}
}

Graph::~Graph() {
//	for(std::map<Edge*, unsigned int>::iterator pseudoEdgeIt = pseudoEdges_correspondingToGapEdgePaths.begin(); pseudoEdgeIt != pseudoEdges_correspondingToGapEdgePaths.end(); pseudoEdgeIt++)
//	{
//		Edge* pseudoEdge = pseudoEdgeIt->first;
//		delete(pseudoEdge);
//	}
//	for(std::set<Node*>::iterator nodeIt = pseudoNodes.begin(); nodeIt != pseudoNodes.end(); nodeIt++)
//	{
//		Node* n = *nodeIt;
//		delete(n);
//	}
}

Graph::Graph() {
	haveComputedGapEdgePaths = false;
}

void Graph::regenerateNodeIncomingOutgoingEdges()
{
	// for after de-serialization
	for(std::set<Edge*>::iterator edgeIt = Edges.begin(); edgeIt != Edges.end(); edgeIt++)
	{
		Edge* e = *edgeIt;
		Node* fromNode = e->From;
		Node* toNode = e->To;
		assert(Nodes.count(fromNode));
		assert(Nodes.count(toNode));
		fromNode->Outgoing_Edges.insert(e);
		toNode->Incoming_Edges.insert(e);
	}
	
	for(std::set<Node*>::iterator nodeIt = Nodes.begin(); nodeIt != Nodes.end(); nodeIt++)
	{
		Node* n = *nodeIt;
		n->g = this;
	}
	
	checkStructure();
}
void Graph::checkStructure()
{
	std::set<Node*> nodesInLevels;
	for(unsigned int lI = 0; lI < NodesPerLevel.size(); lI++)
	{
		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(lI).begin(); nodeIt != NodesPerLevel.at(lI).end(); nodeIt++)
		{
			Node* n = *nodeIt;
			assert(Nodes.count(n));
			nodesInLevels.insert(n);
			assert(n->g == this);
		}
	}
	
	assert(Nodes.size() == nodesInLevels.size());
	
	std::set<Edge*> edgesAtNodes;
	for(set<Node*>::iterator nodeIt = Nodes.begin(); nodeIt != Nodes.end(); nodeIt++)
	{
		Node* n = *nodeIt;
		for(std::set<Edge*>::iterator edgeIt = n->Outgoing_Edges.begin(); edgeIt != n->Outgoing_Edges.end(); edgeIt++)
		{
			Edge* e = *edgeIt;
			assert(e != 0);
			assert(Edges.count(e));
			Node* fromNode = e->From;
			Node* toNode = e->To;
			assert(Nodes.count(fromNode));
			assert(Nodes.count(toNode));
			edgesAtNodes.insert(e);
			
			assert(fromNode->Outgoing_Edges.count(e));
			assert(toNode->Incoming_Edges.count(e));
		}
	}
	
	assert(edgesAtNodes.size() == Edges.size());
}

std::set<Edge*> Graph::getEdgesEmanatingFromLevel(unsigned int lI)
{
	std::set<Edge*> forReturn;
	for(std::set<Node*>::iterator nodeIt = NodesPerLevel.at(lI).begin(); nodeIt != NodesPerLevel.at(lI).end(); nodeIt++)
	{
		Node* n = *nodeIt;
		forReturn.insert(n->Outgoing_Edges.begin(), n->Outgoing_Edges.end());
	}
	return forReturn;
}

void Graph::buildFromHaplotypes(HaplotypePanel hp, bool wantPGFprotection, int want_suffix_length)
{
	Node* n0 = new Node();
	n0->level = 0;
	n0->terminal = false;
	registerNode(n0, 0);

	vector<int> realHaplotypeIndices;
	map<int, set<int> > snpHaplotypeIndices;

	map<Node*, set<int> > NodeToHaplotype;

	int last_real_haplotype_index = -1;
	int pgf_haplotype_id = -1;
	for(unsigned int haplotypeI = 0; haplotypeI < hp.HaplotypeIDs.size(); haplotypeI++)
	{
		string haplotypeID = hp.HaplotypeIDs.at(haplotypeI);
		if(haplotypeID.substr(0,4) == "SNPs")
		{
			assert(last_real_haplotype_index != -1);
			snpHaplotypeIndices[last_real_haplotype_index].insert(haplotypeI);
		}
		else
		{
			realHaplotypeIndices.push_back(haplotypeI);
			NodeToHaplotype[n0].insert(haplotypeI);
			last_real_haplotype_index = haplotypeI;
		}

		if((haplotypeID == "pgf_____1") || ((haplotypeID.length() > 10) && (haplotypeID.substr(0, 10) == "pgfallele_")))
		{
			pgf_haplotype_id = haplotypeI;
			cout << "PGF haplotype index: " << pgf_haplotype_id << "\n";
		}
	}

	if(wantPGFprotection)
	{
		if(pgf_haplotype_id == -1)
		{
			cout << "Haplotype IDs:\n";
			for(unsigned int haplotypeI = 0; haplotypeI < hp.HaplotypeIDs.size(); haplotypeI++)
			{
				string haplotypeID = hp.HaplotypeIDs.at(haplotypeI);
				cout << haplotypeID << "\n";
			}
		}
		assert(pgf_haplotype_id != -1);
	}

	vector<string> loci = hp.getOrderedLoci();
	positionsSorter sortClass;
	sortClass.p = &(hp.LocusPositions);
	sort(loci.begin(), loci.end(), sortClass);

	string lastLocus = loci.at(loci.size()-1);
	int lastLocusPosition = hp.LocusPositions[lastLocus];

	string pufferLocusID = "END_PUFFER";
	int pufferLocusLocation = lastLocusPosition+1;
	loci.push_back(pufferLocusID);
	hp.LocusPositions[pufferLocusID] = pufferLocusLocation;
	hp.LocusStrands[pufferLocusID] = "+";
	unsigned char puffer_symbol = 'N';
	for(unsigned int haplotypeI = 0; haplotypeI < hp.HaplotypeIDs.size(); haplotypeI++)
	{
		hp.HaplotypesByLoci[pufferLocusID].push_back(puffer_symbol);
	}

	for(int level = 0; level < (int)loci.size(); level++)
	{

		if(level > 5)
		{
			//exit(0);
		}

		if((level % 100000) == 0)
		{
		}

		if((level % 1000) == 0)
		{
			// cout << "\rLevel " << level << flush;
		}

		string locusID = loci.at(level);
		unsigned char star_symbol = '*';

		// cout << "\n";
		// cout << "\n\nLevel " << level << " locusID " << locusID << flush;

		map<Edge*, set<int> > looseEdges;
		map<Edge*, set<Edge*> > snpEdges;
		set<Edge*> snpEdges_to_catch;

		set<int> seen_haplotypes;

		bool level_protected_PGF = false;
		Edge* pgfEdge = 0;


		for(map< Node*, set<int> >::iterator nIt = NodeToHaplotype.begin(); nIt != NodeToHaplotype.end(); nIt++)
		{
			// cout << "\tNode " << nIt->first << "\n";

			map<unsigned char, set<int> > emission2Haplotypes;
			map<unsigned char, set<unsigned char> > emissions_possible_SNPs;

			for(set<int>::iterator hIt = nIt->second.begin(); hIt != nIt->second.end(); hIt++)
			{
				int haplotype_index = *hIt;

				assert(hp.HaplotypesByLoci.count(locusID) > 0);

				unsigned char haplo_symbol = hp.HaplotypesByLoci[locusID].at(haplotype_index);
				emission2Haplotypes[haplo_symbol].insert(haplotype_index);

				for(set<int>::iterator snpHaplotypesIndicesIt = snpHaplotypeIndices[haplotype_index].begin(); snpHaplotypesIndicesIt != snpHaplotypeIndices[haplotype_index].end(); snpHaplotypesIndicesIt++)
				{
					unsigned char haplo_symbol_possibleSNP = hp.HaplotypesByLoci[locusID].at(*snpHaplotypesIndicesIt);

					if(haplo_symbol_possibleSNP != star_symbol)
					{
						emissions_possible_SNPs[haplo_symbol].insert(haplo_symbol_possibleSNP);
						assert(haplo_symbol != star_symbol);
					}
				}

				seen_haplotypes.insert(haplotype_index);

				// cout << "\t\t attached haplotype  " << hp.HaplotypeIDs.at(haplotype_index) << "char: " << haplo_symbol << " real symbol: " << g->CODE.deCode(locusID, haplo_symbol) << "\n";
			}

			/*
			if(emission2Haplotypes.count(star_symbol) > 0)
			{
				if(emission2Haplotypes.size() == 2)
				{
					set<int> haplosWithStar = emission2Haplotypes[star_symbol];
					emission2Haplotypes.erase(star_symbol);
					unsigned char remainingNonStar = (emission2Haplotypes.begin()->first);
					emission2Haplotypes[remainingNonStar].insert(haplosWithStar.begin(), haplosWithStar.end());
				}
			}
			*/

			for(map<unsigned char, set<int> >::iterator charIt = emission2Haplotypes.begin(); charIt != emission2Haplotypes.end(); charIt++)
			{
				assert(charIt->second.size() > 0);
				unsigned char symbol = charIt->first;

				//cout << "level " << level << " node " << nIt->first << " char " << symbol << " value " << g->CODE.deCode(locusID, symbol) << "\n";
				Edge* newE = new Edge();

				newE->count = 1;
				newE->emission = symbol;
				newE->locus_id = locusID;
				if((locusID.substr(0,2) == "rs") || (locusID.substr(0,3) == "HLA"))
				{
					newE->label = locusID+"*";
					newE->label.push_back(symbol);
				}

				newE->pgf_protect = false;

				if((pgf_haplotype_id != -1) && (charIt->second.count(pgf_haplotype_id) > 0))
				{
					newE->pgf_protect = true;
					level_protected_PGF = true;

					/*
					if(level > 0)
					{
						bool found_preceding_PGF = false;
						Node* fromNode = nIt->first;
						for(set<Edge*>::iterator InEdgeIt = fromNode->Incoming_Edges.begin(); InEdgeIt != fromNode->Incoming_Edges.end(); InEdgeIt++)
						{
							Edge *preEdge = *InEdgeIt;
							if((bool)preEdge->pgf_protect)
							{
								found_preceding_PGF = true;
							}
						}
						assert(found_preceding_PGF);
					}
					*/

					assert(pgfEdge == 0);
					pgfEdge = newE;
				}

				registerEdge(newE);
				newE->From = nIt->first;
				nIt->first->Outgoing_Edges.insert(newE);

				looseEdges[newE] = charIt->second;

				if(emissions_possible_SNPs.count(symbol) > 0)
				{
					assert(symbol != star_symbol);
					for(set<unsigned char>::iterator possibleSNPit = emissions_possible_SNPs[symbol].begin(); possibleSNPit != emissions_possible_SNPs[symbol].end(); possibleSNPit++ )
					{
						unsigned char alternativeEdgeSymbol = *possibleSNPit;
						if(alternativeEdgeSymbol != symbol)
						{
							Edge* newE_SNP = new Edge();

							newE_SNP->count = 1;
							newE_SNP->emission = alternativeEdgeSymbol;
							newE_SNP->locus_id = locusID;
							newE_SNP->pgf_protect = false;

							if((locusID.substr(0,2) == "rs") || (locusID.substr(0,3) == "HLA"))
							{
								newE_SNP->label = locusID+"*";
								newE->label.push_back(alternativeEdgeSymbol);
							}
							registerEdge(newE_SNP);
							newE_SNP->From = nIt->first;
							nIt->first->Outgoing_Edges.insert(newE_SNP);

							snpEdges[newE].insert(newE_SNP);

							snpEdges_to_catch.insert(newE_SNP);
						}
					}
				}
			}


		}

		// PGF protection
		assert((! wantPGFprotection) || level_protected_PGF);

		vector< set<int> > groupingHaplotypes;
		vector< set<Edge*> > groupingEdges;

		for(map<Edge*, set<int> >::iterator looseEdgeIt = looseEdges.begin(); looseEdgeIt != looseEdges.end(); looseEdgeIt++)
		{
			assert(looseEdgeIt->second.size() > 0);
			groupingHaplotypes.push_back(looseEdgeIt->second);
			set<Edge*> edgeSet;
			edgeSet.insert(looseEdgeIt->first);
			groupingEdges.push_back(edgeSet);
		}

		vector< vector<int> > pairsToCheck;
		for(int i = 0; i < (int)groupingHaplotypes.size(); i++)
		{
			for(int j = i+1; j < (int)groupingHaplotypes.size(); j++)
			{
				vector<int> p;
				p.push_back(i);
				p.push_back(j);
				pairsToCheck.push_back(p);
			}
		}
		set<int> groupsDeleted;

		while(pairsToCheck.size() != 0)
		{
			vector<int> thisPair = pairsToCheck.at(0);
			assert(thisPair.size() == 2);

			pairsToCheck.erase(pairsToCheck.begin());

			int g1 = thisPair.at(0);
			int g2 = thisPair.at(1);
			if((groupsDeleted.count(g1) > 0) || (groupsDeleted.count(g2) > 0))
			{
				continue;
			}
			else
			{

				bool joinPair = true;

				int local_need_suffix_length = want_suffix_length;
				int got_suffix_length = 0;

				while((joinPair == true) && (want_suffix_length != got_suffix_length))
				{
					//cout << "want_suffix_length: " << want_suffix_length << " local_need_suffix_length: " << local_need_suffix_length << " got_suffix_length: " << got_suffix_length << "\n";

					map< basic_string<unsigned char>, int> suffixes_g1;
					map< basic_string<unsigned char>, int> suffixes_g2;
					set< basic_string<unsigned char > > suffixes;

					if(local_need_suffix_length > (want_suffix_length*10))
					{
						joinPair = false;
						continue;
					}

					if((level + local_need_suffix_length) > (int)(loci.size() - 1))
					{
						joinPair = false;
					}
					else
					{
						vector<string> takeSuffixLoci;

						for(int lI = 1; lI <= local_need_suffix_length; lI++)
						{
							int suffixL = level + lI;
							takeSuffixLoci.push_back(loci.at(suffixL));
						}

						vector<unsigned char> takeSuffixLociStars;
						for(vector<string>::iterator locusIt = takeSuffixLoci.begin(); locusIt != takeSuffixLoci.end(); locusIt++)
						{
							string locus = *locusIt;
							unsigned char starSymbol = '*';
							takeSuffixLociStars.push_back(starSymbol);
						}

						for(set<int>::iterator hIt = groupingHaplotypes.at(g1).begin(); hIt != groupingHaplotypes.at(g1).end(); hIt++)
						{
							// normal haplotype

							basic_string<unsigned char> suffix = hp.getIndividualHaplotype(takeSuffixLoci, hp.HaplotypeIDs.at(*hIt), false);
							assert((int)suffix.length() == local_need_suffix_length);
							if(suffixes_g1.count(suffix) == 0)
							{
								suffixes_g1[suffix] = 0;
							}
							suffixes_g1[suffix]++;
							suffixes.insert(suffix);

							// SNPs

							for(set<int>::iterator snpHaplotypeIt = snpHaplotypeIndices[*hIt].begin(); snpHaplotypeIt != snpHaplotypeIndices[*hIt].end(); snpHaplotypeIt++)
							{
								basic_string<unsigned char> suffix_2 = hp.getIndividualHaplotype(takeSuffixLoci, hp.HaplotypeIDs.at(*snpHaplotypeIt), false);
								assert((int)suffix_2.length() == local_need_suffix_length);
								for(int pos = 0; pos < (int)suffix_2.length(); pos++)
								{
									if(suffix_2.at(pos) == takeSuffixLociStars.at(pos))
									{
										suffix_2.at(pos) = suffix.at(pos);
									}
								}
								if(suffixes_g1.count(suffix_2) == 0)
								{
									suffixes_g1[suffix_2] = 0;
								}
								suffixes_g1[suffix_2]++;
								suffixes.insert(suffix_2);
							}

						}

						for(set<int>::iterator hIt = groupingHaplotypes.at(g2).begin(); hIt != groupingHaplotypes.at(g2).end(); hIt++)
						{
							// normal haplotypes

							basic_string<unsigned char> suffix = hp.getIndividualHaplotype(takeSuffixLoci, hp.HaplotypeIDs.at(*hIt), false);
							assert((int)suffix.length() == local_need_suffix_length);

							if(suffixes_g2.count(suffix) == 0)
							{
								suffixes_g2[suffix] = 0;
							}
							suffixes_g2[suffix]++;
							suffixes.insert(suffix);

							// SNPs

							for(set<int>::iterator snpHaplotypeIt = snpHaplotypeIndices[*hIt].begin(); snpHaplotypeIt != snpHaplotypeIndices[*hIt].end(); snpHaplotypeIt++)
							{
								basic_string<unsigned char> suffix_2 = hp.getIndividualHaplotype(takeSuffixLoci, hp.HaplotypeIDs.at(*snpHaplotypeIt), false);
								assert((int)suffix_2.length() == local_need_suffix_length);

								for(int pos = 0; pos < (int)suffix_2.length(); pos++)
								{
									if(suffix_2.at(pos) == takeSuffixLociStars.at(pos))
									{
										suffix_2.at(pos) = suffix.at(pos);
									}
								}

								if(suffixes_g2.count(suffix_2) == 0)
								{
									suffixes_g2[suffix_2] = 0;
								}
								suffixes_g2[suffix_2]++;
								suffixes.insert(suffix_2);
							}
						}

						got_suffix_length =  local_need_suffix_length;


						bool one_node_all_stars = false;
						for(set< basic_string<unsigned char> >::iterator suffixIt = suffixes.begin(); suffixIt != suffixes.end(); suffixIt++)
						{
							assert((int)suffixIt->size() == local_need_suffix_length);
							bool all_stars = true;
							for(int i = 0; i < local_need_suffix_length; i++)
							{
								unsigned char suffixC = suffixIt->at(i);
								string suffixCharLocus = takeSuffixLoci.at(i);
								string decodedSuffix;
								decodedSuffix.push_back(suffixC);
								all_stars = (all_stars && (decodedSuffix == "*"));
							}

							if(all_stars)
							{
								if(((suffixes_g1.count(*suffixIt) > 0) && (suffixes_g1.size() == 1)) || ((suffixes_g2.count(*suffixIt) > 0) && (suffixes_g2.size() == 1)))
								{
									one_node_all_stars = true;
								}
							}
						}


						for(set< basic_string<unsigned char> >::iterator suffixIt = suffixes.begin(); suffixIt != suffixes.end(); suffixIt++)
						{
							int gapCounter = 0;
							assert((int)suffixIt->size() == local_need_suffix_length);
							bool all_stars = true;
							for(int i = 0; i < local_need_suffix_length; i++)
							{
								unsigned char suffixC = suffixIt->at(i);
								string suffixCharLocus = takeSuffixLoci.at(i);
								string decodedSuffix;
								decodedSuffix.push_back(suffixC);
								all_stars = (all_stars && (decodedSuffix == "*"));

								if(decodedSuffix == "_")
								{
									gapCounter++;
									if(i == 0)
									{
										joinPair = false;
									}
								}
							}

							if((local_need_suffix_length - gapCounter) < got_suffix_length)
							{
								got_suffix_length = (local_need_suffix_length - gapCounter);
							}

							if(! one_node_all_stars)
							{
								if(! all_stars)
								{
									if(!((suffixes_g1.count(*suffixIt) > 0) && (suffixes_g2.count(*suffixIt) > 0)))
									{
										joinPair = false;
									}
								}
							}
						}
						assert(got_suffix_length <= local_need_suffix_length);
						local_need_suffix_length = want_suffix_length + (local_need_suffix_length - got_suffix_length);
					}
				}

				if(joinPair)
				{
					set<int> newGroupHaplos;
					set<Edge*> newGroupEdges;

					newGroupHaplos.insert(groupingHaplotypes.at(g1).begin(), groupingHaplotypes.at(g1).end());
					newGroupHaplos.insert(groupingHaplotypes.at(g2).begin(), groupingHaplotypes.at(g2).end());
					groupingHaplotypes.at(g1).clear();
					groupingHaplotypes.at(g2).clear();

					newGroupEdges.insert(groupingEdges.at(g1).begin(), groupingEdges.at(g1).end());
					newGroupEdges.insert(groupingEdges.at(g2).begin(), groupingEdges.at(g2).end());
					groupingEdges.at(g1).clear();
					groupingEdges.at(g2).clear();

					groupingHaplotypes.push_back(newGroupHaplos);
					groupingEdges.push_back(newGroupEdges);
					assert(groupingHaplotypes.size() == groupingEdges.size());

					groupsDeleted.insert(g1);
					groupsDeleted.insert(g2);

					int index_of_new_element = groupingHaplotypes.size()-1;
					assert(groupingHaplotypes.at(index_of_new_element).size() > 0);
					for(int i = 0; i < index_of_new_element; i++)
					{
						if(groupsDeleted.count(i) > 0)
						{
							continue;
						}
						else
						{
							vector<int> newP;
							newP.push_back(i);
							newP.push_back(index_of_new_element);
							pairsToCheck.push_back(newP);
						}
					}
				}
			}
		}

		for(unsigned int haplotypeII = 0; haplotypeII < realHaplotypeIndices.size(); haplotypeII++)
		{
			int haplotypeI = realHaplotypeIndices.at(haplotypeII);
			if(!(seen_haplotypes.count(haplotypeI) > 0))
			{
				cout << "level: " << level << ", haplotypeI: " << haplotypeI << ", seen_haplotypes.size(): " << seen_haplotypes.size() << ", hp.HaplotypeIDs.size(): " << hp.HaplotypeIDs.size() << "\n";
				cout << flush;
			}
			assert(seen_haplotypes.count(haplotypeI) > 0);
		}

		NodeToHaplotype.clear();

		int addedNodeCounter = 0;
		for(int i = 0; i < (int)groupingEdges.size(); i++)
		{
			if(groupsDeleted.count(i) > 0)
			{
				continue;
			}
			else
			{
				addedNodeCounter++;

				Node* newN = new Node();
				if(level == (int)(loci.size()-1))
				{
					newN->terminal = true;
				}
				else
				{
					newN->terminal = false;
				}

				newN->level = level + 1;
				registerNode(newN, newN->level);

				for(set<Edge*>::iterator eIt = groupingEdges.at(i).begin(); eIt != groupingEdges.at(i).end(); eIt++)
				{
					newN->Incoming_Edges.insert(*eIt);
					(*eIt)->To = newN;

					if(snpEdges.count(*eIt) > 0)
					{
						for(set<Edge*>::iterator snpEdgeIt = snpEdges[*eIt].begin(); snpEdgeIt != snpEdges[*eIt].end(); snpEdgeIt++)
						{
							Edge* snpEdge = *snpEdgeIt;
							newN->Incoming_Edges.insert(snpEdge);
							snpEdge->To = newN;
							snpEdges_to_catch.erase(snpEdge);
						}
					}
				}


				assert(groupingHaplotypes.size() > 0);
				NodeToHaplotype[newN] = groupingHaplotypes.at(i);
				// cout << "New node memory address: " << newN << "\n";
				for(set<int>::iterator hIndex = groupingHaplotypes.at(i).begin(); hIndex != groupingHaplotypes.at(i).end(); hIndex++)
				{
					// cout << "\t associated haplotype index : " << *hIndex << "\n";
				}
				// cout << flush;
			}
		}

		assert(snpEdges_to_catch.size() == 0);

	}

	// todo re-activate
	// g->removeStarPaths();

	if(wantPGFprotection)
	{
		for(unsigned int level = 0; level < (NodesPerLevel.size()-1); level++)
		{
			bool foundPGF = false;
			for(set<Node*>::iterator nodeIt = NodesPerLevel.at(level).begin(); nodeIt != NodesPerLevel.at(level).end(); nodeIt++)
			{
				Node* node = *nodeIt;

				for(set<Edge*>::iterator eIt = node->Outgoing_Edges.begin(); eIt != node->Outgoing_Edges.end(); eIt++)
				{
					Edge* e = *eIt;
					foundPGF = (foundPGF || (bool)e->pgf_protect);
					if((level > 0) && (level < (NodesPerLevel.size()-2)))
					{
						if((bool)e->pgf_protect)
						{
							bool found_preceding_PGF = false;

							for(set<Edge*>::iterator beforeIt = e->From->Incoming_Edges.begin(); beforeIt != e->From->Incoming_Edges.end(); beforeIt++)
							{
								if((bool)(*beforeIt)->pgf_protect)
								{
									found_preceding_PGF = true;
									break;
								}
							}
							assert(found_preceding_PGF);

							bool found_following_PGF = false;
							for(set<Edge*>::iterator afterIt = e->To->Outgoing_Edges.begin(); afterIt != e->To->Outgoing_Edges.end(); afterIt++)
							{
								if((bool)(*afterIt)->pgf_protect)
								{
									found_following_PGF = true;
									break;
								}
							}

							if(! found_following_PGF)
							{
								cout << "Level " << level << " of " << (NodesPerLevel.size()-2) << ", edge " << e << " pgf protect: " << (bool)e->pgf_protect << " starting at level " << e->From->level << "\n";
								for(set<Edge*>::iterator afterIt = e->To->Outgoing_Edges.begin(); afterIt != e->To->Outgoing_Edges.end(); afterIt++)
								{
									cout << "\t following " << (*afterIt) << ": " << (bool)(*afterIt)->pgf_protect << "\n" << flush;
									assert( e->To->level == (level + 1));
								}

								cout << "We now look at all edges from the this level " << (level) << ":\n";
								for(set<Node*>::iterator nodeIt2 = NodesPerLevel.at(level).begin(); nodeIt2 != NodesPerLevel.at(level).end(); nodeIt2++)
								{
									Node* node2 = *nodeIt2;

									for(set<Edge*>::iterator eIt2 = node2->Outgoing_Edges.begin(); eIt2 != node2->Outgoing_Edges.end(); eIt2++)
									{
										Edge* e2 = *eIt2;
										cout << "\t general " << e2 << ": " << (bool)e2->pgf_protect << "\n" << flush;
									}
								}

								cout << "We now look at all edges from the next level " << (level+1) << ":\n";
								for(set<Node*>::iterator nodeIt2 = NodesPerLevel.at(level+1).begin(); nodeIt2 != NodesPerLevel.at(level+1).end(); nodeIt2++)
								{
									Node* node2 = *nodeIt2;

									for(set<Edge*>::iterator eIt2 = node2->Outgoing_Edges.begin(); eIt2 != node2->Outgoing_Edges.end(); eIt2++)
									{
										Edge* e2 = *eIt2;
										cout << "\t general " << e2 << ": " << (bool)e2->pgf_protect << "\n" << flush;
									}
								}
							}

							assert(found_following_PGF);
						}
					}
				}
			}
			assert(foundPGF);
		}
	}

//	for(set<Node*>::iterator nodeIt = NodesPerLevel.at(NodesPerLevel.size()-1).begin(); nodeIt != NodesPerLevel.at(NodesPerLevel.size()-1).end(); nodeIt++)
//	{
//		set<Edge*> nodeIncomingEdges = (*nodeIt)->Incoming_Edges;
//		for(set<Edge*>::iterator edgeIt = nodeIncomingEdges.begin(); edgeIt != nodeIncomingEdges.end(); edgeIt++)
//		{
//			Edge* edge = *edgeIt;
//			string emission = edge->getEmission();
//			assert(emission == "*");
//		}
//	}

	for(set<Node*>::iterator nodeIt = NodesPerLevel.at(NodesPerLevel.size()-1).begin(); nodeIt != NodesPerLevel.at(NodesPerLevel.size()-1).end(); nodeIt++)
	{
		set<Edge*> nodeIncomingEdges = (*nodeIt)->Incoming_Edges;
		for(set<Edge*>::iterator edgeIt = nodeIncomingEdges.begin(); edgeIt != nodeIncomingEdges.end(); edgeIt++)
		{
			Edge* edge = *edgeIt;
			string emission = edge->getEmission();
			assert(emission != "_");
		}
	}

	// cout << "\n\nBuilding nucleotide graph done!\n\n" << flush;
}


std::string Graph::getOneLocusIDforLevel(unsigned int level)
{
	std::string locus;
	
	if(level != (NodesPerLevel.size()-1))
	{	
		assert(NodesPerLevel.at(level).size() > 0);
		Node* n = *(NodesPerLevel.at(level).begin());

		assert(n->Outgoing_Edges.size() > 0);
		locus = (*(n->Outgoing_Edges.begin()))->locus_id;
	}

	return locus;
}


void Graph::graphViz(int level_start, int level_stop, std::string output_filename)
{
	std::ofstream output;
	output.open(output_filename.c_str());
	assert(output.is_open());

	output << "digraph G {\n";

	assert(level_stop > level_start);

	std::vector<std::map<Node*, std::string> > _node_2_int;
	_node_2_int.resize((level_stop - level_start + 1));

	auto getIDForNode = [&](int level, Node* n) -> std::string {
		assert(level >= level_start);
		assert(level <= level_stop);
		if(_node_2_int.at(level-level_start).count(n) == 0)
		{
			int existingNodes = _node_2_int.at(level-level_start).size();
			int thisNode = existingNodes + 1;
			std::string nodeID = "L" + Utilities::ItoStr(level) + "N" + Utilities::ItoStr(thisNode);
			_node_2_int.at(level-level_start)[n] = nodeID;
		}
		return _node_2_int.at(level-level_start).at(n);
	};

	for(int level = level_start; level < level_stop; level++)
	{
		std::set<Node*> nodes_thisLevel_set = NodesPerLevel.at(level);
		std::vector<Node*> nodes_thisLevel(nodes_thisLevel_set.begin(), nodes_thisLevel_set.end());

		for(unsigned int nI = 0; nI < nodes_thisLevel.size(); nI++)
		{
			Node* n = nodes_thisLevel.at(nI);
			std::string n_ID = getIDForNode(level, n);
			for(std::set<Edge*>::iterator eIt = n->Outgoing_Edges.begin(); eIt != n->Outgoing_Edges.end(); eIt++)
			{
				Edge* e = *eIt;
				Node* n2 = e->To;
				std::string n2_ID = getIDForNode(level+1, n2);

				std::string e_emission = e->getEmission();
				output << "\t" << n_ID << " -> " << n2_ID << "[label=\"" << e_emission <<  "\"]" << "\n";
			}
		}
	}

	output << "\n" << "}" << "\n";
	output.close();

}


void Graph::graphViz2(std::string locus_string, std::string output_filename)
{
	std::ofstream output;
	output.open(output_filename.c_str());
	assert(output.is_open());

	std::set<unsigned int> levels;
	for(unsigned int l = 0; l < NodesPerLevel.size(); l++)
	{
		std::set<Node*> nodes_at_l = NodesPerLevel.at(l);
		for(std::set<Node*>::iterator nIt = nodes_at_l.begin(); nIt != nodes_at_l.end(); nIt++)
		{
			Node* n = *nIt;
			assert(n->level == l);
			for(std::set<Edge*>::iterator edgeIt = n->Outgoing_Edges.begin(); edgeIt != n->Outgoing_Edges.end(); edgeIt++ )
			{
				Edge* e = *edgeIt;
				if(e->locus_id.find(locus_string) != std::string::npos)
				{
					if(e->locus_id.find("buffer") == std::string::npos)
					{
						if(
							(e->locus_id.find("exon/2_") != std::string::npos)
//							|| (e->locus_id.find("exon/3_") != std::string::npos)
						)
						{
							levels.insert(l);
						}
					}
				}
			}
		}
	}

	assert(levels.size() > 1);

	int level_start = (int)(*levels.begin());
	int level_stop = (int)(*levels.rbegin());
	int levels_extract = level_stop - level_start + 1;
	
	std::cout << "Graph::graphViz2(..): search string " << locus_string << ", extract from " << level_start << " to " << level_stop << " (" << levels_extract << ") levels\n" << std::flush;
	
	output << "digraph G {\n";
	output << "\tgraph [rotate=90]\n";
	output << "\tnode [label=\"\", height = 0.1, width = 0.1]\n";

	assert(level_stop > level_start);
	
	std::vector<std::map<Node*, std::string> > _node_2_int;
	_node_2_int.resize((level_stop - level_start + 1));

	auto getIDForNode = [&](int level, Node* n) -> std::string {
		assert(level >= level_start);
		assert(level <= level_stop);
		if(_node_2_int.at(level-level_start).count(n) == 0)
		{
			int existingNodes = _node_2_int.at(level-level_start).size();
			int thisNode = existingNodes + 1;
			std::string nodeID = "L" + Utilities::ItoStr(level) + "N" + Utilities::ItoStr(thisNode);
			_node_2_int.at(level-level_start)[n] = nodeID;
		}
		return _node_2_int.at(level-level_start).at(n);
	};

	for(int level = level_start; level < level_stop; level++)
	{
		std::set<Node*> nodes_thisLevel_set = NodesPerLevel.at(level);
		std::vector<Node*> nodes_thisLevel(nodes_thisLevel_set.begin(), nodes_thisLevel_set.end());

		for(unsigned int nI = 0; nI < nodes_thisLevel.size(); nI++)
		{
			Node* n = nodes_thisLevel.at(nI);
			std::string n_ID = getIDForNode(level, n);
			for(std::set<Edge*>::iterator eIt = n->Outgoing_Edges.begin(); eIt != n->Outgoing_Edges.end(); eIt++)
			{
				Edge* e = *eIt;
				Node* n2 = e->To;
				std::string n2_ID = getIDForNode(level+1, n2);

				std::string e_emission = e->getEmission();
				output << "\t" << n_ID << " -> " << n2_ID << "[label=\"" << e_emission <<  "\"]" << "\n";
			}
		}
	}

	output << "\n" << "}" << "\n";
	output.close();

}

void Graph::unRegisterNode(Node* n)
{
	assert(Nodes.count(n) > 0);
	int l = n->level;
	Nodes.erase(n);
	NodesPerLevel.at(l).erase(n);
	delete(n);
}

void Graph::makeEdgesGaps(double proportion)
{
	srand ( time(NULL) );

	assert(proportion >= 0);
	assert(proportion <= 1);

	for(set<Edge*>::iterator edgeIt = Edges.begin(); edgeIt != Edges.end(); edgeIt++)
	{
		Edge* e = (*edgeIt);
		double f = (double)rand() / RAND_MAX;
		if(f <= proportion)
		{
			unsigned char newEmission = '_';
			e->emission = newEmission;
		}
	}
}

std::vector<std::string> Graph::simulateHaplotypes(int number, bool includeGaps)
{
	// srand ( time(NULL) ); // todo activate
	assert(NodesPerLevel.at(0).size()==1);

	std::vector<std::string> forReturn;

	for(int iteration = 0; iteration < number; iteration++)
	{
		vector<string> symbols;
		int levels = NodesPerLevel.size();
		Node* currentNode = *(NodesPerLevel.at(0).begin());
		while((int)currentNode->level != (levels-1))
		{
			double To_Outgoing_Sum = currentNode->Sum_Outgoing();
			assert(To_Outgoing_Sum > 0);
			vector<double> stateP;
			vector<Edge*> outEdges( currentNode->Outgoing_Edges.begin(), currentNode->Outgoing_Edges.end());
			assert(outEdges.size() > 0);
			for(int eI = 0; eI < (int)outEdges.size(); eI++)
			{
				stateP.push_back(outEdges.at(eI)->count/To_Outgoing_Sum);
			}
			int chosenEdgeI = Utilities::chooseFromVector(stateP);
			Edge* chosenEdge = outEdges.at(chosenEdgeI);
			string symbolToEmit = chosenEdge->getEmission();

			if((symbolToEmit != "_") || includeGaps)
				symbols.push_back(symbolToEmit);

			currentNode = chosenEdge->To;
		}

		forReturn.push_back(Utilities::join(symbols, ""));
		// cout << iteration << " " << Utilities::join(symbols, "") << "\n";
	}

	return forReturn;
}


diploidEdgePointerPath Graph::simulateRandomDiploidPath()
{
	vector< vector<Edge*> > edgePaths;
	for(int pI = 0; pI < 2; pI++)
	{
		vector<Edge*> currentEdgePath;

		Node* currentNode;
		Node* n0 = *(NodesPerLevel.at(0).begin());

		currentNode = n0;
		while(currentNode->Outgoing_Edges.size() != 0)
		{
			int n_edges = currentNode->Outgoing_Edges.size();
			vector<Edge*> currentEdges (currentNode->Outgoing_Edges.begin(), currentNode->Outgoing_Edges.end());

			double f = (double)rand() / RAND_MAX;

			f = f * n_edges;
			int selected_edge = (int) f;
			if(selected_edge == n_edges)
			{
				selected_edge = n_edges - 1;
			}

			assert(selected_edge >= 0);
			assert(selected_edge < n_edges);

			Edge* selectedEdge = currentEdges.at(selected_edge);

			currentEdgePath.push_back(selectedEdge);
			currentNode = selectedEdge->To;
		}

		assert(currentEdgePath.size() == (NodesPerLevel.size() - 1));

		edgePaths.push_back(currentEdgePath);
	}

	diploidEdgePointerPath forReturn;
	forReturn.h1 = edgePaths.at(0);
	forReturn.h2 = edgePaths.at(1);

	return forReturn;
}

int Graph::trimGraph(bool remove2DHLA)
{
	int removeNodes = 0;
	int removeEdges = 0;

	set<Node*> leadToEnd;
	set<Node*> leadToStart;

	set<Edge*> zeroEdges;

	for(int l = NodesPerLevel.size()-1; l >= 0; l--)
	{
		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(l).begin(); nodeIt !=  NodesPerLevel.at(l).end(); nodeIt++)
		{
			Node* n = *nodeIt;
			if(l == (int)NodesPerLevel.size()-1)
			{
				leadToEnd.insert(n);
			}
			else
			{
				double To_Outgoing_Sum = n->Sum_Outgoing();
				if(To_Outgoing_Sum > 0)
				{
					for(set<Edge*>::iterator outgoingIt = n->Outgoing_Edges.begin(); outgoingIt != n->Outgoing_Edges.end(); outgoingIt++)
					{
						Edge* e = (*outgoingIt);
						double p = e->count / To_Outgoing_Sum;
						if(p > 0)
						{
							if(leadToEnd.count(e->To) > 0)
							{
								leadToEnd.insert(n);
								assert(1 == 0);
								// when this is called, check that what happens is intentional!
							}
						}
						else
						{
							zeroEdges.insert(e);
						}
					}
				}
			}
		}
	}


	for(int l = 0; l < (int)NodesPerLevel.size(); l++)
	{
		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(l).begin(); nodeIt !=  NodesPerLevel.at(l).end(); nodeIt++)
		{
			Node* n = *nodeIt;
			if(l == 0)
			{
				leadToStart.insert(n);
			}

			double To_Outgoing_Sum = n->Sum_Outgoing();
			if(To_Outgoing_Sum > 0)
			{
				for(set<Edge*>::iterator outgoingIt = n->Outgoing_Edges.begin(); outgoingIt != n->Outgoing_Edges.end(); outgoingIt++)
				{
					Edge* e = (*outgoingIt);
					double p = e->count / To_Outgoing_Sum;
					if(p > 0)
					{
						if(leadToStart.count(n) > 0)
						{
							leadToStart.insert(e->To);
							assert(1 == 0);
							// when this happens, again, check that this is intentional!
						}
					}
					else
					{
						zeroEdges.insert(e);
					}
				}
			}
		}
	}


	for(set<Edge*>::iterator edgeRemoveIt = zeroEdges.begin(); edgeRemoveIt != zeroEdges.end(); edgeRemoveIt++)
	{
		Edge* e = *edgeRemoveIt;
		e->From->Outgoing_Edges.erase(e);
		e->To->Incoming_Edges.erase(e);
		unRegisterEdge(e);
		removeEdges++;
	}

	for(set<Node*>::iterator nodeIt = Nodes.begin(); nodeIt != Nodes.end(); nodeIt++)
	{
		Node* n = *nodeIt;
		if((leadToEnd.count(n) == 0) || (leadToStart.count(n) == 0))
		{
			for(set<Edge*>::iterator edgeRemoveIt = n->Outgoing_Edges.begin(); edgeRemoveIt != n->Outgoing_Edges.end(); edgeRemoveIt++)
			{
				Edge* e = *edgeRemoveIt;
				if(Edges.count(e) > 0)
				{
					e->From->Outgoing_Edges.erase(e);
					e->To->Incoming_Edges.erase(e);
					unRegisterEdge(e);
					removeEdges++;
				}
			}
			for(set<Edge*>::iterator edgeRemoveIt = n->Incoming_Edges.begin(); edgeRemoveIt != n->Incoming_Edges.end(); edgeRemoveIt++)
			{
				Edge* e = *edgeRemoveIt;
				if(Edges.count(e) > 0)
				{
					e->From->Outgoing_Edges.erase(e);
					e->To->Incoming_Edges.erase(e);
					unRegisterEdge(e);
					removeEdges++;
				}
			}
			unRegisterNode(n);
			removeNodes++;
		}
	}

	checkConsistency(false);

	cerr << "Trimming: removed " << removeEdges << " edges and " << removeNodes << " nodes.\n";
	return removeEdges+removeNodes;
}

void Graph::unRegisterEdge(Edge* e)
{
	assert(Edges.count(e) > 0);
	Edges.erase(e);
	delete(e);
}

void Graph::registerNode(Node* n, unsigned int level)
{
	assert(Nodes.count(n) == 0);
	assert(n->level == level);
	if((level+1) > NodesPerLevel.size())
	{
		NodesPerLevel.resize(level+1);
	}
	NodesPerLevel.at(level).insert(n);
	Nodes.insert(n);
	n->g = this;
}

vector<string> Graph::getAssignedLoci()
{
	vector<string> loci(NodesPerLevel.size()-1, "");
	for(set<Edge*>::iterator eIt = Edges.begin(); eIt != Edges.end(); eIt++)
	{
		Edge* e = *eIt;
		Node* n1 = e->From;
		//Node* n2 = e->To;
		
		int level = n1->level;
		string locusID = e->locus_id;
		
		if(loci.at(level) != "")
		{
			assert(loci.at(level) == locusID);
		}
		else
		{
			loci.at(level) = locusID;
		}
	}
	return loci;
}

void Graph::checkLocusOrderConsistency(vector<string> loci)
{
	assert((NodesPerLevel.size()-1) == loci.size());
	bool all_ok = true;
	for(set<Edge*>::iterator eIt = Edges.begin(); eIt != Edges.end(); eIt++)
	{
		Edge* e = *eIt;
		Node* n1 = e->From;
		//Node* n2 = e->To;
		
		int level = n1->level;
		
		if(loci.at(level) != e->locus_id)
		{
			all_ok = false;
			cerr << "Level " << level << " edge says locus " << e->locus_id << ", but locus list says " << loci.at(level) << "\n"; 
		}
	}
	
	if(! all_ok)
	{
		cerr << "Locus list:\n";
		for(int j = 0; j < (int)loci.size(); j++)
		{
			cerr << j << " " << loci.at(j) << "\n";
		}
	}

	assert(all_ok == true);
}


void Graph::printComplexity (string filename)
{
    using boost::lexical_cast;
    using boost::bad_lexical_cast;
    
	ofstream output;
	output.open (filename.c_str(), ios::out | ios::trunc);
	if (output.is_open())
	{
		output << "Level\tNodes\tEdges\n";

		for(unsigned int level = 0; level < NodesPerLevel.size(); level++)
		{
			int nodes = NodesPerLevel.at(level).size();
			int edges = 0;
			for(set<Node*>::iterator nodeIt = NodesPerLevel.at(level).begin(); nodeIt != NodesPerLevel.at(level).end(); nodeIt++)
			{
				Node* node = *nodeIt;
				for(set<Edge*>::iterator eIt = node->Outgoing_Edges.begin(); eIt != node->Outgoing_Edges.end(); eIt++)
				{
					edges++;
				}
			}

			output << level << "\t total nodes: " << nodes << "\t total edges: " << edges << "\n";

			for(set<Node*>::iterator nodeIt = NodesPerLevel.at(level).begin(); nodeIt != NodesPerLevel.at(level).end(); nodeIt++)
			{
				Node* node = *nodeIt;
				output << "\t\t Outgoing edges: " << node->Outgoing_Edges.size() << "\t Incoming edges: " << node->Incoming_Edges.size() << "\n";
				for(set<Edge*>::iterator eIt = node->Outgoing_Edges.begin(); eIt != node->Outgoing_Edges.end(); eIt++)
				{
					Edge* e = *eIt;
					output << "\t\t\t\t" << e->getEmission() << "          from " << e->From << "to " << e->To << "    [edge " << e << "]\n";
				}
			}
		}
	}
	else
	{
		throw std::runtime_error("Cannot open output file for graph serialization: "+filename);
	}



}



void Graph::checkConsistency(bool terminalCheck)
{
    using boost::lexical_cast;
    using boost::bad_lexical_cast;

	// all nodes implied by the Edges set are in the Nodes set
	// all incoming/outgoing relations implied by the Edges are reflected in the nodes
	// all edges have incoming and outgoing nodes
	for(set<Edge*>::iterator eIt = Edges.begin(); eIt != Edges.end(); eIt++)
	{
		Edge* e = *eIt;
		Node* n1 = e->From;
		Node* n2 = e->To;

		string emission = lexical_cast<string>(e->emission);
		if(emission.length() == 0)
		{
			cerr << "Problem transforming emission to string -- emission:" << e->emission << "(end)\n";
			throw std::runtime_error("Unsigned Char Conversion problem!");
		}

		assert(e->count >= 0);


		assert(n1 != NULL);
		if(n1 != NULL)
		{
			assert(count(n1->Outgoing_Edges.begin(), n1->Outgoing_Edges.end(), e) > 0);
		}

		assert(n2 != NULL);
		if(n2 != NULL)
		{
			assert(count(n2->Incoming_Edges.begin(), n2->Incoming_Edges.end(), e) > 0);
		}

		assert(n1->level == (n2->level-1));

		assert(Nodes.count(n1) > 0);
		assert(Nodes.count(n2) > 0);
	}

	// The NodesPerLevel set and the Nodes set contain the same nodes
	// The edges implied by the Nodes are in the Edges set
	set<Node*> saw_node;
	for(unsigned int level = 0; level < NodesPerLevel.size(); level++)
	{
		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(level).begin(); nodeIt != NodesPerLevel.at(level).end(); nodeIt++)
		{
			Node* node = *nodeIt;
			saw_node.insert(node);

			assert(node->level == level);
			assert(Nodes.count(node) > 0);

			for(set<Edge*>::iterator eIt = node->Incoming_Edges.begin(); eIt != node->Incoming_Edges.end(); eIt++)
			{
				Edge* e = *eIt;
				assert(Edges.count(e) > 0);
			}

			for(set<Edge*>::iterator eIt = node->Outgoing_Edges.begin(); eIt != node->Outgoing_Edges.end(); eIt++)
			{
				Edge* e = *eIt;
				assert(Edges.count(e) > 0);
			}
		}
	}


	// there are no loneley nodes
	for(set<Node*>::iterator nodeIt = Nodes.begin(); nodeIt != Nodes.end(); nodeIt++)
	{
		Node* node = *nodeIt;
		assert(saw_node.count(node) > 0);
		assert((! terminalCheck) || (node->terminal) || (node->Outgoing_Edges.size() > 0));
		assert((node->level == 0) || (node->Incoming_Edges.size() > 0));
	}
}

void Graph::freeMemory()
{
	vector<Edge*> edges_for_free(Edges.begin(), Edges.end());
	vector<Node*> nodes_for_free(Nodes.begin(), Nodes.end());

	for(unsigned int i = 0; i < edges_for_free.size(); i++)
	{
		Edge* e = edges_for_free.at(i);
		delete(e);
	}

	for(unsigned int i = 0; i < nodes_for_free.size(); i++)
	{
		Node* n = nodes_for_free.at(i);
		delete(n);
	}
}


vector<levelInfo> Graph::getLevelInfo()
{
	vector<levelInfo> forReturn;

	for(int l = 0; l < (int)NodesPerLevel.size(); l++)
	{
		int nodes = NodesPerLevel.at(l).size();
		int edges = 0;
		int symbols = 0;

		std::set<std::string> edge_emissions;

		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(l).begin(); nodeIt !=  NodesPerLevel.at(l).end(); nodeIt++)
		{
			Node* n = *nodeIt;
			edges = edges + n->Outgoing_Edges.size();

			if(l != (int)(NodesPerLevel.size()-1))
			{
				string locus = (*(n->Outgoing_Edges.begin()))->locus_id;


				for(set<Edge*>::iterator outgoingIt = n->Outgoing_Edges.begin(); outgoingIt != n->Outgoing_Edges.end(); outgoingIt++)
				{
					assert((*outgoingIt)->locus_id == locus);
					edge_emissions.insert((*outgoingIt)->getEmission());
				}
			}
		}

		symbols = edge_emissions.size();

		levelInfo lI;
		lI.nodes = nodes;
		lI.edges = edges;
		lI.symbols = symbols;

		forReturn.push_back(lI);
	}

	return forReturn;
}


void Graph::registerEdge(Edge* e)
{
	assert(Edges.count(e) == 0);
	Edges.insert(e);
}


void Graph::removeStarPaths()
{

	// Find edges that are the start of paths that we potentially may want to remove
	// Such paths are defined as: there is an edge with a star symbol, and there are
	// other edges with non-star symbols emanating from the same node

	map<Edge*, bool> edgesPotentialPathStarts;
	for(int level = 0; level < (int)NodesPerLevel.size(); level++)
	{
		string locusID = "";
		Node* firstNode = *(NodesPerLevel.at(level).begin());

		if(firstNode->Outgoing_Edges.size() == 0)
		{
			assert(level == (int)(NodesPerLevel.size() - 1));
		}

		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(level).begin(); nodeIt != NodesPerLevel.at(level).end(); nodeIt++)
		{
			Node* n = *nodeIt;
			set<Node*> has_non_star_edge_incoming;
			set<Edge*> star_edges_outgoing;
			for(set<Edge*>::iterator edgeIt = n->Outgoing_Edges.begin(); edgeIt != n->Outgoing_Edges.end(); edgeIt++)
			{
				assert(locusID != "");
				Edge* e = *edgeIt;
				if(e->getEmission() == "*")
				{
					star_edges_outgoing.insert(e);
				}
				else
				{
					has_non_star_edge_incoming.insert(e->To);
				}
			}

			if(n->Outgoing_Edges.size() > star_edges_outgoing.size())
			{
				for(set<Edge*>::iterator edgeIt = star_edges_outgoing.begin(); edgeIt != star_edges_outgoing.end(); edgeIt++)
				{
					Edge* e = *edgeIt;
					if(! e->pgf_protect)
					{
						edgesPotentialPathStarts[e] = true;
					}
				}
			}
		}
	}

	// Find out whether these potential start nodes actually
	// define a path

	set<Node*> delete_Nodes;
	set<Edge*> delete_Edges;
	set<Edge*> edgeFromOtherPath;

	for(map<Edge*, bool>::iterator edgeIt = edgesPotentialPathStarts.begin(); edgeIt != edgesPotentialPathStarts.end(); edgeIt++)
	{
		Edge* e = edgeIt->first;

		if(edgeIt->second != true)
		{
			continue;
		}
		if(delete_Edges.count(e) > 0)
		{
			continue;
		}

		// openPaths contains the paths that we may still want to elongate
		// closedPaths contains the paths that we are still processing
		// the path criterion is: on the way only starred edges + ending in a
		// node that has other incoming edges than the one from the path

		vector< vector<Edge*> > openPaths;
		vector< vector<Edge*> > closedPaths;

		vector<Edge*> firstPath;

		firstPath.push_back(e);
		openPaths.push_back(firstPath);
		set<Edge*> edgeFromPath;
		edgeFromPath.insert(e);

		while(openPaths.size() > 0)
		{
			int extendedPathLength = -1;

			int openPathsSize = openPaths.size();
			for(int openPathI = 0; openPathI < openPathsSize; openPathI++)
			{
				int lastEdgeIndex = openPaths.at(openPathI).size()-1;
				Edge* lastEdge = openPaths.at(openPathI).at(lastEdgeIndex);
				Node* lastEdgeTo = lastEdge->To;

				// legitimate stop nodes are only those which have another
				// incoming edge which does not come this path

				// if we find a legitimate stop node,
				// we end this path and break
				bool stopNode = false;
				for(set<Edge*>::iterator edgeIt = lastEdgeTo->Incoming_Edges.begin(); edgeIt != lastEdgeTo->Incoming_Edges.end(); edgeIt++)
				{
					if((edgeFromPath.count(*edgeIt) == 0) && (edgeFromOtherPath.count(*edgeIt) == 0))
					{
						stopNode = true;
						closedPaths.push_back(openPaths.at(openPathI));
						openPaths.at(openPathI).clear();
						break;
					}
				}

				// if the last node of this path is not legitimate stop node,
				// we extend
				if(! stopNode)
				{
					// find edges to iterate over... if any of the edges is no star symbol, we
					// need to abort ...
					bool pursueEdges = true;
					for(set<Edge*>::iterator edgeIt = lastEdgeTo->Outgoing_Edges.begin(); edgeIt != lastEdgeTo->Outgoing_Edges.end(); edgeIt++)
					{
						Edge* e = *edgeIt;
						if(e->getEmission() != "*")
						{
							pursueEdges = false;
						}
						if(e->pgf_protect)
						{
							pursueEdges = false;
						}
					}
					// ... same thing if we have no further edges to follow.
					if(lastEdgeTo->Outgoing_Edges.size() == 0)
					{
						pursueEdges = false;
					}

					// if we abort, we clear all found paths ... none is valid; and we break.
					if(! pursueEdges)
					{
						openPaths.clear();
						closedPaths.clear();
						break;
					}

					// on the other hand, if we do not have to abort, we follow the edges
					vector<Edge*> pursueEdgesVector(lastEdgeTo->Outgoing_Edges.begin(), lastEdgeTo->Outgoing_Edges.end());
					vector<Edge*> originalPath = openPaths.at(openPathI);
					for(int pursueEdgesI = 0; pursueEdgesI < (int)pursueEdgesVector.size(); pursueEdgesI++)
					{
						if(pursueEdgesI == 0)
						{
							openPaths.at(openPathI).push_back(pursueEdgesVector.at(pursueEdgesI));
							if(extendedPathLength == -1)
							{
								extendedPathLength = openPaths.at(openPathI).size();
							}
							else
							{
								assert(extendedPathLength == (int)openPaths.at(openPathI).size());
							}
						}
						else
						{
							openPaths.push_back(originalPath);
							openPaths.at(openPaths.size()-1).push_back(pursueEdgesVector.at(pursueEdgesI));
							assert(extendedPathLength == (int)openPaths.at(openPaths.size()-1).size());
						}

						edgeFromPath.insert(pursueEdgesVector.at(pursueEdgesI));
					}
				}
			}

			// if there are any open paths left, we will have to follow them
			for(int openPathI = openPaths.size(); openPathI > 0; openPathI--)
			{
				if(openPaths.at(openPathI-1).size() == 0)
				{
					openPaths.erase(openPaths.begin()+openPathI-1);
				}
				else
				{
					assert(extendedPathLength != -1);
					assert((int)openPaths.at(openPathI-1).size() == extendedPathLength);
				}
			}
		}

		assert(openPaths.size() == 0);

		for(size_t pI = 0; pI < closedPaths.size(); pI++)
		{
			for(size_t eI = 0; eI < closedPaths.at(pI).size(); eI++)
			{
				Edge* e = closedPaths.at(pI).at(eI);
				delete_Edges.insert(e);
				edgeFromOtherPath.insert(e);
				if(eI == 0)
				{
					assert(e->From->Outgoing_Edges.size() > 1);
				}
				if(eI != 0)
				{
					delete_Nodes.insert(e->From);
				}
				if(eI != (closedPaths.at(pI).size() - 1))
				{
					delete_Nodes.insert(e->To);
				}
				if(eI == (closedPaths.at(pI).size() - 1))
				{
					assert(e->To->Incoming_Edges.size() > 1);
				}
			}
		}
	}

	int removedEdges = 0;
	int removedNodes = 0;
	//map<Node*, vector<string> > deletedKmers_cache;
	for(set<Edge*>::iterator edgeRemoveIt = delete_Edges.begin(); edgeRemoveIt != delete_Edges.end(); edgeRemoveIt++)
	{
		Edge* e = *edgeRemoveIt;
		assert(e->getEmission() == "*");
		//deletedKmers_cache[e->To].push_back(CODE.deCode(e->locus_id, e->emission));
		e->From->Outgoing_Edges.erase(e);
		e->To->Incoming_Edges.erase(e);
		unRegisterEdge(e);
		removedEdges++;
	}
	for(set<Node*>::iterator nodeRemoveIt = delete_Nodes.begin(); nodeRemoveIt != delete_Nodes.end(); nodeRemoveIt++)
	{
		Node* n = *nodeRemoveIt;
		assert(n->Outgoing_Edges.size() == 0);
		assert(n->Incoming_Edges.size() == 0);
		unRegisterNode(n);
		removedNodes++;
	}

	for(set<Node*>::iterator nodeIt = Nodes.begin(); nodeIt != Nodes.end(); nodeIt++)
	{
		Node* node = *nodeIt;
		if(!(((node->level == 0) || (node->Incoming_Edges.size() > 0))))
		{
			cout << node << "\n";
//			assert(deletedKmers_cache.count(node) > 0);
//			for(int i = 0; i < (int)deletedKmers_cache[node].size(); i++)
//			{
//				cout << "\t" << deletedKmers_cache[node].at(i) << "\n";
//			}
		}

		assert((node->level == 0) || (node->Incoming_Edges.size() > 0));

	}

	cout << "Star removal: removed " << removedNodes << " nodes and " << removedEdges << " edges.\n";

	checkConsistency(false);
}

LocusCodeAllocation Graph::constructCODEfromGraph()
{
	LocusCodeAllocation CODE_forReturn;

	for(size_t l = 0; l < NodesPerLevel.size(); l++)
	{
		for(set<Node*>::iterator nodeIt = NodesPerLevel.at(l).begin(); nodeIt !=  NodesPerLevel.at(l).end(); nodeIt++)
		{
			Node* n = *nodeIt;

			if(l != (NodesPerLevel.size()-1))
			{
				string locus = (*(n->Outgoing_Edges.begin()))->locus_id;

				for(set<Edge*>::iterator outgoingIt = n->Outgoing_Edges.begin(); outgoingIt != n->Outgoing_Edges.end(); outgoingIt++)
				{
					assert((*outgoingIt)->locus_id == locus);
					CODE_forReturn.doCode(locus, (*outgoingIt)->getEmission());
				}
			}
		}
	}

	return CODE_forReturn;
}


void Graph::writeToFile(string filename)
{
    using boost::lexical_cast;
    using boost::bad_lexical_cast;

	checkConsistency(false);

	LocusCodeAllocation localCODE = constructCODEfromGraph();

	vector<string> linesFromCode = localCODE.serializeIntoVector();
	vector<string> linesForNodes;
	vector<string> linesForEdges;

	map<Edge*, int> EdgeCounter;
	map<Node*, int> NodeCounter;

	for(set<Node*>::iterator nodeIt = Nodes.begin(); nodeIt != Nodes.end(); nodeIt++)
	{
		Node* n = *nodeIt;

		if(NodeCounter.count(n) == 0)
		{
			int s = NodeCounter.size();
			NodeCounter[n] = s + 1;
		}
		int n_index = NodeCounter[n];

		vector<string> line_fields;
		line_fields.push_back(lexical_cast<string>(n_index));
		line_fields.push_back(lexical_cast<string>(n->level));
		line_fields.push_back(lexical_cast<string>(n->terminal));

		linesForNodes.push_back(boost::join(line_fields, separatorForSerialization));
	}

	for(set<Edge*>::iterator edgeIt = Edges.begin(); edgeIt != Edges.end(); edgeIt++)
	{
		Edge* e = *edgeIt;

		if(EdgeCounter.count(e) == 0)
		{
			int s = EdgeCounter.size();
			EdgeCounter[e] = s + 1;
		}
		int e_index = EdgeCounter[e];

		vector<string> line_fields;
		line_fields.push_back(lexical_cast<string>(e_index));
		line_fields.push_back(lexical_cast<string>(e->locus_id));
		line_fields.push_back(lexical_cast<string>(e->count));
		string emission = boost::lexical_cast<std::string>(localCODE.doCode(e->locus_id, e->getEmission()));

		if(emission.length() == 0)
		{
			cerr << "Problem transforming emission to string -- emission:" << e->emission << "(end)\n"; 
			throw std::runtime_error("Conversion problem!");
		}
		line_fields.push_back(emission);


		Node* n1 = e->From;
		Node* n2 = e->To;
		assert(NodeCounter.count(n1) > 0);
		assert(NodeCounter.count(n2) > 0);

		int n1_i = NodeCounter[n1];
		int n2_i = NodeCounter[n2];

		assert(n1_i < (int)NodeCounter.size()+1);
		assert(n2_i < (int)NodeCounter.size()+1);


		line_fields.push_back(lexical_cast<string>(n1_i));
		line_fields.push_back(lexical_cast<string>(n2_i));

		//if(e->label != "")
		//{
			line_fields.push_back(e->label);

		//}

		line_fields.push_back(lexical_cast<string>(e->pgf_protect));

			
		linesForEdges.push_back(boost::join(line_fields, separatorForSerialization));
	}

	ofstream output;
	output.open (filename.c_str(), ios::out | ios::trunc);
	if (output.is_open())
	{
		output << "CODE:" << "\n";
		output << boost::join(linesFromCode, "\n") << "\n";
		output << "NODES:" << "\n";
		output << boost::join(linesForNodes, "\n") << "\n";
		output << "EDGES:" << "\n";
		output << boost::join(linesForEdges, "\n");
	}
	else
	{
		throw std::runtime_error("Cannot open output file for graph serialization: "+filename);
	}
}

void Graph::readFromFile(string filename)
{
    using boost::lexical_cast;
    using boost::bad_lexical_cast;

	vector<string> linesForCode;
	vector<string> linesForNodes;
	vector<string> linesForEdges;

	string problematic_part = "|||||||";
	string subsitute_problem = "|||SLASH|||";
	string substitue_indicator = "SLASH";

	int mode = -1;
	ifstream graphStr;
	graphStr.open (filename.c_str(), ios::in);
	if(graphStr.is_open())
	{
		string line;
		int lineCounter = 0;
		while(graphStr.good())
		{
			lineCounter++;

			if((lineCounter % 100000) == 0)
				cout << "\r" << lineCounter << flush;
			
			getline (graphStr, line);
			Utilities::eraseNL(line);

			if(line.length() == 0)
				continue;

			if(line.find(problematic_part) != string::npos)
			{
				line.replace(line.find(problematic_part),problematic_part.length(),subsitute_problem);
			}
			
			if(line == "CODE:")
			{
				mode = 1;
			}
			else if(line == "NODES:")
			{
				mode = 2;
			}
			else if(line == "EDGES:")
			{
				mode = 3;
			}
			else
			{
				assert(mode > 0);
				if(mode == 1)
				{
					linesForCode.push_back(line);
				}
				else if(mode == 2)
				{
					linesForNodes.push_back(line);
				}
				else if(mode == 3)
				{
					linesForEdges.push_back(line);
				}
				else
				{
					throw std::runtime_error("Undefined mode state.");
				}
			}
		}
		graphStr.close();
	}
	else
	{
		throw std::runtime_error("Cannot open graph file: "+filename);
	}
		
	LocusCodeAllocation localCODE;
	localCODE.readFromVector(linesForCode);

	map<int, Node*> idx2Node;
	map<int, Edge*> idx2Edge;
	for(unsigned int i = 0; i < linesForNodes.size(); i++)
	{
		if((i % 100000) == 0)
			cout << "\r" << i << "/" << linesForNodes.size() << "[N]" << flush;
		
		string line = linesForNodes.at(i);
		vector<string> fields;
		boost::iter_split(fields, line, boost::first_finder(separatorForSerialization, boost::is_iequal()));
		if(fields.size() != 3)
		{
			throw std::runtime_error("Cannot node-parse this line (expext 6 fields): "+line);
		}
		int n_idx = lexical_cast<int>(fields.at(0));
		int level = lexical_cast<int>(fields.at(1));
		bool terminal = lexical_cast<bool>(fields.at(2));

		Node* n = new Node();
		n->level = level;
		n->terminal = terminal;
		registerNode(n, level);

		idx2Node[n_idx] = n;
	}
	
	bool encounteredErr = false;
	for(unsigned int i = 0; i < linesForEdges.size(); i++)
	{
		if((i % 100000) == 0)
			cout << "\r" << i << "/" << linesForEdges.size() << "[E]" << flush;
		
		string line = linesForEdges.at(i);
		vector<string> fields;
		boost::iter_split(fields, line, boost::first_finder(separatorForSerialization, boost::is_iequal()));
		if((fields.size() != 6) &&  (fields.size() != 8))
		{
			throw std::runtime_error("Cannot edge-parse this line (expext 6/8 fields): "+line);
		}
		int e_idx = lexical_cast<int>(fields.at(0));
		string locusID = fields.at(1);
		double count = lexical_cast<double>(fields.at(2));
		if(fields.at(3) == substitue_indicator)
		{
			fields.at(3) = "|";
		}
		
		unsigned char emission;
        try
        {
        	emission = lexical_cast<unsigned char>(fields.at(3));
        }
        catch(bad_lexical_cast &)
        {
            cerr << "Cannot cast to unsigned char: " << fields.at(3) << "--" << line << "\n";
            for(int i = 0; i < (int)line.length(); i++)
            {
            	cerr << i << "-" <<  line.at(i) << "-" << int(line.at(i)) << "\n";
            }
            encounteredErr = true;
            
            exit(1);
        }
        

		int from_idx;
		int to_idx;

		try {
		from_idx = lexical_cast<int>(fields.at(4));
		to_idx = lexical_cast<int>(fields.at(5));
		}
		catch(...)
		{
			cout << "Exception at line\n" << line << "\n\n";
			exit(1);
		}

		string label = "";
		bool protected_pgf = false;
		
		if(fields.size()>6)
		{
			label = fields.at(6);
	        try
	        {
	        	protected_pgf = lexical_cast<bool>(fields.at(7));
	        }
	        catch(bad_lexical_cast &)
	        {
	            cerr << "Cannot cast to bool: " << fields.at(7) << "--" << line << "\n";
	            for(int i = 0; i < (int)line.length(); i++)
	            {
	            	cerr << i << "-" <<  line.at(i) << "-" << int(line.at(i)) << "\n";
	            }
	            
				//encounteredErr = true;

				protected_pgf = false;
				
	            //exit(1);
	        }
		
		}

		Edge* e = new Edge();
		e->count = count;
		std::string emissionString = localCODE.deCode(locusID, emission);
		assert(emissionString.length() == 1);
		e->emission = emissionString.at(0);
		e->locus_id = locusID;
		e->label = label;
		e->pgf_protect = protected_pgf;
		idx2Edge[e_idx] = e;
		registerEdge(e);

		assert(idx2Node.count(from_idx) > 0);
		if(!(idx2Node.count(to_idx) > 0))
		{
			cout << "Node index " << to_idx << " is unknown\n";
			cout << "Edge line: " << line << "\n\n" << flush;
		}
		assert(idx2Node.count(to_idx) > 0);

		Node* n1 = idx2Node[from_idx];
		Node* n2 = idx2Node[to_idx];

		e->From = n1;
		e->To = n2;

		n1->Outgoing_Edges.insert(e);
		n2->Incoming_Edges.insert(e);
	}

	if(encounteredErr)
	{
		cerr << "Error!";
		exit(1);
	}
	
	//TODO reactivate
	/*
	cout << "\rCONSISTENCY...             " << flush;
	checkConsistency(false);
	*/

	cout << "\n\n" << flush;
	
	filename_last_read = filename;
}



std::vector<std::string> Graph::readGraphLoci(std::string graphDir)
{
	std::vector<std::string> forReturn;

	std::string segmentsFile = graphDir + "/PRG/segments.txt";

	ifstream segmentsStream;
	segmentsStream.open(segmentsFile.c_str());
	if(! segmentsStream.is_open())
	{
		throw std::runtime_error("Cannot open segments file: "+segmentsFile);
	}

	std::vector<std::string> individualSegmentFiles;
	std::string line;
	while(segmentsStream.good())
	{
		std::getline(segmentsStream, line);
		Utilities::eraseNL(line);
		if(line.length())
		{
			individualSegmentFiles.push_back(graphDir+"/PRG/"+line);
		}
	}
	segmentsStream.close();

	for(unsigned int segmentI = 0; segmentI < individualSegmentFiles.size(); segmentI++)
	{
		std::string F = individualSegmentFiles.at(segmentI);
		ifstream thisSegmentStream;
		thisSegmentStream.open(F.c_str());
		if(! thisSegmentStream.is_open())
		{
			throw std::runtime_error("Cannot open one segment file: "+F);
		}
		assert(thisSegmentStream.good());
		std::string firstLine;
		std::getline(thisSegmentStream, firstLine);
		Utilities::eraseNL(firstLine);

		std::vector<std::string> firstLine_fields = Utilities::split(firstLine, " ");

		for(unsigned int fI = 1; fI < firstLine_fields.size(); fI++)
		{
			forReturn.push_back(firstLine_fields.at(fI));
		}

		thisSegmentStream.close();
	}

	return forReturn;
}
