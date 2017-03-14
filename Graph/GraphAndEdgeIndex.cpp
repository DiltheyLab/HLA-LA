/*
 * GraphAndEdgeIndex.cpp
 *
 *  Created on: 29.07.2013
 *      Author: AlexanderDilthey
 */

#include "GraphAndEdgeIndex.h"

#include "../Utilities.h"

int pointerStrintLength = -1;

vector<kMerInfo> forwardScan(Node* start, int limit, int firstEdgeGap);
vector<kMerInfo> forwardScanRec(Node* currentNode, int depth, int realdepth, int limit, int firstEdgeGap);


GraphAndEdgeIndex::GraphAndEdgeIndex(Graph* graph, int k) {
	kMerSize = k;
	g = graph;
	if(g != 0)
	{
		Index();
	}
}

std::vector<std::string> GraphAndEdgeIndex::getIndexedkMers()
{
	assert(g != 0);
	std::vector<std::string> forReturn;
	for(std::map<std::string, std::vector<kMerInGraphSpec> >::iterator kMerIt = kMers.begin(); kMerIt != kMers.end(); kMerIt++)
	{
		std::string thiskMer = kMerIt->first;
		forReturn.push_back(thiskMer);
	}
	return forReturn;
}

std::vector<kMerEdgeChain*> GraphAndEdgeIndex::findChains(std::string sequence)
{
	assert(g != 0);
	std::vector<kMerEdgeChain> forReturn;
	std::vector<kMerEdgeChain> runningChains;
	bool superquiet = true;

	if(! superquiet)
		std::cout  << Utilities::timestamp() << "GraphAndEdgeIndex::findChains(..): Enter function.\n" << std::flush;

	// std::cout << "Search for sequence " << sequence << "\n";

	std::vector<std::string> kMers = Utilities::partitionStringIntokMers(sequence, kMerSize);

	if(kMers.size() == 0)
	{
		std::vector<kMerEdgeChain*> forReturn_empty;
		return forReturn_empty;
	}
	
	bool verbose = false;
	
	if(! superquiet)
		std::cout  << Utilities::timestamp() << "GraphAndEdgeIndex::findChains(..): Process " << kMers.size() << " kMers.\n" << std::flush;

	size_t positions_unknown_kMers = 0;

	auto archiveChain = [&](kMerEdgeChain chain) -> void {
		forReturn.push_back(chain);
	};

	std::string firstMer = kMers.at(0);
	std::vector<kMerInGraphSpec> firstMer_positions = queryIndex(firstMer);
	for(unsigned int firstMerI = 0; firstMerI < firstMer_positions.size(); firstMerI++)
	{
		kMerEdgeChain firstMerChain;
		firstMerChain.sequence_begin = 0;
		firstMerChain.sequence_end = kMerSize - 1;
		firstMerChain.traversedEdges = firstMer_positions.at(firstMerI).traversedEdges;
		runningChains.push_back(firstMerChain);
	}

	for(unsigned int seqI = kMerSize; seqI < sequence.size(); seqI++)
	{
		std::string seqCharacter = sequence.substr(seqI, 1);
	
		if((! superquiet) && ((seqI % 100000) == 0))
			std::cout << "\r" << "seqI " << seqI << " / " << sequence.size() << ": Search character " << seqCharacter << "\n" << std::flush;

		if(verbose)
			std::cout << "seqI " << seqI << " / " << sequence.size() << ": Search character " << seqCharacter << "\n" << std::flush;


		// extend running chains
		for(int chainI = (runningChains.size() - 1); chainI >= 0; chainI--)
		{
			kMerEdgeChain& chain = runningChains.at(chainI);
			Node* targetNode = chain.traversedEdges.back()->To;
			
			if(verbose)
			{
				std::cout << "\t" << "Process chain " << chainI << "/" <<  runningChains.size() << " [" << (kMerEdgeChain*)(&chain) << "], ending in targetNode " << targetNode << "\n" << std::flush;
				std::cout << "\t\tChains still in: ";
				for(unsigned int cI2 = 0; cI2 < runningChains.size(); cI2++)
				{
					std::cout << (kMerEdgeChain*)(&(runningChains.at(cI2))) << " ";
				}
				std::cout << "\n";
			}

//			if(chain.traversedEdges.at(0)->From->level == 77)
//			{
//				std::cerr << "Chain starting node " << chain.traversedEdges.at(0)->From << ": Enter extension.\n";
//			}
			while(nodes_jumpOverGaps.count(targetNode))
			{

				std::vector<Edge*> moreTraversedEdges = nodes_jumpOverGaps.at(targetNode);
				chain.traversedEdges.insert(chain.traversedEdges.end(), moreTraversedEdges.begin(), moreTraversedEdges.end());
				targetNode = moreTraversedEdges.back()->To;
				
				if(verbose)
				{
					std::cout << "\t" << "Had to jump over target nodes, new final node " << targetNode << "\n" << std::flush;
				}
				
//				if(chain.traversedEdges.at(0)->From->level == 77)
//				{
//					std::cerr << "Chain starting node " << chain.traversedEdges.at(0)->From << ": Jump over gaps.\n";
//				}
			}
			
			std::vector<std::vector<Edge*>> compatibleEdgeSequences;
			std::set<Edge*> availableStartEdges = targetNode->Outgoing_Edges;

			auto forwardScan = [&](Edge* startEdge) -> std::vector<std::vector<Edge*> > {
				std::vector<std::vector<Edge*> > foundEdgePaths;
				std::vector<std::vector<Edge*> > runningEdges;

				std::vector<Edge*> firstV;
				firstV.push_back(startEdge);

				runningEdges.push_back(firstV);

				while(runningEdges.size() > 0)
				{
					for(int eI = (runningEdges.size() - 1); eI >= 0; eI--)
					{
						std::vector<Edge*>& edgeSequence = runningEdges.at(eI);
						Edge* tipEdge = edgeSequence.at(edgeSequence.size() - 1);
						std::string tipEdgeEmission = tipEdge->getEmission();
						if(tipEdgeEmission == "_")
						{
							Node* tipNode = tipEdge->To;
							while(nodes_jumpOverGaps.count(tipNode))
							{
								std::vector<Edge*> moreTraversedEdges = nodes_jumpOverGaps.at(targetNode);
								edgeSequence.insert(edgeSequence.end(), moreTraversedEdges.begin(), moreTraversedEdges.end());
								tipNode = moreTraversedEdges.back()->To;
							}
							std::set<Edge*> nextEdges = tipNode->Outgoing_Edges;
							if(nextEdges.size() == 0)
							{
								runningEdges.erase(runningEdges.begin() + eI);
							}
							else if(nextEdges.size() == 1)
							{
								edgeSequence.push_back(*(nextEdges.begin()));
							}
							else
							{
								assert(nextEdges.size() > 1);
								std::vector<Edge*> templateToCopy = edgeSequence;
								for(std::set<Edge*>::iterator nextEdgeIt = nextEdges.begin(); nextEdgeIt != nextEdges.end(); nextEdgeIt++)
								{
									Edge* nextEdge = *nextEdgeIt;
									if(nextEdgeIt == nextEdges.begin())
									{
										edgeSequence.push_back(*nextEdgeIt);
									}
									else
									{
										std::vector<Edge*> newEdgeSequence = templateToCopy;
										newEdgeSequence.push_back(*nextEdgeIt);
										runningEdges.push_back(newEdgeSequence);
									}
								}
							}
						}
						else
						{
							foundEdgePaths.push_back(edgeSequence);
							runningEdges.erase(runningEdges.begin() + eI);
						}
					}
				}

				return foundEdgePaths;
			};


			for(std::set<Edge*>::iterator availableStartEdgeIt = availableStartEdges.begin(); availableStartEdgeIt != availableStartEdges.end(); availableStartEdgeIt++)
			{
				Edge* availableStartEdge = *availableStartEdgeIt;
				std::vector<std::vector<Edge*> > followingPaths = forwardScan(availableStartEdge);
				
				if(verbose)
				{
					std::cout << "\t\t" << "Scan forward paths from edge " << availableStartEdge << " [coming from node " << availableStartEdge->From << "]" << "\n" << std::flush;
				}
				
				for(unsigned int fI = 0; fI < followingPaths.size(); fI++)
				{
					std::vector<Edge*>& impliedEdgeSequence = followingPaths.at(fI);
					Edge* lastEdge = impliedEdgeSequence.at(impliedEdgeSequence.size() - 1);
					std::string edgeEmission = lastEdge->getEmission();
					
					
					if(verbose)
					{
						std::cout << "\t\t\t" << "Alternative " << fI << "/" << followingPaths.size() << ": path length " << impliedEdgeSequence.size() << " with symbol: " << edgeEmission << " [looking for " << seqCharacter << "]" << "\n" << std::flush;
					}
									
					assert(edgeEmission != "_");
//					if(chain.traversedEdges.at(0)->From->level == 77)
//					{
//						std::cerr << "Chain starting node " << chain.traversedEdges.at(0)->From << ": fI = " << fI << ", edgeEmission = " << edgeEmission << ", seqCharacter = " << seqCharacter << "\n" << std::flush;
//					}
					if(edgeEmission == seqCharacter)
					{
						compatibleEdgeSequences.push_back(impliedEdgeSequence);
					}
				}
			}

//			if(chain.traversedEdges.at(0)->From->level == 77)
//			{
//				std::cerr << "Chain starting node " << chain.traversedEdges.at(0)->From << ": compatibleEdgeSequences.size(): " <<  compatibleEdgeSequences.size() << "\n" << std::flush;
//			}



			if(compatibleEdgeSequences.size() == 0)
			{
				assert(&(chain) == &(runningChains.at(chainI)));
				if(verbose)
				{				
					// std::cout << "\n\t\t\tChains still in: ";
					// for(unsigned int cI2 = 0; cI2 < runningChains.size(); cI2++)
					// {
						// std::cout << (kMerEdgeChain*)(&(runningChains.at(cI2))) << " ";
					// }
					// std::cout << "\n";		
					std::cout << "\t\t Get rid of chainI " << chainI << "/" << (kMerEdgeChain*)(&chain) << "\n" << std::flush;

				}
								
				archiveChain(chain);
				std::vector<kMerEdgeChain>::iterator toDeleteIt = runningChains.begin() + chainI;
				
				assert(&(*toDeleteIt) == &(chain));
				
				runningChains.erase(runningChains.begin() + chainI);
				
				// if(verbose)
				// {				
					// std::cout << "\t\t\tChains still in: ";
					// for(unsigned int cI2 = 0; cI2 < runningChains.size(); cI2++)
					// {
						// std::cout << (kMerEdgeChain*)(&(runningChains.at(cI2))) << " ";
					// }
					// std::cout << "\n";					
				// }				
				
			}
			else if(compatibleEdgeSequences.size() == 1)
			{
				chain.sequence_end = seqI;
				chain.traversedEdges.insert(chain.traversedEdges.end(), compatibleEdgeSequences.at(0).begin(), compatibleEdgeSequences.at(0).end());
			}
			else
			{
				assert(compatibleEdgeSequences.size() > 1);
				kMerEdgeChain templateToCopy = chain;
				for(std::vector<std::vector<Edge*>>::iterator compatibleEdgeIt = compatibleEdgeSequences.begin(); compatibleEdgeIt != compatibleEdgeSequences.end(); compatibleEdgeIt++)
				{
					if(compatibleEdgeIt == compatibleEdgeSequences.begin())
					{
						chain.sequence_end = seqI;
						chain.traversedEdges.insert(chain.traversedEdges.end(), compatibleEdgeIt->begin(), compatibleEdgeIt->end());
					}
					else
					{
						kMerEdgeChain newChain = templateToCopy;
						newChain.sequence_end = seqI;
						newChain.traversedEdges.insert(newChain.traversedEdges.end(), compatibleEdgeIt->begin(), compatibleEdgeIt->end());
						runningChains.push_back(newChain);
					}
				}
			}
		}

		// kMers ending in this positions
		int kMerEndingIndex = seqI - kMerSize + 1;
		assert(kMerEndingIndex >= 0);
		assert(kMerEndingIndex < (int)kMers.size());
		std::string kMerEndingHere = kMers.at(kMerEndingIndex);
		std::vector<kMerInGraphSpec> kMer_positions = queryIndex(kMerEndingHere);

		for(unsigned int positionI = 0; positionI < kMer_positions.size(); positionI++)
		{
			// is this kMer hit represented by an existing chain?

			kMerInGraphSpec& kMerPossiblePosition = kMer_positions.at(positionI);
			bool represented = false;
			for(int chainI = (runningChains.size() - 1); chainI >= 0; chainI--)
			{
				kMerEdgeChain& existingChain = runningChains.at(chainI);
				if(existingChain.traversedEdges.back()->To == kMerPossiblePosition.traversedEdges.back()->To)
				{
					if(existingChain.traversedEdges.size() >= kMerPossiblePosition.traversedEdges.size())
					{
						int offset = kMerPossiblePosition.traversedEdges.size();
						assert(kMerPossiblePosition.traversedEdges.size() > 0);
						if(existingChain.traversedEdges.at(existingChain.traversedEdges.size() - offset)->From == kMerPossiblePosition.traversedEdges.at(0)->From)
						{
							represented = true;
							break;
						}
					}
				}
			}

			if(! represented)
			{
				kMerEdgeChain newChain;
				newChain.sequence_begin = seqI - kMerSize + 1;
				newChain.sequence_end = seqI;
				newChain.traversedEdges = kMerPossiblePosition.traversedEdges;
				runningChains.push_back(newChain);
			}
		}
	}

	// add in remaining existing chains
	for(unsigned int chainI = 0; chainI < runningChains.size(); chainI++)
	{
		archiveChain(runningChains.at(chainI));
	}

	std::vector<kMerEdgeChain*> forReturn_pointers;
	for(unsigned int chainI = 0; chainI < forReturn.size(); chainI++)
	{
		kMerEdgeChain* pointerChain = new kMerEdgeChain(forReturn.at(chainI));
		forReturn_pointers.push_back(pointerChain);
	}
	return forReturn_pointers;
}

void GraphAndEdgeIndex::fillEdgeJumper()
{
	assert(g != 0);
	class runningGapChain {
	public:
		std::vector<Edge*> traversedEdges;
		Node* targetNode;
	};

	std::map<Node*, runningGapChain> runningGapChains;

	auto archiveChain = [&](runningGapChain c) -> void {
		if(c.traversedEdges.size() > 1)
		{
			Node* firstNode = c.traversedEdges.front()->From;
			assert(nodes_jumpOverGaps.count(firstNode) == 0);
			nodes_jumpOverGaps[firstNode] = c.traversedEdges;
		}
	};
	auto nodeExtendsWithSimpleGap = [&](Node* n) -> Edge* {
		if(n->Outgoing_Edges.size() == 1)
		{
			Edge* e = *(n->Outgoing_Edges.begin());
			std::string emission = e->getEmission();
			assert(emission.length() == 1);
			if(emission == "_")
			{
				return e;
			}
		}
		return 0;
	};
	int levels = g->NodesPerLevel.size();
	for(int level = 0; level < levels; level++)
	{
		std::set<Node*> nodesAtLevel_set = g->NodesPerLevel.at(level);
		std::vector<Node*> nodesAtLevel(nodesAtLevel_set.begin(), nodesAtLevel_set.end());
		for(unsigned int nodeI = 0; nodeI < nodesAtLevel.size(); nodeI++)
		{
			Node* n = nodesAtLevel.at(nodeI);
			Edge* simpleExtensionEdgeIfExists = nodeExtendsWithSimpleGap(n);
			if(simpleExtensionEdgeIfExists != 0)
			{
				Node* targetNode = simpleExtensionEdgeIfExists->To;
				runningGapChain chain;
				if(runningGapChains.count(n))
				{
					chain = runningGapChains.at(n);
					runningGapChains.erase(n);
				}
				chain.traversedEdges.push_back(simpleExtensionEdgeIfExists);
				chain.targetNode = targetNode;
				runningGapChains[targetNode] = chain;
			}
			else
			{
				if(runningGapChains.count(n))
				{
					archiveChain(runningGapChains.at(n));
					runningGapChains.erase(n);
				}
			}
		}
	}
	for(std::map<Node*, runningGapChain>::iterator remainingChainIt = runningGapChains.begin(); remainingChainIt != runningGapChains.end(); remainingChainIt++)
	{
		archiveChain(remainingChainIt->second);
	}
}

void GraphAndEdgeIndex::Index()
{
	assert(g != 0);
	bool quiet = true;
	bool superQuiet = false;

	int nonQuiet_levelModulo = 1;

	std::ostringstream nullPointerString;
	nullPointerString << setw(15) << (void const *) 0;
	int pointerStrintLength = nullPointerString.str().length();

	int levels = g->NodesPerLevel.size();
	vector<string> loci = g->getAssignedLoci();

	set<Node*> nodesLastLevel = g->NodesPerLevel[levels-1];
	string lastLocusID = loci.at(loci.size()-1);
	Node* classicalN0 = *(g->NodesPerLevel.at(0).begin());

	std::map<std::string, std::map<std::string, std::pair<int, int> > > _kMerPositions;
	auto kMer_bookkeeping = [&](std::string kMer, int firstLevel, int lastLevel, std::vector<Edge*>& traversedEdges) -> void {
		assert((int)kMer.length() == this->kMerSize);
		kMerInGraphSpec kM;
		kM.traversedEdges = traversedEdges;
		kMers[kMer].push_back(kM);
	};

	// make sure that last graph level has no gaps and that all edges are correctly labelled
	for(set<Node*>::iterator nodeIt = nodesLastLevel.begin(); nodeIt != nodesLastLevel.end(); nodeIt++)
	{
		set<Edge*> nodeIncomingEdges = (*nodeIt)->Incoming_Edges;
		for(set<Edge*>::iterator edgeIt = nodeIncomingEdges.begin(); edgeIt != nodeIncomingEdges.end(); edgeIt++)
		{
			Edge* edge = *edgeIt;
			assert(lastLocusID == edge->locus_id);
			string emission = edge->getEmission();
			assert(emission != "_");
		}
	}

	map<Node*, vector< kMerAtNode > > originalGraphAssignedKMers;


	for(int level = 0; level <= (levels - 1 - kMerSize); level++)
	{
		if(! superQuiet)
		{
			if((level % 10000) == 0)
			{
				std::cout << "\r" << level << " / " << (levels - 1 - kMerSize) << std::flush;
			}
		}

		int min_span_level = -1;
		int max_span_level = -1;

		if((! quiet) && ((level % nonQuiet_levelModulo) == 0))
		{
			cout << "Level " << level << "/" << (levels - 1 - kMerSize) << " of original graph\n";
			Node* firstN = *(g->NodesPerLevel.at(level).begin());
			if(firstN->Outgoing_Edges.size() > 0)
			{
				Edge* e = *(firstN->Outgoing_Edges.begin());
				cout << "\t\t" << e->locus_id << "\t\t" << e->label << "\n";
			}
		}
		string locusID = loci.at(level);

		map<Node*, map<string, vector<Edge*> > > edgeTargetCache;
		map<Edge*, kMerInfo> newEdgeNodeInfos;

		if(level == 0)
		{
			// Level 0: Forward-Scan nach kMers der Laenge x

			vector<kMerInfo> startFromN0 = forwardScan(classicalN0, kMerSize, -1);

			if(! quiet)
				cout << "\tFound " << startFromN0.size() << " initial kMers with no gap \n";

			for(vector<kMerInfo>::iterator kMerIt = startFromN0.begin(); kMerIt != startFromN0.end(); kMerIt++)
			{
				// Fuer alle gefundenen kMers:
				//	- erzeuge eine Edge
				//	- lege einen Eintrag an im Cache, der ueber die zu erzeugenden Nodes bestimmt
				//		dabei verwenden wir die Target-Node und die Memory-Adressen der letzten x-1 Edges als Index

				string kMer_string = Utilities::join(kMerIt->kMer_deCoded, "");
				assert((int)kMer_string.size() == kMerSize);
				assert((int)kMerIt->kMer_deCoded.size() == kMerSize);

				int firstLevel = kMerIt->traverseEdges.at(0)->From->level;
				int lastLevel = kMerIt->traverseEdges.at(kMerIt->traverseEdges.size()-1)->To->level;
				kMer_bookkeeping(kMer_string, firstLevel, lastLevel, kMerIt->traverseEdges);

				Edge* kMerEdge = new Edge();
				assert(kMerIt->traverseEdges.at(0)->getEmission() != "_");

				string km1Mer_index_string = kMerIt->traverseEdges_string.substr(pointerStrintLength);
				edgeTargetCache[kMerIt->traverseEdges.back()->To][km1Mer_index_string].push_back(kMerEdge);
				newEdgeNodeInfos[kMerEdge] = *kMerIt;

				generated_edges.insert(kMerEdge);
			}

			// Level 0: Forward-Scan nach kMers der Laenge x-1 mit Gaps vorne

			vector<kMerInfo> startFromN0_GAP = forwardScan(classicalN0, kMerSize-1, 1);

			if(! quiet)
				cout << "\tFound " << startFromN0_GAP.size() << " initial kMers with gap \n";

			for(vector<kMerInfo>::iterator kMerIt = startFromN0_GAP.begin(); kMerIt != startFromN0_GAP.end(); kMerIt++)
			{
				// Fuer alle gefundenen kMers:
				//	- erzeuge eine Edge
				//	- lege einen Eintrag an im Cache, der ueber die zu erzeugenden Nodes bestimmt
				//		dabei verwenden wir die Target-Node und die Memory-Adressen der letzten x-1 Edges als Index

				string kMer_string = Utilities::join(kMerIt->kMer_deCoded, "");
				assert(kMerIt->kMer_deCoded.size() == (kMerSize -1));
				assert((int)kMer_string.size() == kMerSize-1);

				string newLocusID = "L"+Utilities::ItoStr(level);

				// cout << kMer_string << "\n" << flush;

				Edge* kMerEdge = new Edge();
				kMerIt->gapEdge = true;

				assert(kMerIt->traverseEdges.at(0)->getEmission() == "_");

				string km1Mer_index_string = kMerIt->traverseEdges_string.substr(pointerStrintLength);

				edgeTargetCache[kMerIt->traverseEdges.back()->To][km1Mer_index_string].push_back(kMerEdge);
				newEdgeNodeInfos[kMerEdge] = *kMerIt;

				generated_edges.insert(kMerEdge);
			}
		}
		else
		{
			// alle Nodes dieses Levels
			vector<Node*> originalNodes = vector<Node*>(g->NodesPerLevel.at(level).begin(), g->NodesPerLevel.at(level).end());
			int node_c = 0;
			int node_c_total = originalNodes.size();

			for(vector<Node*>::iterator originalNodeIt = originalNodes.begin(); originalNodeIt != originalNodes.end(); originalNodeIt++)
			{
				node_c++;

				// potentielle Kanten (attachte KMers) durchgehen
				Node* originalNode = *originalNodeIt;
				assert(originalNode->Incoming_Edges.size() > 0);

				// only if next assertion fails...
				if(!(originalGraphAssignedKMers.count(originalNode) > 0))
				{
					cout << "\t Level " << level << ", search for originalNode " << originalNode << " information with result " << (originalGraphAssignedKMers.count(originalNode) > 0) << "\n";
					Edge* edgesToPrevious = *(originalNode->Incoming_Edges.begin());
					Node* previousNode = edgesToPrevious->From;

					vector<kMerInfo> startFromLastNode = forwardScan(previousNode, kMerSize, 0);

					cout << "\t Found " << startFromLastNode.size() << " paths which should traverse our dubious node, taking edge " << edgesToPrevious << "\n";

					for(vector<kMerInfo>::iterator kMerIt = startFromLastNode.begin(); kMerIt != startFromLastNode.end(); kMerIt++)
					{
						string kMer_string = Utilities::join(kMerIt->kMer_deCoded, "");
						cout << "\t\t kMer_string " << kMer_string << "\n";
						cout << "\t\t traversed Edges" << kMerIt->traverseEdges_string << "\n";
					}
				}

				assert(originalGraphAssignedKMers.count(originalNode) > 0);

				int attachedKMer_c = 0;
				int attachedKMer_c_total = originalGraphAssignedKMers[originalNode].size();
				for(vector< kMerAtNode >::iterator kMerBasisEdgeIt = originalGraphAssignedKMers[originalNode].begin();  kMerBasisEdgeIt != originalGraphAssignedKMers[originalNode].end(); kMerBasisEdgeIt++)
				{
					attachedKMer_c++;
					if((! quiet) && ((level % nonQuiet_levelModulo) == 0))
					{
						cout << "\r                                                            ";
						cout << "\r\t Node " << node_c << "/" << node_c_total << ", attached kMer " << attachedKMer_c << "/" << attachedKMer_c_total << flush;
					}
					// fuer jede attachte Kante neue Kanten erstellen, die sich von der letzten Node
					// im Original-Graph ergeben

					kMerAtNode BasisForNewEdges = *kMerBasisEdgeIt;

					assert(BasisForNewEdges.traverseEdges.at(0)->From == originalNode);

					// Wenn das erste Symbol des attachten kMers kein Gap ist, fuegen wir frohlich neue Kanten hinzu
					if(BasisForNewEdges.traverseEdges.at(0)->getEmission() != "_")
					{
						Node* lastNodeInOriginalGraph = BasisForNewEdges.traverseEdges.back()->To;

						// some debug information

						if(min_span_level == -1)
						{
							min_span_level = lastNodeInOriginalGraph->level;
						}
						if(max_span_level == -1)
						{
							max_span_level = lastNodeInOriginalGraph->level;
						}

						if((int)lastNodeInOriginalGraph->level > max_span_level)
						{
							max_span_level = lastNodeInOriginalGraph->level;
						}
						if((int)lastNodeInOriginalGraph->level < min_span_level)
						{
							min_span_level = lastNodeInOriginalGraph->level;
						}

						vector<kMerInfo> startFromLastNode = forwardScan(lastNodeInOriginalGraph, 1, 0);


						if(startFromLastNode.size() == 0)
						{
							if(! quiet)
								cout << " -- reached end of graph" << "\n";

							//assert((int)BasisForNewEdges.traverseEdges.back()->To->level == (int)levels);

							// Wir sind am Ende des normalen Graphen angelangt, haben aber noch zu wenige Positionen im kmer-Graph
							// Das kommt durch Gaps in den letzten X Symbolen
							// Wir fuegen Gap-KMers hinzu

							// Wenn das erste Symbol dieser potentiellen Kante ein Gap ist, muessen wir nur eine Gap-Edge ans Ende hinzufuegen
							// neues kMerInfo-Objekt, GAP!

							kMerInfo newNodeKMerInfo;
							newNodeKMerInfo.gapEdge = true;
							newNodeKMerInfo.kMer_coded = BasisForNewEdges.km1Mer_coded;
							newNodeKMerInfo.kMer_deCoded = BasisForNewEdges.km1Mer_deCoded;
							newNodeKMerInfo.p = 1;
							newNodeKMerInfo.traverseEdges = BasisForNewEdges.traverseEdges;
							newNodeKMerInfo.traverseEdges_string = BasisForNewEdges.traverseEdges_string;
							newNodeKMerInfo.allPGF = BasisForNewEdges.allPGF;

							string kMer_string = "_";

							Edge* kMerEdge = new Edge();

							BasisForNewEdges.lastNewNode->Outgoing_Edges.insert(kMerEdge);

							assert(newNodeKMerInfo.traverseEdges_string.length() == (newNodeKMerInfo.traverseEdges.size()*pointerStrintLength));
							string km1Mer_index_string = newNodeKMerInfo.traverseEdges_string.substr(pointerStrintLength);

							assert(g->Nodes.count(newNodeKMerInfo.traverseEdges.back()->To) > 0);

							edgeTargetCache[newNodeKMerInfo.traverseEdges.back()->To][km1Mer_index_string].push_back(kMerEdge);
							newEdgeNodeInfos[kMerEdge] = newNodeKMerInfo;

							generated_edges.insert(kMerEdge);

						}
						else
						{
							if((! quiet) && ((level % nonQuiet_levelModulo) == 0))
								cout << " -- " << startFromLastNode.size() << " expansions" << "\n";

							for(vector<kMerInfo>::iterator attachKMerIt = startFromLastNode.begin(); attachKMerIt != startFromLastNode.end(); attachKMerIt++)
							{

									assert(attachKMerIt->kMer_coded.size() == 1);
									assert(attachKMerIt->kMer_deCoded.size() == 1);

									// neues kMerInfo-Objekt
									kMerInfo newNodeKMerInfo;
									newNodeKMerInfo.gapEdge = false;
									newNodeKMerInfo.kMer_coded = BasisForNewEdges.km1Mer_coded;
									newNodeKMerInfo.kMer_deCoded = BasisForNewEdges.km1Mer_deCoded;
									newNodeKMerInfo.p = attachKMerIt->p;
									newNodeKMerInfo.traverseEdges = BasisForNewEdges.traverseEdges;
									newNodeKMerInfo.traverseEdges_string = BasisForNewEdges.traverseEdges_string;

									if(!(newNodeKMerInfo.kMer_deCoded.size() == (kMerSize-1)))
									{
										std::cerr << "!(newNodeKMerInfo.kMer_deCoded.size() == (kMerSize-1))\n";
										std::cerr << "newNodeKMerInfo.kMer_deCoded.size(): " << newNodeKMerInfo.kMer_deCoded.size() << "\n";
										std::cerr << "(kMerSize-1): " << (kMerSize-1) << "\n" << std::flush;
										std::cerr << "newNodeKMerInfo.kMer_deCoded: " << Utilities::join(newNodeKMerInfo.kMer_deCoded, " ") << "\n" << std::flush;

									}

									assert(newNodeKMerInfo.kMer_deCoded.size() == (kMerSize-1));


									newNodeKMerInfo.kMer_coded.push_back(attachKMerIt->kMer_coded.at(0));
									newNodeKMerInfo.kMer_deCoded.push_back(attachKMerIt->kMer_deCoded.at(0));
									newNodeKMerInfo.traverseEdges.insert(newNodeKMerInfo.traverseEdges.end(), attachKMerIt->traverseEdges.begin(), attachKMerIt->traverseEdges.end());
									newNodeKMerInfo.traverseEdges_string = newNodeKMerInfo.traverseEdges_string.append(attachKMerIt->traverseEdges_string);
									newNodeKMerInfo.allPGF = false; //(BasisForNewEdges.allPGF && attachKMerIt->allPGF);

									string kMer_string = Utilities::join(newNodeKMerInfo.kMer_deCoded, "");
									assert((int)kMer_string.size() == kMerSize);

									int firstLevel = newNodeKMerInfo.traverseEdges.at(0)->From->level;
									int lastLevel = newNodeKMerInfo.traverseEdges.at(newNodeKMerInfo.traverseEdges.size()-1)->To->level;
									kMer_bookkeeping(kMer_string, firstLevel, lastLevel, newNodeKMerInfo.traverseEdges);

									Edge* kMerEdge = new Edge();
									BasisForNewEdges.lastNewNode->Outgoing_Edges.insert(kMerEdge);

									for(int lI = 0; lI < (int)attachKMerIt->traverseEdges.size(); lI++)
									{
										Edge* traversedEdge = attachKMerIt->traverseEdges.at(lI);
										if(traversedEdge->getEmission() != "_")
										{
											// kMerEdge->levelsNucleotideGraph.push_back(traversedEdge->From->level);
										}
									}


									assert(newNodeKMerInfo.traverseEdges_string.length() == (newNodeKMerInfo.traverseEdges.size()*pointerStrintLength));
									string km1Mer_index_string = newNodeKMerInfo.traverseEdges_string.substr(pointerStrintLength);

									edgeTargetCache[attachKMerIt->traverseEdges.back()->To][km1Mer_index_string].push_back(kMerEdge);
									newEdgeNodeInfos[kMerEdge] = newNodeKMerInfo;

									generated_edges.insert(kMerEdge);
							}
						}
					}
					else
					{

						Node* lastNodeInOriginalGraph = BasisForNewEdges.traverseEdges.back()->To;

						if(min_span_level == -1)
						{
							min_span_level = lastNodeInOriginalGraph->level;
						}
						if(max_span_level == -1)
						{
							max_span_level = lastNodeInOriginalGraph->level;
						}
						if((int)lastNodeInOriginalGraph->level > max_span_level)
						{
							max_span_level = lastNodeInOriginalGraph->level;
						}
						if((int)lastNodeInOriginalGraph->level < min_span_level)
						{
							min_span_level = lastNodeInOriginalGraph->level;
						}

						// Wenn das erste Symbol dieser potentiellen Kante ein Gap ist, muessen wir nur eine Gap-Edge ans Ende hinzufuegen
						// neues kMerInfo-Objekt, GAP!

						if((! quiet) && ((level % nonQuiet_levelModulo) == 0))
							cout << " -- gap length " << BasisForNewEdges.traverseEdges.size() << "\n";

						kMerInfo newNodeKMerInfo;
						newNodeKMerInfo.gapEdge = true;
						newNodeKMerInfo.kMer_coded = BasisForNewEdges.km1Mer_coded;
						newNodeKMerInfo.kMer_deCoded = BasisForNewEdges.km1Mer_deCoded;
						newNodeKMerInfo.p = 1;
						newNodeKMerInfo.traverseEdges = BasisForNewEdges.traverseEdges;
						newNodeKMerInfo.traverseEdges_string = BasisForNewEdges.traverseEdges_string;

						string kMer_string = "_";

						Edge* kMerEdge = new Edge();
						BasisForNewEdges.lastNewNode->Outgoing_Edges.insert(kMerEdge);

						//vector<Edge*> traverseEdges_m1(newNodeKMerInfo.traverseEdges.begin()+1, newNodeKMerInfo.traverseEdges.end());
						//traverseEdges_m1.erase (traverseEdges_m1.begin(),traverseEdges_m1.begin()+1);

						assert(newNodeKMerInfo.traverseEdges_string.length() == (newNodeKMerInfo.traverseEdges.size()*pointerStrintLength));
						string km1Mer_index_string = newNodeKMerInfo.traverseEdges_string.substr(pointerStrintLength);

						edgeTargetCache[newNodeKMerInfo.traverseEdges.back()->To][km1Mer_index_string].push_back(kMerEdge);
						newEdgeNodeInfos[kMerEdge] = newNodeKMerInfo;

						generated_edges.insert(kMerEdge);

					}
				}
			}
		}

		originalGraphAssignedKMers.clear();

		int newNodeC = 0;
		int targetEdgeC = 0;

		for(map<Node*, map<string, vector<Edge*> > >::iterator targetNodeIt = edgeTargetCache.begin(); targetNodeIt != edgeTargetCache.end(); targetNodeIt++)
		{
			for( map< string, vector<Edge*> >::iterator targetStringIt = targetNodeIt->second.begin(); targetStringIt != targetNodeIt->second.end(); targetStringIt++)
			{
				// In der inneren Schleife haben wir alle Edges, die zur selben Node im original Graph fuehren und die
				// das auf den letzten x-1 Edges auf demselben Weg tun -- i.e. die in dieselbe Node fuehren sollen,
				// weil die darauffolgenden kMers identisch sind

				newNodeC++;

				string targetString = targetStringIt->first;
				vector<Edge*> attachedEdges = targetStringIt->second;
				assert(attachedEdges.size()>0);

				// Wir generieren eine Node und haengen sie an alle Kanten
				Node* kMerTargetNode = new Node();
				generated_nodes.insert(kMerTargetNode);

				kMerTargetNode->level = level+1;
				kMerTargetNode->terminal = (level == (levels - 1 - kMerSize)) ? true : false; // TODO korrekt?

				// Aus der ersten Edge, die wir haben, generieren wir die Info fuer zukuenftige kMer-Kanten,
				// die wir an die 2. Node des originals Graphs heften

				Edge* firstEdge = attachedEdges.at(0);
				assert(newEdgeNodeInfos.count(firstEdge) > 0);
				kMerInfo firstEdgeKMerInfo = newEdgeNodeInfos[firstEdge];
				assert((int)firstEdgeKMerInfo.kMer_deCoded.size() <= kMerSize);

				vector<unsigned char> km1Mer_coded = firstEdgeKMerInfo.kMer_coded;
				vector<string> km1Mer_deCoded  = firstEdgeKMerInfo.kMer_deCoded;

				//vector<Edge*> traverseEdges_m1 = firstEdgeKMerInfo.traverseEdges;
				vector<Edge*> traverseEdges_m1(firstEdgeKMerInfo.traverseEdges.begin()+1, firstEdgeKMerInfo.traverseEdges.end());
				string traverseEdges_m1_string = firstEdgeKMerInfo.traverseEdges_string.substr(pointerStrintLength);
				assert(traverseEdges_m1_string.length() == (traverseEdges_m1.size()*pointerStrintLength));

				if(firstEdgeKMerInfo.gapEdge == false)
				{
					km1Mer_coded.erase (km1Mer_coded.begin(),km1Mer_coded.begin()+1);
					km1Mer_deCoded.erase (km1Mer_deCoded.begin(),km1Mer_deCoded.begin()+1);
				}
				else
				{
					//assert(attachedEdges.size()==1);
				}


				kMerAtNode infoForNode2;
				infoForNode2.km1Mer_coded = km1Mer_coded;
				infoForNode2.km1Mer_deCoded = km1Mer_deCoded;
				infoForNode2.traverseEdges = traverseEdges_m1;
				infoForNode2.lastNewNode = kMerTargetNode;
				infoForNode2.traverseEdges_string = traverseEdges_m1_string;

				Node* originalGraphNode2 = firstEdgeKMerInfo.traverseEdges.at(0)->To;
				assert(g->Nodes.count(originalGraphNode2) > 0);

				if((! quiet) && ((level % nonQuiet_levelModulo) == 0) && (level > 6000000))
				{
					//cout << "\tFound originalGraphNode2 " << originalGraphNode2 << " with " << attachedEdges.size() << " attached edges!\n";
				}

				if(! quiet)
				{
					std::cout << "Attach to node " << originalGraphNode2 << " at level " << originalGraphNode2->level << " new kMerAtNode object!\n";
					std::cout << "\tkm1Mer_deCoded: " << Utilities::join(infoForNode2.km1Mer_deCoded, " ") << "\n";
					std::cout << "\ttraverseEdges.size(): " << infoForNode2.traverseEdges.size() << "\n";
					std::cout << "\tinfoForNode2.traverseEdges.back()->To->level: " << infoForNode2.traverseEdges.back()->To->level << "\n" << std::flush;
					std::cout << "\tfirstEdgeKMerInfo.gapEdge: " << firstEdgeKMerInfo.gapEdge << "\n" << std::flush;

				}


				originalGraphAssignedKMers[originalGraphNode2].push_back(infoForNode2);

				/*
				for(unsigned int idx = 0; idx < attachedEdges.size(); idx++)
				{
					Edge* e = attachedEdges.at(idx);
					Node* shouldbeNode2 = newEdgeNodeInfos[e].traverseEdges.at(0)->To;
					assert(shouldbeNode2 == originalGraphNode2);
				}
				*/

			}
		}

		assert(min_span_level <= max_span_level);

		string newLocusID = "L"+Utilities::ItoStr(level);
		if((! quiet) && ((level % nonQuiet_levelModulo) == 0))
		{
			cout << "\t\t added " << newNodeC << " nodes, " << targetEdgeC << " edges, and have " << kMers.size() << " encoded kMers \n";
			cout << "\t\t\t target node span original graph: " << min_span_level << " - " << max_span_level << "\n";
		}
	}

	if(! superQuiet)
	{
		std::cout << "\n" << std::flush;
	}


	for(std::set<Node*>::iterator nIt = generated_nodes.begin(); nIt != generated_nodes.end(); nIt++)
	{
		delete(*nIt);
	}

	// make sure kMer list is clean.
	for(std::map<std::string, std::vector<kMerInGraphSpec> >::iterator kMerIt = kMers.begin(); kMerIt != kMers.end(); kMerIt++)
	{
		std::string kMer = kMerIt->first;
		std::vector<kMerInGraphSpec>& graphPositions = kMerIt->second;
		for(unsigned int pI = 0; pI < (graphPositions.size() - 1); pI++)
		{
			for(unsigned int pI2 = pI+1; pI2 < graphPositions.size(); pI2++)
			{
				bool identical = true;
				kMerInGraphSpec& kMer_1_details = graphPositions.at(pI);
				kMerInGraphSpec& kMer_2_details = graphPositions.at(pI2);
				if(kMer_1_details.traversedEdges.size() != kMer_2_details.traversedEdges.size())
				{
					identical = false;
				}
				else
				{
					for(unsigned int edgeI = 0; edgeI < kMer_1_details.traversedEdges.size(); edgeI++)
					{
						if(kMer_1_details.traversedEdges.at(edgeI) != kMer_2_details.traversedEdges.at(edgeI))
						{
							identical = false;
							break;
						}
					}
				}
				assert(! identical);
			}
		}
	}
}


void GraphAndEdgeIndex::printIndex()
{
	assert(g != 0);
	std::cout << "GraphAndEdgeIndex::printIndex(): " << kMers.size() << " kMers.\n" << std::flush;

	for(std::map<std::string, std::vector<kMerInGraphSpec > >::iterator kMerPosIt = kMers.begin(); kMerPosIt != kMers.end(); kMerPosIt++)
	{
		std::string kMer = kMerPosIt->first;

		std::cout << "\tkMer " << kMer << "\n";

		std::vector< kMerInGraphSpec >& kMerPos = kMers.at(kMer);
		for(std::vector<kMerInGraphSpec >::iterator posIt = kMerPos.begin(); posIt != kMerPos.end(); posIt++)
		{
			std::vector<Edge*>& traversedEdges = posIt->traversedEdges;
			Node* firstNode = traversedEdges.front()->From;
			Node* lastNode = traversedEdges.back()->To;
			std::cout << "\t\tfrom node " << firstNode << " / " << firstNode->level << " to node " << lastNode << " / " << lastNode->level << "\n";
		}
	}

	std::cout << std::flush;
}

std::vector<kMerInGraphSpec> GraphAndEdgeIndex::queryIndex(std::string kMer)
{
	assert(g != 0);
	if(kMers.count(kMer) == 0)
	{
		return std::vector<kMerInGraphSpec>();
	}
	else
	{
		return kMers.at(kMer);
	}
}

GraphAndEdgeIndex::~GraphAndEdgeIndex()
{
//	for(std::set<Node*>::iterator nodeIt = generated_nodes.begin(); nodeIt != generated_nodes.end(); nodeIt++)
//	{
//		Node* n = *nodeIt;
//		delete(n);
//	}
//	generated_nodes.clear();

	for(std::set<Edge*>::iterator edgeIt = generated_edges.begin(); edgeIt != generated_edges.end(); edgeIt++)
	{
		Edge* e = *edgeIt;
		delete(e);
	}
	generated_edges.clear();
}

vector<kMerInfo> forwardScan(Node* start, int limit, int firstEdgeGap)
{
	if(pointerStrintLength == -1)
	{
		std::ostringstream nullPointerString;
		nullPointerString << setw(15) << (void const *) 0;
		pointerStrintLength = nullPointerString.str().length();
	}

	vector<kMerInfo> forReturn = forwardScanRec(start, 0, 0, limit, firstEdgeGap);
	for(vector<kMerInfo>::iterator kMerIt = forReturn.begin(); kMerIt != forReturn.end(); kMerIt++)
	{
		assert((int)(kMerIt->kMer_coded.size()) == limit);
		assert((int)(kMerIt->kMer_deCoded.size()) == limit);
		reverse(kMerIt->kMer_coded.begin(), kMerIt->kMer_coded.end());
		reverse(kMerIt->kMer_deCoded.begin(), kMerIt->kMer_deCoded.end());
		reverse(kMerIt->traverseEdges.begin(), kMerIt->traverseEdges.end());

		std::ostringstream o_traverseEdges_string;
		bool allPGF = true;
		for(int i = 0; i < (int)kMerIt->traverseEdges.size(); i++)
		{
			o_traverseEdges_string << setw(15) << (void const * ) kMerIt->traverseEdges.at(i);
			allPGF = (allPGF && (bool)kMerIt->traverseEdges.at(i)->pgf_protect);
		}
		kMerIt->allPGF = allPGF;

		string traverseEdges_string = o_traverseEdges_string.str();
		if(traverseEdges_string.length() != (kMerIt->traverseEdges.size()*pointerStrintLength))
		{
			cout << "One pointer as string: " << setw(15) << (void const * ) kMerIt->traverseEdges.at(0) << "\n";
			cout << "kMerIt->traverseEdges.size(): " << kMerIt->traverseEdges.size() << "\n";
			cout << "pointerStrintLength: " << pointerStrintLength << "\n";
			cout << "traverseEdges_string.length(): " << traverseEdges_string.length() << "\n";
			cout << "traverseEdges_string: " << traverseEdges_string << "\n";
		}
		assert(traverseEdges_string.length() == (kMerIt->traverseEdges.size()*pointerStrintLength));


		kMerIt->traverseEdges_string = traverseEdges_string;

	}
	return forReturn;
}

vector<kMerInfo> forwardScanRec(Node* currentNode, int depth, int realdepth, int limit, int firstEdgeGap)
{
	vector<kMerInfo> returnMers;

	if(limit == depth)
	{
		kMerInfo i;
		i.p = 1;
		i.gapEdge = false;
		returnMers.push_back(i);
		return returnMers;
	}
	else
	{
		for(set<Edge*>::iterator eIt = currentNode->Outgoing_Edges.begin(); eIt != currentNode->Outgoing_Edges.end(); eIt++)
		{
			Edge* e = *eIt;
			Node* targetNode = e->To;

			// TODO maybe we have such nodes, then we need to ignore them
			assert(currentNode->Sum_Outgoing() != 0);
			double edgeP = e->count/currentNode->Sum_Outgoing();
			assert(edgeP >= 0);
			assert(edgeP <= 1);

			vector<kMerInfo> localReturn;
			bool neverTrue = false;
			if(realdepth == 0)
			{
				if(firstEdgeGap == 1)
				{
					if(e->getEmission() != "_")
					{
						neverTrue = true;
						continue;
					}
				}
				else if (firstEdgeGap == -1)
				{
					if(e->getEmission() == "_")
					{
						neverTrue = true;
						continue;
					}
				}
			}
			assert(neverTrue == false);

			if(e->getEmission() != "_")
			{
				localReturn = forwardScanRec(targetNode, depth+1, realdepth + 1, limit, firstEdgeGap);
				for(vector<kMerInfo>::iterator kMerIt = localReturn.begin(); kMerIt != localReturn.end(); kMerIt++)
				{
					kMerIt->kMer_coded.push_back(e->emission);
					kMerIt->kMer_deCoded.push_back(e->getEmission());

					kMerIt->p = kMerIt->p*edgeP;
					kMerIt->traverseEdges.push_back(e);
				}
			}
			else
			{
				// this is a simple look-forward to check whether we can non-recursively extend
				if((targetNode->Outgoing_Edges.size()==1) && ((*targetNode->Outgoing_Edges.begin())->getEmission() == "_"))
				{
					vector<Edge*> traversedEdges;
					Node* currentTargetNode = targetNode;
					traversedEdges.push_back(e);
					while((currentTargetNode->Outgoing_Edges.size()==1) && ((*currentTargetNode->Outgoing_Edges.begin())->getEmission() == "_"))
					{
						traversedEdges.push_back((*currentTargetNode->Outgoing_Edges.begin()));
						currentTargetNode = (*currentTargetNode->Outgoing_Edges.begin())->To;
					}
					reverse(traversedEdges.begin(),traversedEdges.end());

					localReturn = forwardScanRec(currentTargetNode, depth, realdepth + 1, limit, firstEdgeGap);
					for(vector<kMerInfo>::iterator kMerIt = localReturn.begin(); kMerIt != localReturn.end(); kMerIt++)
					{
						kMerIt->p = kMerIt->p*edgeP;
						kMerIt->traverseEdges.insert(kMerIt->traverseEdges.end(), traversedEdges.begin(), traversedEdges.end());
					}
				}
				else
				{
					localReturn = forwardScanRec(targetNode, depth, realdepth + 1, limit, firstEdgeGap);
					for(vector<kMerInfo>::iterator kMerIt = localReturn.begin(); kMerIt != localReturn.end(); kMerIt++)
					{
						kMerIt->p = kMerIt->p*edgeP;
						kMerIt->traverseEdges.push_back(e);
					}
				}
			}

			returnMers.insert(returnMers.end(), localReturn.begin(), localReturn.end());

		}

		return returnMers;
	}
}


