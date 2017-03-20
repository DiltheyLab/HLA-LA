/*
 * processBAM.cpp
 *
 *  Created on: 21.09.2015
 *      Author: AlexanderDilthey
 */


#include <assert.h>
#include <string>
#include <exception>
#include <stdexcept>
#include <chrono>
#include <iostream>
#include "omp.h"
#include <tuple>

#include "processBAM.h"

#include "api/BamAux.h"
#include "utils/bamtools_utilities.h"
#include "../Utilities.h"

#include "reads/oneReadPair.h"

#include "../simulator/trueReadLevels.h"

namespace mapper {

processBAM::processBAM(std::string graphDir_, unsigned int maxThreads) {
	graphDir = graphDir_;

    // graph
    std::string graphFile = graphDir + "/PRG/graph.txt";
    assert(Utilities::fileExists(graphFile));
	
	std::string graph_serialized_fn = graphDir + "/serializedGRAPH";

	if(Utilities::fileExists(graph_serialized_fn) && (Utilities::fileLastWrite(graph_serialized_fn) > Utilities::fileLastWrite(graphFile)))
	{
		std::cout << Utilities::timestamp() << "Graph serialization existing and newer than graph file; read from " << graph_serialized_fn << "\n" << std::flush;
		
		std::ifstream serialization_istream;
		serialization_istream.open(graph_serialized_fn.c_str());
		assert(serialization_istream.is_open());
		boost::archive::text_iarchive archive_in(serialization_istream);
		archive_in >> g;
		serialization_istream.close();
		g->regenerateNodeIncomingOutgoingEdges();
		g->checkStructure();

		std::cout << Utilities::timestamp() << "\tdone.\n" << std::flush;		
	}
	else
	{
		g = new Graph();
		g->readFromFile(graphFile);
	}
	
	graphLevel_2_underlyingSequencePositions.clear();
	graphLevel_2_underlyingSequencePositions.resize(g->NodesPerLevel.size());

	for(unsigned int levelI = 0; levelI < g->NodesPerLevel.size(); levelI++)
	{
		g_level_names.push_back(g->getOneLocusIDforLevel(levelI));
	}
    // extended reference genome

	std::string extendedReferenceGenomePath;
	if(Utilities::fileExists(graphDir + "/extendedReferenceGenomePath.txt"))
	{
		extendedReferenceGenomePath  = Utilities::getFirstLine(graphDir + "/extendedReferenceGenomePath.txt");
	}
	else
	{
		extendedReferenceGenomePath  = graphDir + "/extendedReferenceGenome/extendedReferenceGenome.fa";
		assert(Utilities::fileExists(extendedReferenceGenomePath));
	}
    assert(Utilities::fileExists(extendedReferenceGenomePath));

    std::string PRGonlyReferenceGenomePath = graphDir + "/mapping_PRGonly/referenceGenome.fa";

    assert(Utilities::fileExists(PRGonlyReferenceGenomePath));

    extendedReferenceGenomeSequences = Utilities::readFASTA(extendedReferenceGenomePath, false);
    PRGonlyReferenceGenomeSequences = Utilities::readFASTA(PRGonlyReferenceGenomePath, false);
    assert(PRGonlyReferenceGenomeSequences.count("PRG_5") == 0);
    PRGonlyReferenceGenomeSequences["PRG_5"] = "N";

	assert(g->NodesPerLevel.size() > 1);
	inGraphGapStretch.resize(g->NodesPerLevel.size() - 1, false);
	size_t graphGapStretches = 0;
	int gapStretchMinimumLength = 3;
	auto addGapStretch = [&](int stretchStart_, int stretchStop_) -> void {
		int stretchL = stretchStop_ - stretchStart_ + 1;
		assert(stretchL >= 1);

		if(stretchL >= gapStretchMinimumLength)
		{
			for(int lII = stretchStart_; lII <= stretchStop_; lII++)
			{
				inGraphGapStretch.at(lII) = true;
			}

			graphGapStretches++;

			// std::cerr << "processBAM::processBAM(..) graph gap analysis: added long stretch from " << stretchStart_ << " to " << stretchStop_ << "\n" << std::flush;
		}
	};

	std::cerr << Utilities::timestamp() << "processBAM::processBAM(..): Start graph gap analysis.\n" << std::flush;
	
	int stretchStart = -1;
	for(unsigned int lI = 0; lI < (g->NodesPerLevel.size() - 1); lI++)
	{
		std::set<Edge*> edges = g->getEdgesEmanatingFromLevel(lI);
		bool haveGapEdge = false;
		for(std::set<Edge*>::iterator edgeIt = edges.begin(); edgeIt != edges.end(); edgeIt++)
		{
			Edge* e = *edgeIt;
			if(e->getEmission() == "_")
			{
				haveGapEdge++;
			}
		}

		if(haveGapEdge)
		{
			if(stretchStart == -1)
			{
				stretchStart = lI;
			}
		}
		else
		{
			if(stretchStart != -1)
			{
				int stretchStop = lI - 1;
				addGapStretch(stretchStart, stretchStop);
				stretchStart = -1;
			}
		}
	}

	if(stretchStart != -1)
	{
		int stretchStop = (g->NodesPerLevel.size() - 2);
		addGapStretch(stretchStart, stretchStop);
	}

	std::cerr << Utilities::timestamp() << "processBAM::processBAM(..) graph gap analysis: have " << graphGapStretches << " graph gap stretches; criterion length >= " << gapStretchMinimumLength << "\n" << std::flush;

	eA = new aligner::extensionAligner(g);

	threads_for_HLAtyping = maxThreads;
	
	paranoid = true;
}

processBAM::~processBAM() {
	delete(g);
	for(std::set<Edge*>::iterator eIt = edges_for_deletion.begin(); eIt != edges_for_deletion.end(); eIt++)
	{
		delete(*eIt);
	}
	delete(eA);
}

std::vector<std::string> processBAM::getReadIDs()
{
	std::set<std::string> forReturn;

	assert(_currentBAM.length());

	interestingIntervals_iterator = interestingIntervals.begin();
	interestingIntervals_i = 0;
	assert(R.Rewind());

	bool tryAgain = true;
	size_t examined_reads = 0;
	long long reads_included = 0;

	while(tryAgain)
	{
		if(interestingIntervals_iterator == interestingIntervals.end())
		{
			tryAgain = false;
			continue;
		}
		else if(interestingIntervals_i == (*interestingIntervals_iterator).second.size())
		{
			interestingIntervals_iterator++;
			interestingIntervals_i = 0;
			tryAgain = true;
			continue;
		}
		else
		{
			updateBAMSelector();
			std::string regionID =(*interestingIntervals_iterator).first;

			// std::cout << "interestingIntervals_i: " << interestingIntervals_i << " / " << (*interestingIntervals_iterator).second.size() <<  "\n" << std::flush;
			// std::cout << "\t" << "examined_reads: " << examined_reads << "\n" << std::flush;

			int leftBoundary_0based = (*interestingIntervals_iterator).second.at(interestingIntervals_i).start_1based - 1;
			int rightBoundary_0based = (*interestingIntervals_iterator).second.at(interestingIntervals_i).stop_1based - 1;

			while(R.GetNextAlignment(currentAlignment))
			{

				examined_reads++;

				if(currentAlignment.CigarData.size() == 0)
				{
					continue;
				}

				int alignmentStart = currentAlignment.Position;
				int alignmentStop = currentAlignment.GetEndPosition(false, true);

				assert(alignmentStart < alignmentStop); // might also be <=


				if(	((alignmentStart >= leftBoundary_0based) && (alignmentStart <= rightBoundary_0based)) &&
					((alignmentStop >= leftBoundary_0based) && (alignmentStop <= rightBoundary_0based))
				)
				{
					forReturn.insert(currentAlignment.Name);
				}
			}
			interestingIntervals_i++;
			tryAgain = true;
			continue;
		}
	}

	std::cout << "processBAM::extractSeeds(): getReadIDs " << examined_reads << " reads, collected " << forReturn.size() << " read IDs.\n" << std::flush;
	return std::vector<std::string>(forReturn.begin(), forReturn.end());
}

std::map<std::string, reads::protoSeeds> processBAM::extractSeeds_2BAMs(long long maximumIncludedReads)
{
	assert(_currentBAM1.length());
	assert(_currentBAM2.length());

	interestingIntervals_iterator = interestingIntervals.begin();
	interestingIntervals_i = 0;
	assert(R1.Rewind());
	assert(R2.Rewind());

	std::map<std::string, reads::protoSeeds> seeds;
	bool tryAgain = true;
	size_t examined_reads = 0;
	long long reads_included = 0;

	while(tryAgain)
	{
		if(interestingIntervals_iterator == interestingIntervals.end())
		{
			tryAgain = false;
			continue;
		}
		else if(interestingIntervals_i == (*interestingIntervals_iterator).second.size())
		{
			interestingIntervals_iterator++;
			interestingIntervals_i = 0;
			tryAgain = true;
			continue;
		}
		else
		{
			updateBAMSelectors();

			// std::cout << "interestingIntervals_i: " << interestingIntervals_i << " / " << (*interestingIntervals_iterator).second.size() <<  "\n" << std::flush;
			// std::cout << "\t" << "examined_reads: " << examined_reads << "\n" << std::flush;

			int leftBoundary_0based = (*interestingIntervals_iterator).second.at(interestingIntervals_i).start_1based - 1;
			int rightBoundary_0based = (*interestingIntervals_iterator).second.at(interestingIntervals_i).stop_1based - 1;
			std::string regionID =(*interestingIntervals_iterator).first;

			for(int readerI = 0; readerI <= 1; readerI++)
			{
				BamTools::BamReader& thisR = (readerI == 0) ? R1 : R2;

				while(thisR.GetNextAlignment(currentAlignment))
				{
					examined_reads++;

					if(currentAlignment.CigarData.size() == 0)
					{
						// std::cout << "Weird Read " << currentAlignment.Name << " - no CIGAR DATA\n";
						continue;
					}

					int alignmentStart = currentAlignment.Position;
					int alignmentStop = currentAlignment.GetEndPosition(false, true);
					if(!(alignmentStart < alignmentStop))
					{
						std::cerr << "alignmentStart: " << alignmentStart << "\n";
						std::cerr <<  "alignmentStart: " << alignmentStop << "\n" << std::flush;

						std::cerr << "\tPosition: " << currentAlignment.Position << "\n";
						std::cerr << "\talignedBases: " << currentAlignment.AlignedBases << "\n";
						std::cerr << "\tqueryBases: " << currentAlignment.QueryBases << "\n";
						std::cerr << "\tCIGAR: ";
						for(unsigned int cigarI = 0; cigarI < currentAlignment.CigarData.size(); cigarI++)
						{
							std::cerr << "cigarI: " << cigarI << "\n";
							std::cerr <<  currentAlignment.CigarData.at(cigarI).Type << currentAlignment.CigarData.at(cigarI).Length << " ";
						}
						std::cerr <<"\tEND CIGAR\n" << std::flush;
						std::cerr  << "\n" << std::flush;
						std::cerr << "\n" << std::flush;

					}
					assert(alignmentStart < alignmentStop); // might also be <=


					if(	((alignmentStart >= leftBoundary_0based) && (alignmentStart <= rightBoundary_0based)) &&
						((alignmentStop >= leftBoundary_0based) && (alignmentStop <= rightBoundary_0based))
					)
					{
						assert(! currentAlignment.IsPaired());
						assert(currentAlignment.BuildCharData());

						/*
						if(! currentAlignment.AlignedBases.length())
						{
							reads::protoSeeds::printAlignmentInfo(1, currentAlignment);
						}
						assert(currentAlignment.AlignedBases.length());
						*/

						std::string readID_no12 = currentAlignment.Name;
						int whichMate =  (readerI == 0) ? 1 : 2;

						seeds[readID_no12].takeAlignment(currentAlignment, whichMate, regionID, leftBoundary_0based, 0);

						_loadMapping((*interestingIntervals_iterator).second.at(interestingIntervals_i).PRGid);

						reads_included++;
						if((maximumIncludedReads > 0) && (reads_included >= maximumIncludedReads))
						{
							tryAgain = false;
							break;
						}
					}
					else
					{
						/*
						std::cout << "Read not completely contained " << currentAlignment.Name << "\n";
						std::cout << "\t" << alignmentStart << " - " << alignmentStop << "\n";
						std::cout << "\t" << leftBoundary_0based << " - " << rightBoundary_0based << "\n\n";
						*/
					}
				}
			}
			interestingIntervals_i++;
			tryAgain = true;
			continue;
		}
	}

	std::cout << "processBAM::extractSeeds_2BAms(): examined " << examined_reads << " reads, transformed into " << seeds.size() << " seeds.\n" << std::flush;
	return seeds;
}


void processBAM::extractSeedsInto(std::map<std::string, reads::protoSeeds>& seeds, int whichReader)
{
	assert(_currentBAM.length());

	interestingIntervals_iterator = interestingIntervals.begin();
	interestingIntervals_i = 0;
	assert(R.Rewind());

	bool tryAgain = true;
	size_t examined_reads = 0;
	long long reads_included = 0;

	while(tryAgain)
	{
		if(interestingIntervals_iterator == interestingIntervals.end())
		{
			tryAgain = false;
			continue;
		}
		else if(interestingIntervals_i == (*interestingIntervals_iterator).second.size())
		{
			interestingIntervals_iterator++;
			interestingIntervals_i = 0;
			tryAgain = true;
			continue;
		}
		else
		{
			updateBAMSelector();
			std::string regionID =(*interestingIntervals_iterator).first;

			// std::cout << "interestingIntervals_i: " << interestingIntervals_i << " / " << (*interestingIntervals_iterator).second.size() <<  "\n" << std::flush;
			// std::cout << "\t" << "examined_reads: " << examined_reads << "\n" << std::flush;

			int leftBoundary_0based = (*interestingIntervals_iterator).second.at(interestingIntervals_i).start_1based - 1;
			int rightBoundary_0based = (*interestingIntervals_iterator).second.at(interestingIntervals_i).stop_1based - 1;

			while(R.GetNextAlignment(currentAlignment))
			{
				examined_reads++;

				if(currentAlignment.CigarData.size() == 0)
				{
					// std::cout << "Weird Read " << currentAlignment.Name << " - no CIGAR DATA\n";
					continue;
				}

				int alignmentStart = currentAlignment.Position;
				int alignmentStop = currentAlignment.GetEndPosition(false, true);

				/*
				if((currentAlignment.Name == "H781JADXX131217:1:2115:10484:82360") || (!(alignmentStart < alignmentStop)))
				{
					std::cerr << "ReadID: " << currentAlignment.Name << "\n";
					std::cerr << "Reference: " << R.GetReferenceData().at(currentAlignment.RefID).RefName << "\n";
					std::cerr << "alignmentStart: " << alignmentStart << "\n";
					std::cerr <<  "alignmentStart: " << alignmentStop << "\n" << std::flush;

					std::cerr << "\tPosition: " << currentAlignment.Position << "\n";
					std::cerr << "\talignedBases: " << currentAlignment.AlignedBases << "\n";
					std::cerr << "\tqueryBases: " << currentAlignment.QueryBases << "\n";
					std::cerr << "\tCIGAR: ";
					// for(unsigned int cigarI = 0; cigarI < currentAlignment.CigarData.size(); cigarI++)
					// {
						// std::cerr << "cigarI: " << cigarI << "\n";
						// std::cerr <<  currentAlignment.CigarData.at(cigarI).Type << currentAlignment.CigarData.at(cigarI).Length << " ";
					// }
					std::cerr <<"\tEND CIGAR\n" << std::flush;
					std::cerr  << "\n" << std::flush;
					std::cerr << "\t" <<
					std::cerr << "\n" << std::flush;

				}
				*/
				assert(alignmentStart < alignmentStop); // might also be <=


				if(	((alignmentStart >= leftBoundary_0based) && (alignmentStart <= rightBoundary_0based)) &&
					((alignmentStop >= leftBoundary_0based) && (alignmentStop <= rightBoundary_0based))
				)
				{
					if((currentAlignment.Name == "H781JADXX131217:1:2115:10484:82360"))
					{
						std::cerr << "Take this seed!\n";
					}

					assert(currentAlignment.IsPaired());
					assert(currentAlignment.BuildCharData());

					std::string CIGARstring;
					unsigned int alignment_clipping_deletions = 0;
					unsigned int alignment_length = 0;
					for(auto cigarElement : currentAlignment.CigarData)
					{
						if((cigarElement.Type == 'H') || (cigarElement.Type == 'S') || (cigarElement.Type == 'I') || (cigarElement.Type == 'D'))
						{
							alignment_clipping_deletions += cigarElement.Length;
						}
						alignment_length += cigarElement.Length;
						CIGARstring += Utilities::ItoStr(cigarElement.Length);
						CIGARstring.push_back(cigarElement.Type );
					}

					double prop_clipping_deletions = (double)alignment_clipping_deletions / (double) alignment_length;
					assert(prop_clipping_deletions >= 0);
					assert(prop_clipping_deletions <= 1);

					//if((! currentAlignment.IsPrimaryAlignment()) || (prop_clipping_deletions <= 0.2))
					{
						/*
						if(! currentAlignment.AlignedBases.length())
						{
							reads::protoSeeds::printAlignmentInfo(1, currentAlignment);
						}
						assert(currentAlignment.AlignedBases.length());
						*/

						std::string readID_no12 = currentAlignment.Name;
						int whichMate =  (currentAlignment.IsFirstMate()) ? 1 : 2;

						seeds[readID_no12].takeAlignment(currentAlignment, whichMate, regionID, leftBoundary_0based, whichReader);

						_loadMapping((*interestingIntervals_iterator).second.at(interestingIntervals_i).PRGid);

						reads_included++;
					}
					/*
					else
					{
						std::cout << "Read - " << CIGARstring << " " << prop_clipping_deletions << "\n";
					}
					*/
				}
				else
				{
					/*
					std::cout << "Read not completely contained " << currentAlignment.Name << "\n";
					std::cout << "\t" << alignmentStart << " - " << alignmentStop << "\n";
					std::cout << "\t" << leftBoundary_0based << " - " << rightBoundary_0based << "\n\n";
					*/
				}
			}
			interestingIntervals_i++;
			tryAgain = true;
			continue;
		}
	}

	std::cout << "processBAM::extractSeeds(): examined " << examined_reads << " reads, transformed into " << seeds.size() << " seeds.\n" << std::flush;
}


std::map<std::string, reads::protoSeeds> processBAM::extractSeeds(long long maximumIncludedReads, std::set<std::string> limitToReadIDs)
{
	assert(_currentBAM.length());

	interestingIntervals_iterator = interestingIntervals.begin();
	interestingIntervals_i = 0;
	assert(R.Rewind());

	std::map<std::string, reads::protoSeeds> seeds;
	bool tryAgain = true;
	size_t examined_reads = 0;
	long long reads_included = 0;
	long long reads_included_complete = 0;

	while(tryAgain)
	{
		if(interestingIntervals_iterator == interestingIntervals.end())
		{
			tryAgain = false;
			continue;
		}
		else if(interestingIntervals_i == (*interestingIntervals_iterator).second.size())
		{
			interestingIntervals_iterator++;
			interestingIntervals_i = 0;
			tryAgain = true;
			continue;
		}
		else
		{
			updateBAMSelector();
			std::string regionID =(*interestingIntervals_iterator).first;

			// std::cout << "interestingIntervals_i: " << interestingIntervals_i << " / " << (*interestingIntervals_iterator).second.size() <<  "\n" << std::flush;
			// std::cout << "\t" << "examined_reads: " << examined_reads << "\n" << std::flush;

			int leftBoundary_0based = (*interestingIntervals_iterator).second.at(interestingIntervals_i).start_1based - 1;
			int rightBoundary_0based = (*interestingIntervals_iterator).second.at(interestingIntervals_i).stop_1based - 1;

			while(R.GetNextAlignment(currentAlignment))
			{
				if(limitToReadIDs.size() && (limitToReadIDs.count(currentAlignment.Name) == 0))
				{
					continue;
				}

				examined_reads++;
				
				if(currentAlignment.CigarData.size() == 0)
				{
					// std::cout << "Weird Read " << currentAlignment.Name << " - no CIGAR DATA\n";
					continue;
				}

				int alignmentStart = currentAlignment.Position;
				int alignmentStop = currentAlignment.GetEndPosition(false, true);
				
				if((currentAlignment.Name == "H781JADXX131217:1:2115:10484:82360") || (!(alignmentStart < alignmentStop)))
				{
					std::cerr << "ReadID: " << currentAlignment.Name << "\n";
					std::cerr << "Reference: " << R.GetReferenceData().at(currentAlignment.RefID).RefName << "\n";
					std::cerr << "alignmentStart: " << alignmentStart << "\n";
					std::cerr <<  "alignmentStart: " << alignmentStop << "\n" << std::flush;

					std::cerr << "\tPosition: " << currentAlignment.Position << "\n";
					std::cerr << "\talignedBases: " << currentAlignment.AlignedBases << "\n";
					std::cerr << "\tqueryBases: " << currentAlignment.QueryBases << "\n";
					std::cerr << "\tCIGAR: ";
					// for(unsigned int cigarI = 0; cigarI < currentAlignment.CigarData.size(); cigarI++)
					// {
						// std::cerr << "cigarI: " << cigarI << "\n";
						// std::cerr <<  currentAlignment.CigarData.at(cigarI).Type << currentAlignment.CigarData.at(cigarI).Length << " ";
					// }
					std::cerr <<"\tEND CIGAR\n" << std::flush;
					std::cerr  << "\n" << std::flush;
					std::cerr << "\n" << std::flush;

				}
				assert(alignmentStart < alignmentStop); // might also be <=


				if(	((alignmentStart >= leftBoundary_0based) && (alignmentStart <= rightBoundary_0based)) &&
					((alignmentStop >= leftBoundary_0based) && (alignmentStop <= rightBoundary_0based))
				)
				{
					if((currentAlignment.Name == "H781JADXX131217:1:2115:10484:82360"))
					{
						std::cerr << "Take this seed!\n";
					}
					
					assert(currentAlignment.IsPaired());
					assert(currentAlignment.BuildCharData());

					std::string CIGARstring;
					unsigned int alignment_clipping_deletions = 0;
					unsigned int alignment_length = 0;
					for(auto cigarElement : currentAlignment.CigarData)
					{
						if((cigarElement.Type == 'H') || (cigarElement.Type == 'S') || (cigarElement.Type == 'I') || (cigarElement.Type == 'D'))
						{
							alignment_clipping_deletions += cigarElement.Length;
						}
						alignment_length += cigarElement.Length;
						CIGARstring += Utilities::ItoStr(cigarElement.Length);
						CIGARstring.push_back(cigarElement.Type );
					}

					double prop_clipping_deletions = (double)alignment_clipping_deletions / (double) alignment_length;
					assert(prop_clipping_deletions >= 0);
					assert(prop_clipping_deletions <= 1);
					
					//if((! currentAlignment.IsPrimaryAlignment()) || (prop_clipping_deletions <= 0.2))
					{
						/*
						if(! currentAlignment.AlignedBases.length())
						{
							reads::protoSeeds::printAlignmentInfo(1, currentAlignment);
						}
						assert(currentAlignment.AlignedBases.length());
						*/

						std::string readID_no12 = currentAlignment.Name;
						int whichMate =  (currentAlignment.IsFirstMate()) ? 1 : 2;

						bool thisIsANewRead = (seeds.count(readID_no12) == 0);
						seeds[readID_no12].takeAlignment(currentAlignment, whichMate, regionID, leftBoundary_0based, 0);

						// int score = getAlignmentScore(std::get<2>(al));
						// 		read1_alignments.push_back(make_tuple(referenceID, reference2level_offset_0based, a, whichReader));
						// 		std::pair<int, int> al_startStop_PRG = alignment_get_startstop_PRGcoordinates(al);

						_loadMapping((*interestingIntervals_iterator).second.at(interestingIntervals_i).PRGid);

						if(thisIsANewRead)
						{
							reads_included++;
						}
						
						if(seeds.at(readID_no12).isComplete())
						{
							reads_included_complete++;
						}
						
						if((maximumIncludedReads > 0) && (reads_included >= maximumIncludedReads) && (reads_included_complete >= (maximumIncludedReads/2)))
						{
							tryAgain = false;
							break;
						}
					} 
					/*
					else
					{						
						std::cout << "Read - " << CIGARstring << " " << prop_clipping_deletions << "\n";
					}
					*/
				}
				else
				{
					/*
					std::cout << "Read not completely contained " << currentAlignment.Name << "\n";
					std::cout << "\t" << alignmentStart << " - " << alignmentStop << "\n";
					std::cout << "\t" << leftBoundary_0based << " - " << rightBoundary_0based << "\n\n";
					*/
				}
			}
			interestingIntervals_i++;
			tryAgain = true;
			continue;
		}
	}

	std::cout << "processBAM::extractSeeds(): examined " << examined_reads << " reads, transformed into " << seeds.size() << " seeds.\n" << std::flush;
	return seeds;
}

std::map<std::string, reads::protoSeeds> processBAM::extractSeeds2(std::set<std::string> limitToReadIDs)
{
	std::cout << Utilities::timestamp() << "\t\t\tStart extractSeeds2\n" << std::flush;

	assert(_currentBAM.length());

	interestingIntervals_iterator = interestingIntervals.begin();
	interestingIntervals_i = 0;
	assert(R.Rewind());

	std::map<std::string, reads::protoSeeds> seeds;
	bool tryAgain = true;
	size_t examined_reads = 0;
	long long reads_included = 0;

	const std::vector<BamTools::RefData> bamInternalID2String = R.GetReferenceData();

	while(R.GetNextAlignmentCore(currentAlignment))
	{
		if(! currentAlignment.IsMapped() )
		{
			continue;
		}

		size_t alignment_refID = currentAlignment.RefID;
		const std::string& alignment_refID_string = bamInternalID2String.at(alignment_refID).RefName;

		if(interestingIntervals.count(alignment_refID_string))
		{
			for(auto interestingInterval : interestingIntervals.at(alignment_refID_string))
			{
				int leftBoundary_0based = interestingInterval.start_1based - 1;
				int rightBoundary_0based = interestingInterval.stop_1based - 1;

				assert(currentAlignment.BuildCharData());

				if(limitToReadIDs.size() && (limitToReadIDs.count(currentAlignment.Name) == 0))
				{
					continue;
				}

				examined_reads++;

				if(currentAlignment.CigarData.size() == 0)
				{
					// std::cout << "Weird Read " << currentAlignment.Name << " - no CIGAR DATA\n";
					continue;
				}


				int alignment_start_0based = currentAlignment.Position;
				int alignment_stop_0based = currentAlignment.GetEndPosition(false, true);

				if(	((alignment_start_0based >= leftBoundary_0based) && (alignment_start_0based <= rightBoundary_0based)) &&
					((alignment_stop_0based >= leftBoundary_0based) && (alignment_stop_0based <= rightBoundary_0based))
				)
				{
					if((currentAlignment.Name == "H781JADXX131217:1:2115:10484:82360"))
					{
						std::cerr << "Take this seed!\n";
					}

					assert(currentAlignment.IsPaired());

					std::string CIGARstring;
					unsigned int alignment_clipping_deletions = 0;
					unsigned int alignment_length = 0;
					for(auto cigarElement : currentAlignment.CigarData)
					{
						if((cigarElement.Type == 'H') || (cigarElement.Type == 'S') || (cigarElement.Type == 'I') || (cigarElement.Type == 'D'))
						{
							alignment_clipping_deletions += cigarElement.Length;
						}
						alignment_length += cigarElement.Length;
						CIGARstring += Utilities::ItoStr(cigarElement.Length);
						CIGARstring.push_back(cigarElement.Type );
					}

					double prop_clipping_deletions = (double)alignment_clipping_deletions / (double) alignment_length;
					assert(prop_clipping_deletions >= 0);
					assert(prop_clipping_deletions <= 1);

					//if((! currentAlignment.IsPrimaryAlignment()) || (prop_clipping_deletions <= 0.2))
					{
						/*
						if(! currentAlignment.AlignedBases.length())
						{
							reads::protoSeeds::printAlignmentInfo(1, currentAlignment);
						}
						assert(currentAlignment.AlignedBases.length());
						*/

						std::string readID_no12 = currentAlignment.Name;
						int whichMate =  (currentAlignment.IsFirstMate()) ? 1 : 2;

						seeds[readID_no12].takeAlignment(currentAlignment, whichMate, alignment_refID_string, leftBoundary_0based, 0);

						// int score = getAlignmentScore(std::get<2>(al));
						// 		read1_alignments.push_back(make_tuple(referenceID, reference2level_offset_0based, a, whichReader));
						// 		std::pair<int, int> al_startStop_PRG = alignment_get_startstop_PRGcoordinates(al);

						_loadMapping((*interestingIntervals_iterator).second.at(interestingIntervals_i).PRGid);

						reads_included++;
						//if((maximumIncludedReads > 0) && (reads_included >= maximumIncludedReads))
						//{
						//	tryAgain = false;
						//	break;
						//}
					}
					/*
					else
					{
						std::cout << "Read - " << CIGARstring << " " << prop_clipping_deletions << "\n";
					}
					*/
				}

			}
		}

	}

	std::cout << Utilities::timestamp() << "\t\t\tDone extractSeeds2\n" << std::flush;
	
	std::cout << "processBAM::extractSeeds2(): examined " << examined_reads << " reads, transformed into " << seeds.size() << " seeds.\n" << std::flush;
	return seeds;
}

std::pair<double, double> processBAM::estimateInsertSize_noGraph(std::string BAM)
{
	std::set<std::string> usedReadIDs;

	std::map<int, double> IS_combined_counts;

	BamTools::BamReader localReader;
	assert(localReader.Open(BAM));
	localReader.LocateIndex();
	if ( ! localReader.HasIndex() )
	{
		throw std::runtime_error("File "+BAM+" does not seem to be indexed - please specify indexed BAM!");
	}
		
	BamTools::BamAlignment localAlignment;
	std::map<int, int> referenceIDs_byLength;
	for(int refID = 0; refID < (int)localReader.GetReferenceData().size(); refID++)
	{
		int L = localReader.GetReferenceData().at(refID).RefLength;
		referenceIDs_byLength[refID] = L;
	}

	size_t attemptedAccesses = 0;
	size_t wantReads = 5000;
	bool lastAttemptNoAlignments = false;
	while(usedReadIDs.size() < wantReads)
	{
		if(((usedReadIDs.size()  % 100) == 0) || lastAttemptNoAlignments)
		{
			int nextRefID = Utilities::choose_from_nonnormalized_map(referenceIDs_byLength);
			int nextRefID_length = referenceIDs_byLength.at(nextRefID);
			int nextPosition = Utilities::randomNumber(nextRefID_length-1);
			
			bool jumpSuccess = localReader.Jump(nextRefID, nextPosition);
			if(! jumpSuccess)			
			{
				std::cerr << "Could not jump to " << nextRefID << " - " << nextPosition << " (length: )" << referenceIDs_byLength.at(nextRefID) << "\n" << std::flush;
			}
			// std::cerr << "Jumped to " << nextRefID << " - " << nextPosition << " (length: )" << referenceIDs_byLength.at(nextRefID) << "\n" << std::flush;			
			assert(jumpSuccess);
		}
		
		// std::cerr << "usedReadIDs.size(): " << usedReadIDs.size() << "   --   attemptedAccesses: " << attemptedAccesses << "\n" << std::flush; 

		attemptedAccesses++;
		if(attemptedAccesses > (1000 * wantReads))
		{
			std::cerr << "Leave because " << attemptedAccesses << " > " << (5 * usedReadIDs.size()) << "\n" << std::flush;
			break;
		}

		if(! localReader.GetNextAlignment(localAlignment))
		{
			lastAttemptNoAlignments = true;
			continue;
		}
		lastAttemptNoAlignments = false;
					
		if(usedReadIDs.count(localAlignment.Name))
		{
			continue;
		}

		if(! localAlignment.IsPaired())
		{
			continue;
		}
		
		if(! localAlignment.IsMapped())
		{
			continue;
		}	

		if(localAlignment.IsFailedQC())
		{
			continue;
		}
		
		if(! localAlignment.IsPrimaryAlignment())
		{
			continue;
		}
		
		if(! localAlignment.IsProperPair())
		{
			continue;
		}
		
		if(! localAlignment.IsFirstMate())
		{
			continue;
		}
		
		if(localAlignment.IsReverseStrand())
		{
			continue;
		}
		
		usedReadIDs.insert(localAlignment.Name);

		assert(localAlignment.RefID == localAlignment.MateRefID);

		int lastPosition_thisAlignment = localAlignment.GetEndPosition();
		int distance_begin_nextAlignment = localAlignment.MatePosition - lastPosition_thisAlignment;

		if(IS_combined_counts.count(distance_begin_nextAlignment) == 0)
		{
			IS_combined_counts[distance_begin_nextAlignment] = 0;
		}

		IS_combined_counts[distance_begin_nextAlignment]++;
	}

	if(usedReadIDs.size() != wantReads)
	{
		std::cerr << "estimateInsertSize_noGraph(..): Wanted to collect " << wantReads << ", but have only " << usedReadIDs.size() << "\n";
		std::cerr << "BAM: " << BAM << "\n\n" << std::flush;
	}

	assert(usedReadIDs.size() == wantReads);


	return calculateInsertSizeFromHistogram(IS_combined_counts, false);
}

std::pair<double, double> processBAM::calculateInsertSizeFromHistogram(const std::map<int, double>& IS_combined_counts, bool verbose)
{

	std::set<int> IS_keys;
	double IS_total_size = 0;
	for(auto ISentry : IS_combined_counts)
	{
		IS_keys.insert(ISentry.first);
		assert(ISentry.second >= 0);
		IS_total_size += ISentry.second;
	}

	double cumulative_sum = 0;
	double weighted_median = 0;
	double weighted_20 = 0;
	double weighted_80 = 0;
	bool set_median = false;
	bool set_weighted_20 = false;
	bool set_weighted_80 = false;

	if(verbose)
		std::cout << "\n\nIS histogram over " << IS_total_size << " read pairs:\n";

	for(std::set<int>::iterator ISit = IS_keys.begin(); ISit != IS_keys.end(); ISit++)
	{
		int d = *ISit;
		double c = IS_combined_counts.at(d);

		if(verbose)
			std::cout << "\t" << d << ": " << c << "\n";

		cumulative_sum += IS_combined_counts.at(d);
		if((set_median == false) && (cumulative_sum >= (IS_total_size * 0.5)))
		{
			weighted_median = d;
			set_median = true;
		}
		if((set_weighted_20 == false) && (cumulative_sum >= (IS_total_size * 0.2)))
		{
			weighted_20 = d;
			set_weighted_20 = true;
		}
		if((set_weighted_80 == false) && (cumulative_sum >= (IS_total_size * 0.8)))
		{
			weighted_80 = d;
			set_weighted_80 = true;
		}
	}

	if(verbose)
		std::cout << "\n" << std::flush;

	assert(set_weighted_80 && set_weighted_20 && set_median);

	if(verbose)
	{
		std::cout << "Summary statistics:\n";
		std::cout << "\t" << "Median: " << weighted_median << "\n";
		std::cout << "\t" << "20%: " << weighted_20 << "\n";
		std::cout << "\t" << "80%: " << weighted_80 << "\n";
		std::cout << "\n" << std::flush;
	}

	double f_mean_ret = weighted_median;
	weighted_20 = abs(weighted_median - weighted_20);
	weighted_80 = abs(weighted_median - weighted_80);
	double f_sd_ret = (weighted_20 > weighted_80) ? weighted_20 : weighted_80;

	if(verbose)
	{
		std::cout << "\n\nAssumed values:\n";
		std::cout << "\t" << "Mean: " << f_mean_ret << "\n";
		std::cout << "\t" << "SD: " << f_sd_ret << "\n";
		std::cout << "\n" << std::flush;
		std::cout << "processBAM::calculateInsertSizeFromHistogram(): Utilized " << IS_total_size << " entries.\n";
	}

	return make_pair(f_mean_ret, f_sd_ret);
}

std::pair<double, double> processBAM::estimateInsertSize(std::string BAM, bool extendedReferenceGenome)
{
	initBAM(BAM, extendedReferenceGenome, false);

	std::map<std::string, reads::protoSeeds> seeds = extractSeeds(4000);

	std::set<std::string> incompleteSeeds;
	std::set<std::string> completeProtoSeeds;

	for(std::map<std::string, reads::protoSeeds>::iterator seedIt = seeds.begin(); seedIt != seeds.end(); seedIt++)
	{
		if(! seedIt->second.isComplete())
		{
			incompleteSeeds.insert(seedIt->first);
			// seedIt->second.printDebug(R);
		}
		else
		{
			completeProtoSeeds.insert(seedIt->first);
		}
	}
	std::cout << "processBAM::estimateInsertSize(): Deal " << seeds.size() << " total proto-seeds (i.e. read pairs), of which " << incompleteSeeds.size() << " are incomplete.\n" << std::flush;


	std::map<int, double> IS_combined_counts;

	size_t used_proto_seeds = 0;
	size_t skipped_proto_seeds = 0;
	for(std::set<std::string>::iterator seedIt = completeProtoSeeds.begin(); seedIt != completeProtoSeeds.end(); seedIt++)
	{
		const reads::protoSeeds& protoSeed = seeds.at(*seedIt);

		bool examinedOneCombination = false;

		for(unsigned int alignmentIndex1 = 0; alignmentIndex1 < protoSeed.read1_alignments.size(); alignmentIndex1++)
		{
			if(std::get<2>(protoSeed.read1_alignments.at(alignmentIndex1)).IsPrimaryAlignment() == false)
			{
				continue;
			}

			for(unsigned int alignmentIndex2 = 0; alignmentIndex2 < protoSeed.read2_alignments.size(); alignmentIndex2++)
			{
				if(std::get<2>(protoSeed.read2_alignments.at(alignmentIndex2)).IsPrimaryAlignment() == false)
				{
					continue;
				}

				examinedOneCombination = true;

				used_proto_seeds++;

				const std::tuple<std::string, int,  BamTools::BamAlignment, int>& al_r1 = protoSeed.read1_alignments.at(alignmentIndex1);
				const std::tuple<std::string, int,  BamTools::BamAlignment, int>& al_r2 = protoSeed.read2_alignments.at(alignmentIndex2);

				reads::verboseSeedChain al_r1_chain = alignment2Chain(al_r1, std::get<2>(protoSeed.read1_alignments.at(alignmentIndex1)).QueryBases, std::get<2>(protoSeed.read1_alignments.at(alignmentIndex1)).Qualities);
				reads::verboseSeedChain al_r2_chain = alignment2Chain(al_r2, std::get<2>(protoSeed.read2_alignments.at(alignmentIndex2)).QueryBases, std::get<2>(protoSeed.read2_alignments.at(alignmentIndex2)).Qualities);

				if(paranoid)
				{
					al_r1_chain.checkChainConcordanceWithSequence(std::get<2>(protoSeed.read1_alignments.at(alignmentIndex1)).QueryBases);
					al_r2_chain.checkChainConcordanceWithSequence(std::get<2>(protoSeed.read2_alignments.at(alignmentIndex2)).QueryBases);
				}

				reads::verboseSeedChain al_r1_chain_extended = eA->extendSeedChain(std::get<2>(protoSeed.read1_alignments.at(alignmentIndex1)).QueryBases, al_r1_chain);
				reads::verboseSeedChain al_r2_chain_extended = eA->extendSeedChain(std::get<2>(protoSeed.read2_alignments.at(alignmentIndex2)).QueryBases, al_r2_chain);

				bool strandsValid = eA->alignedReadPair_strandsValid(al_r1_chain_extended, al_r2_chain_extended);

				if(strandsValid)
				{
					std::set<int> underlyingSequencesDistances = eA->alignedReadPair_pairsDistancesUnderlyingSequences(al_r1_chain_extended, al_r2_chain_extended, graphLevel_2_underlyingSequencePositions);

					for(std::set<int>::iterator ISiterator = underlyingSequencesDistances.begin(); ISiterator != underlyingSequencesDistances.end(); ISiterator++)
					{
						int IS = *ISiterator;
						if(IS_combined_counts.count(IS) == 0)
						{
							IS_combined_counts[IS] = 0;
						}

						double weight = 1.0/(double)underlyingSequencesDistances.size();
						IS_combined_counts[IS] += weight;
					}
				}
				else
				{
					skipped_proto_seeds++;
				}

				if(examinedOneCombination)
				{
					break;
				}
			}

			if(examinedOneCombination)
			{
				break;
			}
		}

	}

	std::pair<double, double> forReturn = calculateInsertSizeFromHistogram(IS_combined_counts, false);

	std::cout << "processBAM::estimateInsertSize(): Utilized " << used_proto_seeds << " proto-seeds, skipped " << skipped_proto_seeds << "\n" << std::flush;
	std::cout << "\t IS estimate " << forReturn.first << " (mean) / " << forReturn.second << " (sd).\n" << std::flush;

	return forReturn;
}

void processBAM::initBAM(std::string BAM, bool extendedReferenceGenome, bool multiSampleMapping)
{
	if((_currentBAM != BAM) || (_currentBAM_isExtendedReferenceGenome != extendedReferenceGenome))
	{
		assert(Utilities::fileExists(BAM));
		R.Open(BAM);

		R.LocateIndex();
		if ( ! R.HasIndex() )
		{
			throw std::runtime_error("File "+BAM+" does not seem to be indexed - please specify indexed BAM!");
		}

		std::string extendedReferenceGenomePath;
		if(Utilities::fileExists(graphDir + "/extendedReferenceGenomePath.txt"))
		{
			extendedReferenceGenomePath  = Utilities::getFirstLine(graphDir + "/extendedReferenceGenomePath.txt");
		}
		else
		{
			extendedReferenceGenomePath  = graphDir + "/extendedReferenceGenome/extendedReferenceGenome.fa";
			assert(Utilities::fileExists(extendedReferenceGenomePath));
		}
	    assert(Utilities::fileExists(extendedReferenceGenomePath));
	    std::string PRGonlyReferenceGenomePath = graphDir + "/mapping_PRGonly/referenceGenome.fa";

		// sequences
		std::string sequencesFile = graphDir + "/sequences.txt";

		std::ifstream sequencesStream;
		sequencesStream.open(sequencesFile.c_str());
		assert(sequencesStream.is_open());

		std::string headerLine;
		std::string line;
		assert(sequencesStream.good());
		std::getline(sequencesStream, headerLine);
		Utilities::eraseNL(headerLine);
		std::vector<std::string> headerFields = Utilities::split(headerLine, "\t");

		size_t interestingIntervals_n = 0;

		PRGid_2_BAMid.clear();
		BAMid_2_PRGid.clear();
		interestingIntervals.clear();


		while(sequencesStream.good())
		{
			std::getline(sequencesStream, line);
			Utilities::eraseNL(line);
			if(line.length() == 0)
				continue;

			std::vector<std::string> lineFields = Utilities::split(line, "\t");
			assert(lineFields.size() == headerFields.size());

			std::map<std::string, std::string> L;
			for(unsigned int i = 0; i < headerFields.size(); i++)
			{
				L[headerFields.at(i)] = lineFields.at(i);
			}

			std::string IDforBAMIndex;
			int PRGid = Utilities::StrtoI(L.at("SequenceID"));

			bool inThisBAM = true;
			int startIndex_1based;
			int stopIndex_1based;
			if(L.at("Chr") != "")
			{
				IDforBAMIndex = L.at("Chr");
				int internalBAMid = R.GetReferenceID(IDforBAMIndex);

				if(extendedReferenceGenomeSequences.count(IDforBAMIndex) == 0)
				{
					throw std::runtime_error(IDforBAMIndex + " cannot be found in the extended reference genome "+extendedReferenceGenomePath);
				}
				if(PRGonlyReferenceGenomeSequences.count(IDforBAMIndex) == 0)
				{
					throw std::runtime_error(IDforBAMIndex + " cannot be found in the PRG-only reference genome "+PRGonlyReferenceGenomePath);
				}

				if(R.GetReferenceID(IDforBAMIndex) == -1)
				{
					inThisBAM = false;
					if(! multiSampleMapping)
					{
						if(! Utilities::sequence_all_Ns(extendedReferenceGenomeSequences.at(IDforBAMIndex)))
						{
							throw std::runtime_error(IDforBAMIndex + " cannot be found in the index of "+BAM);
						}
						else
						{
							continue;
						}
					}
				}


				if(inThisBAM)
				{
					BamTools::RefData internalBAMSpecification = R.GetReferenceData().at(internalBAMid);
					assert(internalBAMSpecification.RefName == IDforBAMIndex);

					if(extendedReferenceGenome)
					{
						startIndex_1based = Utilities::StrtoI(L.at("Start_1based"));
						stopIndex_1based = Utilities::StrtoI(L.at("Stop_1based"));
						assert(startIndex_1based <= internalBAMSpecification.RefLength);
						assert(stopIndex_1based <= internalBAMSpecification.RefLength);
					}
					else
					{
						startIndex_1based = 1;
						stopIndex_1based = PRGonlyReferenceGenomeSequences.at(IDforBAMIndex).length();
						assert(stopIndex_1based == internalBAMSpecification.RefLength);
					}
				}
				else
				{
					if(extendedReferenceGenome)
					{
						startIndex_1based = Utilities::StrtoI(L.at("Start_1based"));
						stopIndex_1based = Utilities::StrtoI(L.at("Stop_1based"));
					}
					else
					{
						startIndex_1based = 1;
						stopIndex_1based = PRGonlyReferenceGenomeSequences.at(IDforBAMIndex).length();
					}
				}
			}
			else
			{
				IDforBAMIndex = "PRG_" + L.at("SequenceID");
				int internalBAMid = R.GetReferenceID(IDforBAMIndex);

				if(extendedReferenceGenomeSequences.count(IDforBAMIndex) == 0)
				{
					throw std::runtime_error(IDforBAMIndex + " cannot be found in the extended reference genome "+extendedReferenceGenomePath);
				}
				if(PRGonlyReferenceGenomeSequences.count(IDforBAMIndex) == 0)
				{
					throw std::runtime_error(IDforBAMIndex + " cannot be found in the PRG-only reference genome "+PRGonlyReferenceGenomePath);
				}


				if(R.GetReferenceID(IDforBAMIndex) == -1)
				{
					inThisBAM = false;

					if(! multiSampleMapping)
					{
						if(! Utilities::sequence_all_Ns(extendedReferenceGenomeSequences.at(IDforBAMIndex)))
						{
							throw std::runtime_error(IDforBAMIndex + " cannot be found in the index of "+BAM);
						}
						else
						{
							continue;
						}
					}
				}


				if(inThisBAM)
				{
					BamTools::RefData internalBAMSpecification = R.GetReferenceData().at(internalBAMid);
					assert(internalBAMSpecification.RefName == IDforBAMIndex);

					startIndex_1based = 1;
					stopIndex_1based = extendedReferenceGenomeSequences.at(IDforBAMIndex).length();
					assert(stopIndex_1based > startIndex_1based);
					if(!(stopIndex_1based == internalBAMSpecification.RefLength))
					{
						std::cerr << "Problem in sequences " << sequencesFile << "! for " << "PRG_" + L.at("SequenceID") << "\n";
						std::cerr << "\t" << "BAM" << ": " <<  BAM << "\n";
						std::cerr << "\t" << "stopIndex_1based" << ": " << stopIndex_1based << "\n";
						std::cerr << "\t" << "startIndex_1based" << ": " << startIndex_1based << "\n";
						std::cerr << "\t" << "internalBAMSpecification.RefLength" << ": " << internalBAMSpecification.RefLength << "\n";
						std::cerr << "\t" << "IDforBAMIndex" << ": " << IDforBAMIndex << "\n";
						std::cerr << "\t" << "extendedReferenceGenomePath" << ": " << Utilities::getFirstLine(graphDir + "/extendedReferenceGenomePath.txt") << "\n";
						std::cerr << std::flush;
					}
					assert(stopIndex_1based == internalBAMSpecification.RefLength);
				}
				else
				{
					startIndex_1based = 1;
					stopIndex_1based = extendedReferenceGenomeSequences.at(IDforBAMIndex).length();
					assert(stopIndex_1based > startIndex_1based);
				}

			}

			if(inThisBAM)
			{
				assert(R.GetReferenceID(IDforBAMIndex) != -1);
			}
			assert(extendedReferenceGenomeSequences.count(IDforBAMIndex));

			PRGid_2_BAMid[PRGid] = IDforBAMIndex;
			BAMid_2_PRGid[IDforBAMIndex] = PRGid;

			if(inThisBAM)
			{
				interestingIntervals[IDforBAMIndex].push_back(oneInterestingInterval(PRGid, IDforBAMIndex, startIndex_1based, stopIndex_1based));
				interestingIntervals_n++;
			}
		}

		std::cout << "processBAM::initBAM(..): Will examine " << interestingIntervals_n << " PRG-mapped intervals.\n";

		_currentBAM = BAM;
		_currentBAM_isExtendedReferenceGenome = extendedReferenceGenome;
	}
	_currentBAM1 = "";
	_currentBAM2 = "";
}

void processBAM::initBAMs(std::string BAM1, std::string BAM2, bool extendedReferenceGenome)
{
	if((_currentBAM1 != BAM1) || (_currentBAM2 != BAM2) || (_currentBAM_isExtendedReferenceGenome != extendedReferenceGenome))
	{
		assert(Utilities::fileExists(BAM1));
		assert(Utilities::fileExists(BAM2));
		R1.Open(BAM1);
		R2.Open(BAM2);

		R1.LocateIndex();
		if ( ! R1.HasIndex() )
		{
			throw std::runtime_error("File "+BAM1+" does not seem to be indexed - please specify indexed BAM!");
		}

		R2.LocateIndex();
		if ( ! R2.HasIndex() )
		{
			throw std::runtime_error("File "+BAM2+" does not seem to be indexed - please specify indexed BAM!");
		}

		std::string extendedReferenceGenomePath;
		if(Utilities::fileExists(graphDir + "/extendedReferenceGenomePath.txt"))
		{
			extendedReferenceGenomePath  = Utilities::getFirstLine(graphDir + "/extendedReferenceGenomePath.txt");
		}
		else
		{
			extendedReferenceGenomePath  = graphDir + "/extendedReferenceGenome/extendedReferenceGenome.fa";
			assert(Utilities::fileExists(extendedReferenceGenomePath));
		}
	    assert(Utilities::fileExists(extendedReferenceGenomePath));
	    std::string PRGonlyReferenceGenomePath = graphDir + "/mapping_PRGonly/referenceGenome.fa";

		// sequences
		std::string sequencesFile = graphDir + "/sequences.txt";

		std::ifstream sequencesStream;
		sequencesStream.open(sequencesFile.c_str());
		assert(sequencesStream.is_open());

		std::string headerLine;
		std::string line;
		assert(sequencesStream.good());
		std::getline(sequencesStream, headerLine);
		Utilities::eraseNL(headerLine);
		std::vector<std::string> headerFields = Utilities::split(headerLine, "\t");

		size_t interestingIntervals_n = 0;

		PRGid_2_BAMid.clear();
		BAMid_2_PRGid.clear();
		interestingIntervals.clear();


		while(sequencesStream.good())
		{
			std::getline(sequencesStream, line);
			Utilities::eraseNL(line);
			if(line.length() == 0)
				continue;

			std::vector<std::string> lineFields = Utilities::split(line, "\t");
			assert(lineFields.size() == headerFields.size());

			std::map<std::string, std::string> L;
			for(unsigned int i = 0; i < headerFields.size(); i++)
			{
				L[headerFields.at(i)] = lineFields.at(i);
			}

			std::string IDforBAMIndex;
			int PRGid = Utilities::StrtoI(L.at("SequenceID"));

			int startIndex_1based;
			int stopIndex_1based;
			if(L.at("Chr") != "")
			{
				IDforBAMIndex = L.at("Chr");

				if(extendedReferenceGenomeSequences.count(IDforBAMIndex) == 0)
				{
					throw std::runtime_error(IDforBAMIndex + " cannot be found in the extended reference genome "+extendedReferenceGenomePath);
				}
				if(PRGonlyReferenceGenomeSequences.count(IDforBAMIndex) == 0)
				{
					throw std::runtime_error(IDforBAMIndex + " cannot be found in the PRG-only reference genome "+PRGonlyReferenceGenomePath);
				}


				if(R1.GetReferenceID(IDforBAMIndex) == -1)
				{
					if(! Utilities::sequence_all_Ns(extendedReferenceGenomeSequences.at(IDforBAMIndex)))
					{
						throw std::runtime_error(IDforBAMIndex + " cannot be found in the index of "+BAM1);
					}
					else
					{
						continue;
					}
				}

				if(R2.GetReferenceID(IDforBAMIndex) == -1)
				{
					if(! Utilities::sequence_all_Ns(extendedReferenceGenomeSequences.at(IDforBAMIndex)))
					{
						throw std::runtime_error(IDforBAMIndex + " cannot be found in the index of "+BAM2);
					}
					else
					{
						continue;
					}					
				}


				int internalBAMid = R1.GetReferenceID(IDforBAMIndex);
				assert(R1.GetReferenceID(IDforBAMIndex) == R2.GetReferenceID(IDforBAMIndex));

				BamTools::RefData internalBAMSpecification1 = R1.GetReferenceData().at(internalBAMid);
				BamTools::RefData internalBAMSpecification2 = R2.GetReferenceData().at(internalBAMid);
				assert(internalBAMSpecification1.RefName == IDforBAMIndex);
				assert(internalBAMSpecification1.RefLength == internalBAMSpecification2.RefLength);
				assert(internalBAMSpecification1.RefName == internalBAMSpecification2.RefName);

				if(extendedReferenceGenome)
				{
					startIndex_1based = Utilities::StrtoI(L.at("Start_1based"));
					stopIndex_1based = Utilities::StrtoI(L.at("Stop_1based"));
					assert(startIndex_1based <= internalBAMSpecification1.RefLength);
					assert(stopIndex_1based <= internalBAMSpecification1.RefLength);
				}
				else
				{
					startIndex_1based = 1;
					stopIndex_1based = PRGonlyReferenceGenomeSequences.at(IDforBAMIndex).length();
					assert(stopIndex_1based == internalBAMSpecification1.RefLength);
				}
			}
			else
			{
				IDforBAMIndex = "PRG_" + L.at("SequenceID");
				int internalBAMid = R1.GetReferenceID(IDforBAMIndex);
				assert(R1.GetReferenceID(IDforBAMIndex) == R2.GetReferenceID(IDforBAMIndex));

				if(extendedReferenceGenomeSequences.count(IDforBAMIndex) == 0)
				{
					throw std::runtime_error(IDforBAMIndex + " cannot be found in the extended reference genome "+extendedReferenceGenomePath);
				}
				if(PRGonlyReferenceGenomeSequences.count(IDforBAMIndex) == 0)
				{
					throw std::runtime_error(IDforBAMIndex + " cannot be found in the PRG-only reference genome "+PRGonlyReferenceGenomePath);
				}

				if(R1.GetReferenceID(IDforBAMIndex) == -1)
				{
					if(! Utilities::sequence_all_Ns(extendedReferenceGenomeSequences.at(IDforBAMIndex)))
					{
						throw std::runtime_error(IDforBAMIndex + " cannot be found in the index of "+BAM1);
					}
					else
					{
						continue;
					}					
				}

				if(R2.GetReferenceID(IDforBAMIndex) == -1)
				{
					if(! Utilities::sequence_all_Ns(extendedReferenceGenomeSequences.at(IDforBAMIndex)))
					{
						throw std::runtime_error(IDforBAMIndex + " cannot be found in the index of "+BAM2);
					}
				}

				BamTools::RefData internalBAMSpecification1 = R1.GetReferenceData().at(internalBAMid);
				BamTools::RefData internalBAMSpecification2 = R2.GetReferenceData().at(internalBAMid);

				assert(internalBAMSpecification1.RefName == IDforBAMIndex);

				assert(internalBAMSpecification1.RefLength == internalBAMSpecification2.RefLength);
				assert(internalBAMSpecification1.RefName == internalBAMSpecification2.RefName);

				startIndex_1based = 1;
				stopIndex_1based = extendedReferenceGenomeSequences.at(IDforBAMIndex).length();
				assert(stopIndex_1based > startIndex_1based);
				if(!(stopIndex_1based == internalBAMSpecification1.RefLength))
				{
					std::cerr << "Problem in sequences " << sequencesFile << "! for " << "PRG_" + L.at("SequenceID") << "\n";
					std::cerr << "\t" << "BAM" << ": " <<  BAM1 << "\n";
					std::cerr << "\t" << "stopIndex_1based" << ": " << stopIndex_1based << "\n";
					std::cerr << "\t" << "startIndex_1based" << ": " << startIndex_1based << "\n";
					std::cerr << "\t" << "internalBAMSpecification.RefLength" << ": " << internalBAMSpecification1.RefLength << "\n";
					std::cerr << std::flush;
				}
				assert(stopIndex_1based == internalBAMSpecification1.RefLength);
			}

			assert(R1.GetReferenceID(IDforBAMIndex) != -1);
			assert(extendedReferenceGenomeSequences.count(IDforBAMIndex));
			
			PRGid_2_BAMid[PRGid] = IDforBAMIndex;
			BAMid_2_PRGid[IDforBAMIndex] = PRGid;

			interestingIntervals[IDforBAMIndex].push_back(oneInterestingInterval(PRGid, IDforBAMIndex, startIndex_1based, stopIndex_1based));
			interestingIntervals_n++;
		}

		std::cout << "processBAM::initBAM(..): Will examine " << interestingIntervals_n << " PRG-mapped intervals.\n";

		_currentBAM1 = BAM1;
		_currentBAM2 = BAM2;
		_currentBAM_isExtendedReferenceGenome = extendedReferenceGenome;
	}
	_currentBAM = "";
}

bool processBAM::alignment_touches_gene(const std::tuple<std::string, int,  BamTools::BamAlignment, int>& al, const std::vector<std::string>& levels) const
{
	std::set<std::string> alignment_segments = alignment_get_segments(al, levels);
	bool touchesGene = false;
	for(auto S : alignment_segments)
	{
		if((S.find("intron") != std::string::npos) || (S.find("exon") != std::string::npos))
		{
			touchesGene = true;
		}
	}
	return touchesGene;
}
std::set<std::string> processBAM::alignment_get_segments(const std::tuple<std::string, int,  BamTools::BamAlignment, int>& al, const std::vector<std::string>& levels) const
{
	std::set<std::string> forReturn;
	std::pair<int, int> start_stop = alignment_get_startstop_PRGcoordinates(al);
	assert(start_stop.first <= start_stop.second);
	for(unsigned int level = start_stop.first; level <= start_stop.second; level++)
	{
		std::string levelID = levels.at(level);
		std::vector<std::string> parts = Utilities::split(levelID, "_");
		assert(parts.size() >= 3);
		std::vector<std::string> parts_for_print;
		for(unsigned int i = 1; i < (parts.size() - 1); i++)
		{
			parts_for_print.push_back(parts.at(i));
		}
		std::string id_for_print = Utilities::join(parts_for_print, "_");
		forReturn.insert(id_for_print);
	}
	return forReturn;
}

void processBAM::protoSeedStatistics(const std::map<std::string, reads::protoSeeds>& seeds, aligner::statistics* statisticsStore) const
{
	std::set<std::string> incompleteSeeds;

	double complete_protoSeeds_chains_sum = 0;
	double complete_protoSeeds_chainLength_sum = 0;
	int complete_protoSeeds_chainLength_sum_n = 0;

	double complete_protoSeeds_primaryChainLength_sum = 0;
	int complete_protoSeeds_primaryChainLength_sum_n = 0;

	std::vector<std::string> completeProtoSeeds;
	for(std::map<std::string, reads::protoSeeds>::const_iterator seedIt = seeds.begin(); seedIt != seeds.end(); seedIt++)
	{
		if(! seedIt->second.isComplete())
		{
			incompleteSeeds.insert(seedIt->first);
			// seedIt->second.printDebug(R);
		}
		else
		{
			completeProtoSeeds.push_back(seedIt->first);

			int n_chains = 0;
			complete_protoSeeds_chains_sum += (seedIt->second.read1_alignments.size() + seedIt->second.read2_alignments.size());

			size_t r1_primary = seedIt->second.read1_getPrimaryAlignmentI();
			size_t r2_primary = seedIt->second.read2_getPrimaryAlignmentI();

			for(size_t r1_i = 0; r1_i < seedIt->second.read1_alignments.size(); r1_i++)
			{
				const BamTools::BamAlignment& alignment = std::get<2>(seedIt->second.read1_alignments.at(r1_i));
				int start = alignment.Position;
				int stop = alignment.GetEndPosition(false, true);
				assert(stop >= start);
				int L = stop - start + 1;

				complete_protoSeeds_chainLength_sum += L;
				complete_protoSeeds_chainLength_sum_n++;

				if(r1_i == r1_primary)
				{
					complete_protoSeeds_primaryChainLength_sum += L;
					complete_protoSeeds_primaryChainLength_sum_n++;
				}
			}

			for(size_t r2_i = 0; r2_i < seedIt->second.read2_alignments.size(); r2_i++)
			{
				const BamTools::BamAlignment& alignment = std::get<2>(seedIt->second.read2_alignments.at(r2_i));
				int start = alignment.Position;
				int stop = alignment.GetEndPosition(false, true);
				assert(stop >= start);
				int L = stop - start + 1;

				complete_protoSeeds_chainLength_sum += L;
				complete_protoSeeds_chainLength_sum_n++;

				if(r2_i == r2_primary)
				{
					complete_protoSeeds_primaryChainLength_sum += L;
					complete_protoSeeds_primaryChainLength_sum_n++;
				}
			}
		}
	}

	std::cout << Utilities::timestamp() << "Proto-seed statistics:\n";
	std::cout << "\t" << "Incomplete" << ": " << incompleteSeeds.size() << "\n";
	std::cout << "\t" << "Complete" << ": " << completeProtoSeeds.size() << "\n";
	std::cout << "\t\t" << "Average chains per read" << ": " << (complete_protoSeeds_chains_sum / (2.0 * completeProtoSeeds.size())) << "\n";
	std::cout << "\t\t" << "Average chain length" << ": " << (complete_protoSeeds_chainLength_sum / double(complete_protoSeeds_chainLength_sum_n)) << "\n";
	std::cout << "\t\t" << "Average primary chain length" << ": " << (complete_protoSeeds_primaryChainLength_sum / double(complete_protoSeeds_primaryChainLength_sum_n)) << "\n";
	std::cout << "\n" << std::flush;

	if(statisticsStore != 0)
	{
		statisticsStore->incomplete_seeds += incompleteSeeds.size();
		statisticsStore->complete_seeds += completeProtoSeeds.size();
		statisticsStore->complete_seeds_n_chains += complete_protoSeeds_chains_sum;
		statisticsStore->complete_seeds_chains_combined_length += complete_protoSeeds_chainLength_sum;
		statisticsStore->complete_seeds_n_primaryChains += complete_protoSeeds_primaryChainLength_sum_n;
		statisticsStore->complete_seeds_primaryChain_combined_length += complete_protoSeeds_primaryChainLength_sum;
	}
}

void processBAM::alignReads(std::string BAM, simulator::trueReadLevels* trueReadLevels, double insertSize_mean, double insertSize_sd, std::string outputDirectory, bool extendedReferenceGenome, hla::HLATyper* HLATyper, int threads)
{
	initBAM(BAM, extendedReferenceGenome, false);
	std::map<std::string, reads::protoSeeds> seeds = extractSeeds();
	alignReads_postSeedExtraction(seeds, trueReadLevels, insertSize_mean, insertSize_sd, outputDirectory, HLATyper, threads);
}

void processBAM::alignReads_and_inferHLA(std::string BAM, simulator::trueReadLevels* trueReadLevels, double insertSize_mean, double insertSize_sd, std::string outputDirectory, bool extendedReferenceGenome, hla::HLATyper* HLATyper, int threads)
{
	initBAM(BAM, extendedReferenceGenome, false);

	std::vector<std::string> readIDs = getReadIDs();

	size_t chunks_per_segment = 10000;
	unsigned int segments = readIDs.size() / chunks_per_segment;
	if(segments == 0)
		segments = 1;

	std::vector<std::vector<int>> bases_per_level_perThread;
	std::vector<std::vector<mapper::reads::oneReadPair>> HLA_raw_reads_perThread;
	std::vector<std::vector<mapper::reads::verboseSeedChainPair>> HLA_alignments_perThread;

	std::string outputDirectory_for_HLA = outputDirectory + "/hla/";
	Utilities::make_or_clearDirectory(outputDirectory_for_HLA);


	size_t processed_read_pairs = 0;

	bases_per_level_perThread.resize(threads);
	HLA_raw_reads_perThread.resize(threads);
	HLA_alignments_perThread.resize(threads);

	for(int threadI = 0; threadI < threads; threadI++)
	{
		bases_per_level_perThread.at(threadI).resize(g->NodesPerLevel.size() - 1, 0);
	}

	chrono::milliseconds ms_before = chrono::duration_cast< chrono::milliseconds >(
		chrono::system_clock::now().time_since_epoch() );

	aligner::statistics statisticsStore;

	for(unsigned int segmentI = 0; segmentI < segments; segmentI++)
	{
		size_t segment_start = segmentI * chunks_per_segment;
		size_t segment_stop = (segmentI+1)*chunks_per_segment - 1;
		if(segmentI == (segments - 1))
		{
			segment_stop = readIDs.size() - 1;
		}
		assert(segment_stop >= segment_start);
		std::vector<std::string> readIDs_thisSegment(readIDs.begin()+segment_start, readIDs.begin()+segment_stop+1);

		std::cout << Utilities::timestamp() << "Process " << segmentI << "/" << segments << "\n" << std::flush;

		
		std::cout << Utilities::timestamp() << "\tStart seed extraction\n" << std::flush;
				
		std::map<std::string, reads::protoSeeds> seeds = extractSeeds2(std::set<std::string>(readIDs_thisSegment.begin(), readIDs_thisSegment.end()));

		std::cout << Utilities::timestamp() << "\tAlignment\n" << std::flush;

		processed_read_pairs += alignReads_postSeedExtraction_andStoreInto(seeds, trueReadLevels, insertSize_mean, insertSize_sd, outputDirectory, HLATyper, threads, bases_per_level_perThread, HLA_raw_reads_perThread, HLA_alignments_perThread, &statisticsStore);
		std::cout << Utilities::timestamp() << ".. done. Processed  " << processed_read_pairs << " in total (total read IDs: " << readIDs.size() << ".\n" << std::flush;
	}

	statisticsStore.printStatistics();

	std::vector<mapper::reads::oneReadPair> HLA_raw_reads;
	std::vector<mapper::reads::verboseSeedChainPair> HLA_alignments;
	std::vector<int> bases_per_level;
	bases_per_level.resize(g->NodesPerLevel.size() - 1, 0);

	if(HLATyper != 0)
	{
		outputDirectory_for_HLA = outputDirectory + "/hla/";
		Utilities::make_or_clearDirectory(outputDirectory_for_HLA);
	}

	for(int threadI = 0; threadI < threads; threadI++)
	{
		HLA_raw_reads.insert(HLA_raw_reads.end(), HLA_raw_reads_perThread.at(threadI).begin(), HLA_raw_reads_perThread.at(threadI).end());
		HLA_alignments.insert(HLA_alignments.end(), HLA_alignments_perThread.at(threadI).begin(), HLA_alignments_perThread.at(threadI).end());
		for(unsigned int levelI = 0; levelI < (g->NodesPerLevel.size() - 1); levelI++)
		{
			bases_per_level.at(levelI) += bases_per_level_perThread.at(threadI).at(levelI);
		}
	}

	chrono::milliseconds ms_after= chrono::duration_cast< chrono::milliseconds >(
			chrono::system_clock::now().time_since_epoch()
	);

	chrono::milliseconds ms_duration = ms_after - ms_before;
	double protoSeeds_per_ms =  processed_read_pairs / (double)ms_duration.count();
	double protoSeeds_per_s = protoSeeds_per_ms * 1000.0;

	std::cout << Utilities::timestamp() << "Processed " << processed_read_pairs << " protoSeeds (read pairs)" << "\n" << std::flush;
	std::cout << "Speed: " << protoSeeds_per_s << " protoSeeds (read pairs) per s" << "\n" << std::flush;

	//assert(1 == 0);
	
	std::string output_levels_fn = outputDirectory + "/reads_per_level.txt";
	std::ofstream levels_stream;
	levels_stream.open(output_levels_fn.c_str());
	assert(levels_stream.is_open());

	for(unsigned int lI = 0; lI < bases_per_level.size(); lI++)
	{
		std::string levelName = g->getOneLocusIDforLevel(lI);
		levels_stream << lI << "\t"
					  << levelName << "\t"
					  << bases_per_level.at(lI) << "\n";
	}


	std::cout << Utilities::timestamp() << "Initiate HLA typing!\n" << std::flush;
	std::vector<mapper::reads::oneRead> empty_rawUnpairedReads;
	std::vector<mapper::reads::verboseSeedChain> empty_alignedUnpairedReads;
	std::cout << "Call HLA typing with " << HLA_alignments.size() << " alignments.\n" << std::flush;
	omp_set_num_threads(threads_for_HLAtyping);
	HLATyper->HLATypeInference(HLA_raw_reads, HLA_alignments, empty_rawUnpairedReads, empty_alignedUnpairedReads, insertSize_mean, insertSize_sd, outputDirectory_for_HLA);
	omp_set_num_threads(threads);

}


void processBAM::alignReadsMulti(std::vector<std::string> BAMs, simulator::trueReadLevels* trueReadLevels, double insertSize_mean, double insertSize_sd, std::string outputDirectory, bool extendedReferenceGenome, hla::HLATyper* HLATyper, int threads)
{
	std::map<std::string, reads::protoSeeds> seeds;
	for(unsigned int bamI = 0; bamI < BAMs.size(); bamI++)
	{
		std::string BAM = BAMs.at(bamI);
		initBAM(BAM, extendedReferenceGenome, true);
		extractSeedsInto(seeds, bamI);
	}
	alignReads_postSeedExtraction(seeds, trueReadLevels, insertSize_mean, insertSize_sd, outputDirectory, HLATyper, threads);
}

void processBAM::alignReads2(std::string BAM1, std::string BAM2, simulator::trueReadLevels* trueReadLevels, double insertSize_mean, double insertSize_sd, std::string outputDirectory, bool extendedReferenceGenome, hla::HLATyper* HLATyper, int threads)
{
	initBAMs(BAM1, BAM2, extendedReferenceGenome);
	std::map<std::string, reads::protoSeeds> seeds = extractSeeds_2BAMs();
	alignReads_postSeedExtraction(seeds, trueReadLevels, insertSize_mean, insertSize_sd, outputDirectory, HLATyper, threads);
}

void processBAM::sortChainsInSeeds(std::map<std::string, reads::protoSeeds>& seeds) const
{
	for(auto & read : seeds)
	{ 
		reads::protoSeeds& seeds = read.second;

		auto sort_alignments = [&] (std::vector<std::tuple<std::string, int, BamTools::BamAlignment, int>>& alignments) -> void {
			std::sort(alignments.begin(), alignments.end(),
				[&](const std::tuple<std::string, int, BamTools::BamAlignment, int>& a, const std::tuple<std::string, int, BamTools::BamAlignment, int>& b)
				{
					int score_a = getAlignmentScore(std::get<2>(a));
					int score_b = getAlignmentScore(std::get<2>(b));			
					return (score_a < score_b);
				}
			
			);
			std::reverse(alignments.begin(), alignments.end());
		};
		
		sort_alignments(seeds.read1_alignments);
		sort_alignments(seeds.read2_alignments);
	}
}


void processBAM::reduceNonGeneSeeds(std::map<std::string, reads::protoSeeds>& seeds) const
{
	for(auto & read : seeds)
	{ 
		reads::protoSeeds& seeds = read.second;

		std::vector<std::tuple<std::string, int, BamTools::BamAlignment, int>> new_read1_alignments;
		std::vector<std::tuple<std::string, int, BamTools::BamAlignment, int>> new_read2_alignments;

		for(auto alignment : seeds.read1_alignments)
		{
			if((std::get<3>(alignment) == 0) || alignment_touches_gene(alignment, g_level_names))
			{
				new_read1_alignments.push_back(alignment);
			}
		}
		for(auto alignment : seeds.read2_alignments)
		{
			if((std::get<3>(alignment) == 0) || alignment_touches_gene(alignment, g_level_names))
			{
				new_read2_alignments.push_back(alignment);
			} 
		}

		// std::cout << "Reduce from " << seeds.read1_alignments.size() << " to " << new_read1_alignments.size() << " and from " <<  seeds.read2_alignments.size() << " to " << new_read2_alignments.size() << "\n";
		seeds.read1_alignments = new_read1_alignments;
		seeds.read2_alignments = new_read2_alignments;

		seeds.refreshPrimaryStatus();

	}
}


void processBAM::alignReads_postSeedExtraction(std::map<std::string, reads::protoSeeds>& seeds, simulator::trueReadLevels* trueReadLevels, double insertSize_mean, double insertSize_sd, std::string outputDirectory, hla::HLATyper* HLATyper, int threads)
{

	boost::math::normal rnd_InsertSize(insertSize_mean, insertSize_sd);
	double max_insertsize_penalty = boost::math::pdf(rnd_InsertSize, insertSize_mean + 8 * insertSize_sd);
	assert(max_insertsize_penalty > 0);
	assert(max_insertsize_penalty <= 1);
	double max_insertsize_penalty_log = log(max_insertsize_penalty);


	protoSeedStatistics(seeds);

	reduceNonGeneSeeds(seeds);
	
	sortChainsInSeeds(seeds);
	
	std::cout << "Reduced the number of seeds in non-gene areas - statistics:\n" << std::flush;

	protoSeedStatistics(seeds);

	std::set<std::string> incompleteSeeds;
	std::vector<std::string> completeProtoSeeds;
	for(std::map<std::string, reads::protoSeeds>::const_iterator seedIt = seeds.begin(); seedIt != seeds.end(); seedIt++)
	{
		if(! seedIt->second.isComplete())
		{
			incompleteSeeds.insert(seedIt->first);
			// seedIt->second.printDebug((_currentBAM1 != "") ? R1 : R);
			// debug
			// assert(incompleteSeeds.size() <= 10);
		}
		else
		{
			completeProtoSeeds.push_back(seedIt->first);
		}
	}

	std::cout << "processBAM::alignReads(): Deal " << seeds.size() << " total read pairs with seeds, of which " << incompleteSeeds.size() << " are incomplete.\n" << std::flush;


	chrono::milliseconds ms_before = chrono::duration_cast< chrono::milliseconds >(
		chrono::system_clock::now().time_since_epoch() );

	assert(threads >= 1);

	std::cerr << "threads: " << threads << "\n" << std::flush;

	std::string outputDirectory_for_HLA;
	std::vector<std::vector<mapper::reads::oneReadPair>> HLA_raw_reads_perThread;
	std::vector<std::vector<mapper::reads::verboseSeedChainPair>> HLA_alignments_perThread;

	if(HLATyper != 0)
	{
		outputDirectory_for_HLA = outputDirectory + "/hla/";
		Utilities::make_or_clearDirectory(outputDirectory_for_HLA);
	}

	std::vector<std::vector<int>> bases_per_level_perThread;

	bases_per_level_perThread.resize(threads);
	HLA_raw_reads_perThread.resize(threads);
	HLA_alignments_perThread.resize(threads);

	for(int threadI = 0; threadI < threads; threadI++)
	{
		bases_per_level_perThread.at(threadI).resize(g->NodesPerLevel.size() - 1, 0);
	}

	omp_set_num_threads(threads);
	eA->init_for_threads(threads);

	size_t completeProtoSeeds_size = completeProtoSeeds.size();
	// #pragma omp parallel for
	for(size_t protoSeedI = 0; protoSeedI < completeProtoSeeds_size; protoSeedI++)
	{
		// std::cout << Utilities::timestamp() << "Processing " << protoSeedI << " / " << completeProtoSeeds_size << "\n" << std::flush;
		
		if((protoSeedI % 100000) == 0)
		{
			std::cout << Utilities::timestamp() << "Read pair " << protoSeedI << " of " << completeProtoSeeds_size << "\n" << std::flush;
		}
		int threadI = omp_get_thread_num();
		assert(threadI < threads);
		

		std::string seedID = completeProtoSeeds.at(protoSeedI);
		const reads::protoSeeds& protoSeed = seeds.at(seedID);

		//if(protoSeedI < 140870)
		//	continue;
		
		reads::verboseSeedChainPair alignment = alignOneReadPair(protoSeed, rnd_InsertSize, max_insertsize_penalty_log, trueReadLevels);

		for(unsigned int aI = 0; aI < alignment.chains.first.graph_aligned_levels.size(); aI++)
		{
			int level = alignment.chains.first.graph_aligned_levels.at(aI);
			if((level != -1) && (alignment.chains.first.graph_aligned.at(aI) != '_'))
			{
				bases_per_level_perThread.at(threadI).at(level)++;
			}
		} 
		for(unsigned int aI = 0; aI < alignment.chains.second.graph_aligned_levels.size(); aI++)
		{
			int level = alignment.chains.second.graph_aligned_levels.at(aI);
			if((level != -1) && (alignment.chains.second.graph_aligned.at(aI) != '_'))
			{
				bases_per_level_perThread.at(threadI).at(level)++;
			}
		}

		bool includeInHLA = false;

		std::pair<int, int> read1_graphLevels = make_pair(alignment.chains.first.alignment_firstLevel(), alignment.chains.first.alignment_lastLevel());
		std::pair<int, int> read2_graphLevels = make_pair(alignment.chains.second.alignment_firstLevel(), alignment.chains.second.alignment_lastLevel());

		if(HLATyper != 0)
		{
			if(read1_graphLevels.first != -1)
			{
				assert(read1_graphLevels.first <= read1_graphLevels.second);
				includeInHLA = ( includeInHLA || (HLATyper->intervalOverlapsWithGenes(read1_graphLevels.first, read1_graphLevels.second)) );
			}

			if(read2_graphLevels.first != -1)
			{
				assert(read2_graphLevels.first <= read2_graphLevels.second);
				includeInHLA = ( includeInHLA || (HLATyper->intervalOverlapsWithGenes(read2_graphLevels.first, read2_graphLevels.second)) );
			}
		}
		if(includeInHLA)
		{
			size_t r1_primary = protoSeed.read1_getPrimaryAlignmentI();
			size_t r2_primary = protoSeed.read2_getPrimaryAlignmentI();

			reads::oneRead r1(
				std::get<2>(protoSeed.read1_alignments.at(r1_primary)).Name,
				std::get<2>(protoSeed.read1_alignments.at(r1_primary)).QueryBases,
				std::get<2>(protoSeed.read1_alignments.at(r1_primary)).Qualities
			);
			reads::oneRead r2(
				std::get<2>(protoSeed.read2_alignments.at(r2_primary)).Name,
				std::get<2>(protoSeed.read2_alignments.at(r2_primary)).QueryBases,
				std::get<2>(protoSeed.read2_alignments.at(r2_primary)).Qualities
			);

			if(std::get<2>(protoSeed.read1_alignments.at(r1_primary)).IsReverseStrand())
			{
				r1.invert();
			}

			if(std::get<2>(protoSeed.read2_alignments.at(r2_primary)).IsReverseStrand())
			{
				r2.invert();
			}

			reads::oneReadPair rP(r1, r2, 0);

			HLA_raw_reads_perThread.at(threadI).push_back(rP);
			HLA_alignments_perThread.at(threadI).push_back(alignment);
		}
	}

	std::cerr << "Done\n" << std::flush;

	std::vector<mapper::reads::oneReadPair> HLA_raw_reads;
	std::vector<mapper::reads::verboseSeedChainPair> HLA_alignments;
	std::vector<int> bases_per_level;
	bases_per_level.resize(g->NodesPerLevel.size() - 1, 0);

	if(HLATyper != 0)
	{
		outputDirectory_for_HLA = outputDirectory + "/hla/";
		Utilities::make_or_clearDirectory(outputDirectory_for_HLA);
	}

	for(int threadI = 0; threadI < threads; threadI++)
	{
		HLA_raw_reads.insert(HLA_raw_reads.end(), HLA_raw_reads_perThread.at(threadI).begin(), HLA_raw_reads_perThread.at(threadI).end());
		HLA_alignments.insert(HLA_alignments.end(), HLA_alignments_perThread.at(threadI).begin(), HLA_alignments_perThread.at(threadI).end());
		for(unsigned int levelI = 0; levelI < (g->NodesPerLevel.size() - 1); levelI++)
		{
			bases_per_level.at(levelI) += bases_per_level_perThread.at(threadI).at(levelI);
		}
	}

	chrono::milliseconds ms_after= chrono::duration_cast< chrono::milliseconds >(
			chrono::system_clock::now().time_since_epoch()
	);
	chrono::milliseconds ms_duration = ms_after - ms_before;
	double protoSeeds_per_ms = (double)completeProtoSeeds.size() / (double)ms_duration.count();
	double protoSeeds_per_s = protoSeeds_per_ms * 1000.0;

	std::cout << Utilities::timestamp() << "Processed " << completeProtoSeeds.size() << " protoSeeds (read pairs)" << "\n" << std::flush;
	std::cout << "Speed: " << protoSeeds_per_s << " protoSeeds (read pairs) per s" << "\n" << std::flush;

	std::string output_levels_fn = outputDirectory + "/reads_per_level.txt";
	std::ofstream levels_stream;
	levels_stream.open(output_levels_fn.c_str());
	assert(levels_stream.is_open());

	for(unsigned int lI = 0; lI < bases_per_level.size(); lI++)
	{
		std::string levelName = g->getOneLocusIDforLevel(lI);
		levels_stream << lI << "\t"
					  << levelName << "\t"
					  << bases_per_level.at(lI) << "\n";
	}

	if(HLATyper != 0)
	{
		std::cout << Utilities::timestamp() << "Initiate HLA typing!\n" << std::flush;
		std::vector<mapper::reads::oneRead> empty_rawUnpairedReads;
		std::vector<mapper::reads::verboseSeedChain> empty_alignedUnpairedReads;
		std::cout << "Call HLA typing with " << HLA_alignments.size() << " alignments.\n" << std::flush;
		omp_set_num_threads(threads_for_HLAtyping);
		HLATyper->HLATypeInference(HLA_raw_reads, HLA_alignments, empty_rawUnpairedReads, empty_alignedUnpairedReads, insertSize_mean, insertSize_sd, outputDirectory_for_HLA);
		omp_set_num_threads(threads);
	}
}

size_t processBAM::alignReads_postSeedExtraction_andStoreInto(std::map<std::string, reads::protoSeeds>& seeds, simulator::trueReadLevels* trueReadLevels, double insertSize_mean, double insertSize_sd, std::string outputDirectory, const hla::HLATyper* HLATyper, int threads, std::vector<std::vector<int>>& bases_per_level_perThread, std::vector<std::vector<mapper::reads::oneReadPair>>& HLA_raw_reads_perThread, std::vector<std::vector<mapper::reads::verboseSeedChainPair>>& HLA_alignments_perThread, 	aligner::statistics* statisticsStore)
{
	boost::math::normal rnd_InsertSize(insertSize_mean, insertSize_sd);
	double max_insertsize_penalty = boost::math::pdf(rnd_InsertSize, insertSize_mean + 8 * insertSize_sd);
	assert(max_insertsize_penalty > 0);
	assert(max_insertsize_penalty <= 1);
	double max_insertsize_penalty_log = log(max_insertsize_penalty);

	protoSeedStatistics(seeds, statisticsStore);

	// reduceNonGeneSeeds(seeds);

	sortChainsInSeeds(seeds);

	// std::cout << "Reduced the number of seeds in non-gene areas - statistics:\n" << std::flush;

	// protoSeedStatistics(seeds);

	std::set<std::string> incompleteSeeds;
	std::vector<std::string> completeProtoSeeds;
	for(std::map<std::string, reads::protoSeeds>::const_iterator seedIt = seeds.begin(); seedIt != seeds.end(); seedIt++)
	{
		if(! seedIt->second.isComplete())
		{
			incompleteSeeds.insert(seedIt->first);
			// seedIt->second.printDebug((_currentBAM1 != "") ? R1 : R);
			// debug
			// assert(incompleteSeeds.size() <= 10);
		}
		else
		{
			completeProtoSeeds.push_back(seedIt->first);
		}
	}

	std::cout << "processBAM::alignReads_postSeedExtraction_andStoreInto(): Deal " << seeds.size() << " total read pairs with seeds, of which " << incompleteSeeds.size() << " are incomplete.\n" << std::flush;


	assert(threads >= 1);

	std::cerr << "threads: " << threads << "\n" << std::flush;


	omp_set_num_threads(threads);
	eA->init_for_threads(threads);

	size_t completeProtoSeeds_size = completeProtoSeeds.size();
	// #pragma omp parallel for
	for(size_t protoSeedI = 0; protoSeedI < completeProtoSeeds_size; protoSeedI++)
	{
		std::cout << Utilities::timestamp() << "Processing " << protoSeedI << " / " << completeProtoSeeds_size << "\n" << std::flush;
		 
		if((protoSeedI % 10000) == 0)
		{
			std::cout << Utilities::timestamp() << "Read pair " << protoSeedI << " of " << completeProtoSeeds_size << "\n" << std::flush;
		}
		int threadI = omp_get_thread_num();
		assert(threadI < threads);

		std::string seedID = completeProtoSeeds.at(protoSeedI);
		const reads::protoSeeds& protoSeed = seeds.at(seedID);

		//if(protoSeedI < 140870)
		//	continue;

		aligner::statistics beforeCall = *statisticsStore;
		reads::verboseSeedChainPair alignment = alignOneReadPair(protoSeed, rnd_InsertSize, max_insertsize_penalty_log, trueReadLevels, statisticsStore);

		for(unsigned int aI = 0; aI < alignment.chains.first.graph_aligned_levels.size(); aI++)
		{
			int level = alignment.chains.first.graph_aligned_levels.at(aI);
			if((level != -1) && (alignment.chains.first.graph_aligned.at(aI) != '_'))
			{
				bases_per_level_perThread.at(threadI).at(level)++;
			}
		}
		for(unsigned int aI = 0; aI < alignment.chains.second.graph_aligned_levels.size(); aI++)
		{
			int level = alignment.chains.second.graph_aligned_levels.at(aI);
			if((level != -1) && (alignment.chains.second.graph_aligned.at(aI) != '_'))
			{
				bases_per_level_perThread.at(threadI).at(level)++; 
			}
		}

		bool includeInHLA = false;

		std::pair<int, int> read1_graphLevels = make_pair(alignment.chains.first.alignment_firstLevel(), alignment.chains.first.alignment_lastLevel());
		std::pair<int, int> read2_graphLevels = make_pair(alignment.chains.second.alignment_firstLevel(), alignment.chains.second.alignment_lastLevel());

		if(HLATyper != 0)
		{
			if(read1_graphLevels.first != -1)
			{
				assert(read1_graphLevels.first <= read1_graphLevels.second);
				includeInHLA = ( includeInHLA || (HLATyper->intervalOverlapsWithGenes(read1_graphLevels.first, read1_graphLevels.second)) );
			}

			if(read2_graphLevels.first != -1)
			{
				assert(read2_graphLevels.first <= read2_graphLevels.second);
				includeInHLA = ( includeInHLA || (HLATyper->intervalOverlapsWithGenes(read2_graphLevels.first, read2_graphLevels.second)) );
			}
		}

		if(includeInHLA)
		{
			size_t r1_primary = protoSeed.read1_getPrimaryAlignmentI();
			size_t r2_primary = protoSeed.read2_getPrimaryAlignmentI();

			reads::oneRead r1(
				std::get<2>(protoSeed.read1_alignments.at(r1_primary)).Name,
				std::get<2>(protoSeed.read1_alignments.at(r1_primary)).QueryBases,
				std::get<2>(protoSeed.read1_alignments.at(r1_primary)).Qualities
			);
			reads::oneRead r2(
				std::get<2>(protoSeed.read2_alignments.at(r2_primary)).Name,
				std::get<2>(protoSeed.read2_alignments.at(r2_primary)).QueryBases,
				std::get<2>(protoSeed.read2_alignments.at(r2_primary)).Qualities
			);

			if(std::get<2>(protoSeed.read1_alignments.at(r1_primary)).IsReverseStrand())
			{
				r1.invert();
			}

			if(std::get<2>(protoSeed.read2_alignments.at(r2_primary)).IsReverseStrand())
			{
				r2.invert();
			}

			reads::oneReadPair rP(r1, r2, 0);

			HLA_raw_reads_perThread.at(threadI).push_back(rP);
			HLA_alignments_perThread.at(threadI).push_back(alignment);

			statisticsStore->readPairs_used_for_HLAtyping++;

			statisticsStore->takeInHLARelatedDiff(&beforeCall);
		}
	}

	std::cerr << "Done\n" << std::flush;

	return completeProtoSeeds.size();
}


bool processBAM::PRGContigAlignment2Seed(Graph *g, const reads::PRGContigBAMAlignment& PRGcontigAlignment, bool paranoid, reads::verboseSeedChain& graphSeed, reads::verboseSeedChain& sequenceSeed, const std::vector<bool>& inGraphGapStretch)
{
	std::vector<int> graphSpace_graph_aligned_levels;
	std::string graphSpace_graph_aligned;
	std::string graphSpace_sequence_aligned;
	int firstLevel_graph = PRGcontigAlignment.graph_aligned_levels.front();
	int lastLevel_graph = PRGcontigAlignment.graph_aligned_levels.back();
	int totalLevels_graph = lastLevel_graph - firstLevel_graph + 1;
	int reserveLevels = (totalLevels_graph > (int)PRGcontigAlignment.graph_aligned_levels.size()) ? totalLevels_graph : PRGcontigAlignment.graph_aligned_levels.size();
	// std::cout << "reserveLevels: " << reserveLevels << "\n" << std::flush;

	bool verbose = false;
	if(verbose)
	{
		std::cout <<"Initial alignment:\n";
		std::cout << "\t" << "Levels" << ": " <<  Utilities::join(Utilities::ItoStr(PRGcontigAlignment.graph_aligned_levels), ", ") << "\n";
		std::cout << "\t" << "Graph" << ": " <<  PRGcontigAlignment.graph_aligned << "\n";
		std::cout << "\t" << "Sequence" << ": " <<  PRGcontigAlignment.sequence_aligned << "\n";
		std::cout << "\n" << std::flush;
	}
	
	graphSpace_graph_aligned_levels.reserve(reserveLevels * 1.1);
	graphSpace_graph_aligned.reserve(reserveLevels * 1.1);
	graphSpace_sequence_aligned.reserve(reserveLevels * 1.1);
	int graphSpace_sequence_aligned_startInRaw = PRGcontigAlignment.sequence_aligned_startInRaw;
	int graphSpace_sequence_aligned_stopInRaw = PRGcontigAlignment.sequence_aligned_stopInRaw;

	unsigned int firstColumn = 0;
	while(PRGcontigAlignment.graph_aligned_levels.at(firstColumn) == -1)
	{
		firstColumn++;
		graphSpace_sequence_aligned_startInRaw++;
	}
	assert(firstColumn < PRGcontigAlignment.graph_aligned_levels.size());

	unsigned int lastColumn = PRGcontigAlignment.graph_aligned_levels.size() - 1;
	while(PRGcontigAlignment.graph_aligned_levels.at(lastColumn) == -1)
	{
		lastColumn--;
		graphSpace_sequence_aligned_stopInRaw--;
	}
	assert(firstColumn < lastColumn);

	// std::cout << "firstColumn: " << firstColumn << "\n" << std::flush;

//	std::cout << "Build graph-space sequence\n";
	int lastInserted_graphLevel;
	for(unsigned int columnI = firstColumn; columnI <= lastColumn; columnI++)
	{
//		std::cout << "\tcolumnI: " << columnI << " / " << lastColumn << "\n";
//		std::cout << "\t\tlastInserted_graphLevel: " << lastInserted_graphLevel << "\n";

		if(columnI == firstColumn)
		{
			assert(PRGcontigAlignment.graph_aligned_levels.at(columnI) != -1);
			graphSpace_graph_aligned_levels.push_back(PRGcontigAlignment.graph_aligned_levels.at(columnI));
			graphSpace_graph_aligned.push_back(PRGcontigAlignment.graph_aligned.at(columnI));
			graphSpace_sequence_aligned.push_back(PRGcontigAlignment.sequence_aligned.at(columnI));
			lastInserted_graphLevel = PRGcontigAlignment.graph_aligned_levels.at(columnI);
		}
		else
		{
			if(PRGcontigAlignment.graph_aligned_levels.at(columnI) == -1)
			{
				graphSpace_graph_aligned_levels.push_back(PRGcontigAlignment.graph_aligned_levels.at(columnI));
				graphSpace_graph_aligned.push_back(PRGcontigAlignment.graph_aligned.at(columnI));
				graphSpace_sequence_aligned.push_back(PRGcontigAlignment.sequence_aligned.at(columnI));
			}
			else
			{
				if((lastInserted_graphLevel + 1) != PRGcontigAlignment.graph_aligned_levels.at(columnI))
				{
					assert((lastInserted_graphLevel + 1) < PRGcontigAlignment.graph_aligned_levels.at(columnI));
					for(int insertLevelI = lastInserted_graphLevel + 1; insertLevelI <= (PRGcontigAlignment.graph_aligned_levels.at(columnI) - 1); insertLevelI++)
					{
						graphSpace_graph_aligned_levels.push_back(insertLevelI);
						graphSpace_graph_aligned.push_back('_');
						graphSpace_sequence_aligned.push_back('_');
					}
				}

				graphSpace_graph_aligned_levels.push_back(PRGcontigAlignment.graph_aligned_levels.at(columnI));
				graphSpace_graph_aligned.push_back(PRGcontigAlignment.graph_aligned.at(columnI));
				graphSpace_sequence_aligned.push_back(PRGcontigAlignment.sequence_aligned.at(columnI));
				lastInserted_graphLevel = PRGcontigAlignment.graph_aligned_levels.at(columnI);

			}
		}
	}
	
	if(verbose)
	{
		std::cout << "Graph space alignment after initial transformation:\n" << std::flush;
		std::cout << "\t" << Utilities::join(Utilities::ItoStr(graphSpace_graph_aligned_levels), ", ") << "\n";
		std::cout << "\t" << graphSpace_graph_aligned << "\n";
		std::cout << "\t" << graphSpace_sequence_aligned << "\n";
		std::cout << "\n" << std::flush;	
	}
	
	cleanInitialAlignment(graphSpace_graph_aligned_levels, graphSpace_graph_aligned, graphSpace_sequence_aligned);

	if(verbose)
	{
		std::cout << "Graph space alignment after cleaning:\n" << std::flush;
		std::cout << "\t" << Utilities::join(Utilities::ItoStr(graphSpace_graph_aligned_levels), ", ") << "\n";
		std::cout << "\t" << graphSpace_graph_aligned << "\n";
		std::cout << "\t" << graphSpace_sequence_aligned << "\n";
		std::cout << "\n" << std::flush;	
	}
	
	int alignmentLength_before_restriction = graphSpace_graph_aligned_levels.size();
	restrictInitialAlignmentToNoGapAreas(graphSpace_graph_aligned_levels, graphSpace_graph_aligned, graphSpace_sequence_aligned, graphSpace_sequence_aligned_startInRaw, graphSpace_sequence_aligned_stopInRaw, inGraphGapStretch);
	int alignmentLength_after_restriction = graphSpace_graph_aligned_levels.size();
	int removed_columns_noGap_restriction = alignmentLength_before_restriction - alignmentLength_after_restriction;

	if(verbose)
	{
		std::cout << "Graph space alignment after no-gap restriction:\n" << std::flush;
		std::cout << "\t" << Utilities::join(Utilities::ItoStr(graphSpace_graph_aligned_levels), ", ") << "\n";
		std::cout << "\t" << graphSpace_graph_aligned << "\n";
		std::cout << "\t" << graphSpace_sequence_aligned << "\n";
		std::cout << "\n" << std::flush;	
	}
	
	size_t before_backtrace_matches = 0;
	size_t before_backtrace_length = 0;
	size_t after_backtrace_matches = 0;
	size_t after_backtrace_length = 0;
	for(unsigned int posI = 0; posI < graphSpace_graph_aligned.size(); posI++)
	{
		if(graphSpace_graph_aligned.at(posI) == graphSpace_sequence_aligned.at(posI))
		{
			before_backtrace_matches++;
		}
		before_backtrace_length++;
	}

//	std::cout << "graphSpace_graph_aligned:\n" << graphSpace_graph_aligned << "\n" << std::flush;

//	std::cout << "Transformed " << PRGcontigAlignment.graph_aligned_levels.size() << " alignment levels into " << graphSpace_graph_aligned_levels.size() << " PRG levels.\n" << std::flush;
	assert(graphSpace_graph_aligned_levels.front() != -1);
	assert(graphSpace_graph_aligned_levels.back() != -1);

	if(paranoid)
	{
//		std::cout << Utilities::join(Utilities::ItoStr(PRGcontigAlignment.graph_aligned_levels), ", ") << "\n";
//		std::cout << PRGcontigAlignment.graph_aligned << "\n";
//		std::cout << graphSpace_graph_aligned << "\n";
//		std::cout << "firstColumn" << ": " << firstColumn << "\n";
//		std::cout << "lastColumn" << ": " << lastColumn << "\n";
//		std::cout << "PRGcontigAlignment.graph_aligned_levels.size()" << ": " << PRGcontigAlignment.graph_aligned_levels.size() << "\n";
//		std::cout << "\n" << std::flush;

		if((firstColumn == 0) && (lastColumn == (PRGcontigAlignment.graph_aligned_levels.size() - 1)) && (alignmentLength_before_restriction == alignmentLength_after_restriction))
		{
			assert(Utilities::removeGaps(graphSpace_graph_aligned) == Utilities::removeGaps(PRGcontigAlignment.graph_aligned));
			assert(Utilities::removeGaps(graphSpace_sequence_aligned) == Utilities::removeGaps(PRGcontigAlignment.sequence_aligned));
		}
		int lastLevel;
		for(unsigned int cI = 0; cI < graphSpace_graph_aligned_levels.size(); cI++)
		{
			if(cI == 0)
			{
				assert(graphSpace_graph_aligned_levels.at(cI) != -1);
				lastLevel = graphSpace_graph_aligned_levels.at(cI);
			}
			else
			{
				if(graphSpace_graph_aligned_levels.at(cI) != -1)
				{
					assert(graphSpace_graph_aligned_levels.at(cI) == (lastLevel + 1));
					lastLevel = graphSpace_graph_aligned_levels.at(cI);
				}
			}
		}
	}

	/*
	std::cout << "Graph space alignment:\n" << std::flush;
	std::cout << "\t" << Utilities::join(Utilities::ItoStr(graphSpace_graph_aligned_levels), ", ") << "\n";
	std::cout << "\t" << graphSpace_graph_aligned << "\n";
	std::cout << "\t" << graphSpace_sequence_aligned << "\n";
	std::cout << "\n" << std::flush;
	*/

	class _seedChain_bt {
	public:
		double S;
		std::set<Edge*> takenEdges;
	};
	std::vector<std::map<Node*, _seedChain_bt>> seedChain_backtrack_graph;
	std::vector<std::map<Node*, _seedChain_bt>> seedChain_backtrack_sequence;

	seedChain_backtrack_graph.resize(graphSpace_graph_aligned_levels.size()+1);
	seedChain_backtrack_sequence.resize(graphSpace_graph_aligned_levels.size()+1);
	assert(graphSpace_graph_aligned_levels.size() > 0);
	if(!(graphSpace_graph_aligned_levels.at(0) < (int)g->NodesPerLevel.size()))
	{
		std::cerr << "graphSpace_graph_aligned_levels.at(0)" << ": " << graphSpace_graph_aligned_levels.at(0) << "\n";
		std::cerr << "g->NodesPerLevel.size()" << ": " << g->NodesPerLevel.size() << "\n";
		std::cerr << "\n" << std::flush;
	}
	assert(graphSpace_graph_aligned_levels.at(0) < (int)g->NodesPerLevel.size());
	const std::set<Node*>& firstNodes = g->NodesPerLevel.at(graphSpace_graph_aligned_levels.at(0));
	assert(firstNodes.size() > 0);
	for(std::set<Node*>::const_iterator nodeIt = firstNodes.begin(); nodeIt != firstNodes.end(); nodeIt++)
	{
		Node* n = *nodeIt;
		seedChain_backtrack_graph.at(0)[n].S = 0;
		seedChain_backtrack_sequence.at(0)[n].S = 0;
	}

	assert(graphSpace_graph_aligned_levels.at(0) != -1);

//	std::cout << "Path finding:\n" << std::flush;
	unsigned int lastColumnI_nonGap = 0;
	for(unsigned int columnI = 1; columnI <= graphSpace_graph_aligned_levels.size(); columnI++)
	{
//		std::cout << "columnI" << ": " << columnI << " / " <<  (graphSpace_graph_aligned_levels.size()-1) << "\n" << std::flush;
		if(graphSpace_graph_aligned_levels.at(columnI-1) == -1)
		{
			// this is a PRG gap - we will ignore these columns
			continue;
		}
		else
		{
//			std::cout << "lastColumnI_nonGap" << ": " << lastColumnI_nonGap << "\n" << std::flush;

			int graphLevel = graphSpace_graph_aligned_levels.at(columnI-1);
			char sequenceCharacter = graphSpace_sequence_aligned.at(columnI-1);
			char graphCharacter = graphSpace_graph_aligned.at(columnI-1);
			bool seedAlignment_isMatch = (sequenceCharacter == graphCharacter);

			assert(graphLevel != -1);

			std::map<Node*, _seedChain_bt>& lastRealizedColumn_graph = seedChain_backtrack_graph.at(lastColumnI_nonGap);
			std::map<Node*, _seedChain_bt>& lastRealizedColumn_sequence = seedChain_backtrack_sequence.at(lastColumnI_nonGap);
			assert(lastRealizedColumn_graph.size() > 0);
			assert(lastRealizedColumn_sequence.size() > 0);

			for(std::map<Node*, _seedChain_bt>::iterator nodeIt = lastRealizedColumn_graph.begin(); nodeIt != lastRealizedColumn_graph.end(); nodeIt++)
			{
				Node* fromN = nodeIt->first;
				if(! ((int)fromN->level == (graphLevel)))
				{
					std::cout << "columnI" << ": " << columnI << "\n";
					std::cout << "fromN->level" << ": " << fromN->level << "\n";
					std::cout << "graphLevel" << ": " << graphLevel << "\n";
					std::cout << "\n" << std::flush;
				}
				assert((int)fromN->level == (graphLevel));
				for(std::set<Edge*>::iterator edgeIt = fromN->Outgoing_Edges.begin(); edgeIt != fromN->Outgoing_Edges.end(); edgeIt++)
				{
					Edge* e = *edgeIt; 
					if(e->emission == graphCharacter)
					{
						Node* toN = e->To;
						if(seedChain_backtrack_graph.at(columnI).count(toN) == 0)
						{
							seedChain_backtrack_graph.at(columnI)[toN].S = seedChain_backtrack_graph.at(lastColumnI_nonGap).at(fromN).S + 1;
							seedChain_backtrack_graph.at(columnI)[toN].takenEdges.insert(e);
						}
						else
						{
							assert(seedChain_backtrack_graph.at(columnI).at(toN).S == (seedChain_backtrack_graph.at(lastColumnI_nonGap).at(fromN).S + 1));
							seedChain_backtrack_graph.at(columnI)[toN].takenEdges.insert(e);
						}
					}
				}
			}

			if(!(seedChain_backtrack_graph.at(columnI).size() > 0))
			{
				std::cerr << "Cannot find edge with label " << graphCharacter << " in column " << columnI << " graph level " << graphLevel << " (sequence character " << sequenceCharacter << ")\n" << std::flush;
				
				std::cerr << "Graph space alignment:\n" << std::flush;
				std::cerr << "\t" << Utilities::join(Utilities::ItoStr(graphSpace_graph_aligned_levels), ", ") << "\n";
				std::cerr << "\t" << graphSpace_graph_aligned << "\n";
				std::cerr << "\t" << graphSpace_sequence_aligned << "\n";
				std::cerr << "\n" << std::flush;
	
				std::cerr << "Edges: " << "\n";
				std::set<Node*> nodes_at_level = g->NodesPerLevel.at(graphLevel);
				std::set<Edge*> edges_at_level;
				for(std::set<Node*>::iterator nIt = nodes_at_level.begin(); nIt != nodes_at_level.end(); nIt++)
				{
					Node* n = *nIt;
					edges_at_level.insert(n->Outgoing_Edges.begin(), n->Outgoing_Edges.end());
				}
				for(std::set<Edge*>::iterator eIt = edges_at_level.begin(); eIt != edges_at_level.end(); eIt++)
				{
					Edge* e = *eIt;
					std::cerr << "\t" << e << "\t" << e->getEmission() << "\n";
				}
				std::cerr << std::flush;
			}
			assert(seedChain_backtrack_graph.at(columnI).size() > 0);

			for(std::map<Node*, _seedChain_bt>::iterator nodeIt = lastRealizedColumn_sequence.begin(); nodeIt != lastRealizedColumn_sequence.end(); nodeIt++)
			{
				Node* fromN = nodeIt->first;
				if(! ((int)fromN->level == (graphLevel)))
				{
					std::cout << "columnI" << ": " << columnI << "\n";
					std::cout << "fromN->level" << ": " << fromN->level << "\n";
					std::cout << "graphLevel" << ": " << graphLevel << "\n";
					std::cout << "\n" << std::flush;
				}
				assert((int)fromN->level == (graphLevel));
				for(std::set<Edge*>::iterator edgeIt = fromN->Outgoing_Edges.begin(); edgeIt != fromN->Outgoing_Edges.end(); edgeIt++)
				{
					Edge* e = *edgeIt;
					if(seedAlignment_isMatch)
					{
						if(e->emission != sequenceCharacter)
						{
							continue;
						}
					}

					double S = (e->emission == sequenceCharacter) ? 1 : 0;
					Node* toN = e->To;
					if(seedChain_backtrack_sequence.at(columnI).count(toN) == 0)
					{
						seedChain_backtrack_sequence.at(columnI)[toN].S = seedChain_backtrack_sequence.at(lastColumnI_nonGap).at(fromN).S + S;
						seedChain_backtrack_sequence.at(columnI)[toN].takenEdges.insert(e);
					}
					else
					{
						if(seedChain_backtrack_sequence.at(columnI).at(toN).S == (seedChain_backtrack_sequence.at(lastColumnI_nonGap).at(fromN).S + S))
						{
							seedChain_backtrack_sequence.at(columnI)[toN].takenEdges.insert(e);
						}
						else if (seedChain_backtrack_sequence.at(columnI).at(toN).S < (seedChain_backtrack_sequence.at(lastColumnI_nonGap).at(fromN).S + S))
						{
							seedChain_backtrack_sequence.at(columnI)[toN].takenEdges.clear();
							seedChain_backtrack_sequence.at(columnI)[toN].takenEdges.insert(e);
							seedChain_backtrack_sequence.at(columnI)[toN].S = seedChain_backtrack_sequence.at(lastColumnI_nonGap).at(fromN).S + S;
						}
					}
				}
			}
			lastColumnI_nonGap = columnI;
		}
	}


	auto backtrace = [&](std::vector<std::map<Node*, _seedChain_bt>>& bt) -> reads::verboseSeedChain {
		std::vector<int> bt_graph_aligned_levels;
		std::vector<Edge*> bt_graph_aligned_edges;
		std::string bt_graph_aligned;
		std::string bt_sequence_aligned;

		std::set<Edge*> startEdges_left;
		std::set<Edge*> startEdges_right;

		std::vector<double> S;
		std::vector<Node*> n;
		assert(graphSpace_graph_aligned_levels.at(graphSpace_graph_aligned_levels.size()-1) != -1);
		assert(bt.at(graphSpace_graph_aligned_levels.size()).size() > 0);
		for(std::map<Node*, _seedChain_bt>::iterator nIt = bt.at(graphSpace_graph_aligned_levels.size()).begin(); nIt != bt.at(graphSpace_graph_aligned_levels.size()).end(); nIt++)
		{
			S.push_back(nIt->second.S);
			n.push_back(nIt->first);
		}
		double Smax = Utilities::findVectorMax(S).first;
		std::set<Node*> runningN;
		Node* runningN_singular;

		for(unsigned int i = 0; i < S.size(); i++)
		{
			if(S.at(i) == Smax)
			{
				runningN.insert(n.at(i));
			}
		}
		runningN_singular = *(runningN.begin());

//		std::cout << "Backtrace" << "\n";
		for(int columnI = graphSpace_graph_aligned_levels.size(); columnI >= 1; columnI--)
		{
//			std::cout << "columnI" << ": " << columnI << "\n" << std::flush;

			if(graphSpace_graph_aligned_levels.at(columnI-1) == -1)
			{
				bt_graph_aligned_levels.push_back(-1);
				bt_graph_aligned_edges.push_back(0);
				bt_graph_aligned.push_back('_');
				bt_sequence_aligned.push_back(graphSpace_sequence_aligned.at(columnI-1));

				// this is a PRG gap - we will ignore these columns
				continue;
			}
			else
			{
				std::set<Node*> nextLevel_runningN;
				Node* nextLevel_runningN_singular = 0;

				for(std::set<Node*>::iterator nIt = runningN.begin(); nIt != runningN.end(); nIt++)
				{
					_seedChain_bt& thisN_bt = bt.at(columnI).at(*nIt);
					for(std::set<Edge*>::iterator eIt = thisN_bt.takenEdges.begin(); eIt != thisN_bt.takenEdges.end(); eIt++)
					{
						Edge* e = *eIt;
						nextLevel_runningN.insert(e->From);

						if(columnI == ((int)graphSpace_graph_aligned_levels.size() - 1))
						{
							startEdges_right.insert(e);
						}
						if(columnI == 1)
						{
							startEdges_left.insert(e);
						}

						if(((*nIt) == runningN_singular) && (nextLevel_runningN_singular == 0))
						{
							nextLevel_runningN_singular = e->From;

							assert((int)e->From->level == graphSpace_graph_aligned_levels.at(columnI-1));
							bt_graph_aligned_levels.push_back(graphSpace_graph_aligned_levels.at(columnI-1));
							bt_graph_aligned_edges.push_back(e);
							bt_graph_aligned.push_back(e->emission);
							bt_sequence_aligned.push_back(graphSpace_sequence_aligned.at(columnI-1));
						}
					}
				}

				assert(nextLevel_runningN_singular != 0);
				if(columnI == 1)
				{
					for(std::set<Node*>::iterator nIt = nextLevel_runningN.begin(); nIt != nextLevel_runningN.end(); nIt++)
					{
						Node* n = *nIt;
						assert((int)n->level == graphSpace_graph_aligned_levels.front());
					}
				}

				runningN = nextLevel_runningN;
				runningN_singular = nextLevel_runningN_singular;
			}
		}

		if((!((bt_graph_aligned_levels.size() == graphSpace_graph_aligned_levels.size()))))
		{
			std::cout << "bt_graph_aligned_levels.size()" << ": " << bt_graph_aligned_levels.size() << "\n";
			std::cout << "graphSpace_graph_aligned_levels.size()" << ": " << graphSpace_graph_aligned_levels.size() << "\n";
			std::cout << "\n" << std::flush;
		}

		assert(bt_graph_aligned_levels.size() == graphSpace_graph_aligned_levels.size());
		assert(bt_graph_aligned_edges.size() == graphSpace_graph_aligned_levels.size());
		assert(bt_graph_aligned.size() == graphSpace_graph_aligned_levels.size());
		assert(bt_sequence_aligned.size() == graphSpace_graph_aligned_levels.size());

		std::reverse(bt_graph_aligned_levels.begin(), bt_graph_aligned_levels.end());
		std::reverse(bt_graph_aligned_edges.begin(), bt_graph_aligned_edges.end());
		std::reverse(bt_graph_aligned.begin(), bt_graph_aligned.end());
		std::reverse(bt_sequence_aligned.begin(), bt_sequence_aligned.end());

		reads::verboseSeedChain forReturn;
		forReturn.graph_aligned_edges = bt_graph_aligned_edges;
		forReturn.graph_aligned_levels = bt_graph_aligned_levels;
		forReturn.graph_aligned = bt_graph_aligned;
		forReturn.sequence_aligned = bt_sequence_aligned;
		forReturn.sequence_begin = graphSpace_sequence_aligned_startInRaw;
		forReturn.sequence_end = graphSpace_sequence_aligned_stopInRaw;
		forReturn.reverse = PRGcontigAlignment.reverse;
		forReturn.removed_columns_noGap_restriction = removed_columns_noGap_restriction;

		size_t after_backtrace_matches = 0;
		size_t after_backtrace_length = 0;
		for(unsigned int posI = 0; posI < bt_graph_aligned.size(); posI++)
		{
			if(bt_graph_aligned.at(posI) == bt_sequence_aligned.at(posI))
			{
				after_backtrace_matches++;
			}
			after_backtrace_length++;
		}

		double pre_backtrace_prop_matches = (double)before_backtrace_matches / (double)before_backtrace_length;
		double after_backtrace_prop_matches = (double)after_backtrace_matches / (double)after_backtrace_length;
		double bt_improvement = after_backtrace_prop_matches - pre_backtrace_prop_matches;
		forReturn.improvement_through_bt = bt_improvement;

		if(paranoid)
		{
			if((firstColumn == 0) && (lastColumn == (PRGcontigAlignment.graph_aligned_levels.size()-1))&& (alignmentLength_before_restriction == alignmentLength_after_restriction))
			{
				std::string bt_sequence_aligned_noGaps = Utilities::removeGaps(bt_sequence_aligned);
				std::string PRGcontigAlignment_sequence_aligned_noGaps = Utilities::removeGaps(PRGcontigAlignment.sequence_aligned);
				if(!(bt_sequence_aligned_noGaps == PRGcontigAlignment_sequence_aligned_noGaps))
				{
					std::cerr << "! (bt_sequence_aligned_noGaps == PRGcontigAlignment_sequence_aligned_noGaps)" << "\n";

					std::cerr << "\t" << "PRGcontigAlignment.graph_aligned_levels" << ": " << Utilities::join(Utilities::ItoStr(PRGcontigAlignment.graph_aligned_levels), ", ") << "\n\n";
					std::cerr << "\t" << "PRGcontigAlignment.sequence_aligned" << ": " << PRGcontigAlignment.sequence_aligned << "\n\n";
					std::cerr << "\t" << "bt_sequence_aligned" << ": " << bt_sequence_aligned << "\n";

					std::cerr << "\t" << "PRGcontigAlignment_sequence_aligned_noGaps" << ": " << PRGcontigAlignment_sequence_aligned_noGaps << "\n";
					std::cerr << "\t" << "bt_sequence_aligned_noGaps" << ": " << bt_sequence_aligned_noGaps << "\n";

					std::cerr << "\t" << "graphSpace_sequence_aligned_startInRaw" << ": "  << graphSpace_sequence_aligned_startInRaw << "\n";
					std::cerr << "\t" << "graphSpace_sequence_aligned_stopInRaw" << ": "  << graphSpace_sequence_aligned_stopInRaw << "\n";
					std::cerr << "\t" << "firstColumn" << ": "  << firstColumn << "\n";
					std::cerr << "\t" << "lastColumn" << ": "  << lastColumn << "\n";
					std::cerr << "\t" << "Input alignments columns (size): " << ": "  << PRGcontigAlignment.graph_aligned_levels.size() << "\n";

					std::cerr << "\n" << std::flush;
				}

				assert(bt_sequence_aligned_noGaps == PRGcontigAlignment_sequence_aligned_noGaps);
			}
		}
		return forReturn;
	};

//	std::cout << "start backtrace\n" << std::flush;

	graphSeed = backtrace(seedChain_backtrack_graph);
	sequenceSeed = backtrace(seedChain_backtrack_sequence);



	return false;
}

reads::verboseSeedChain processBAM::alignment2Chain(const std::tuple<std::string, int,  BamTools::BamAlignment, int>& al, const std::string& complete_read_sequence, const std::string& complete_read_qualities, bool verbose) const
{
	//const std::string& BAMid = (_currentBAM != "") ? R.GetReferenceData().at(al.second.RefID).RefName : R1.GetReferenceData().at(al.second.RefID).RefName;

	const std::string& BAMid = std::get<0>(al);

	if(_currentBAM == "")
	{
		assert(_currentBAM1 != "");
		assert(R2.GetReferenceData().at(std::get<2>(al).RefID).RefName == BAMid);
	}
	
	if(!(BAMid_2_PRGid.count(BAMid)))
	{
		std::cerr << "! (BAMid_2_PRGid.count(BAMid))" << "\n";
		std::cerr << "\t" << "BAMid" << ": " << BAMid << "\n";
		std::cerr << std::flush;
	}
	assert(BAMid_2_PRGid.count(BAMid));
	
	int PRGid = BAMid_2_PRGid.at(BAMid);

	if(!(extendedReferenceGenome_levelTranslation.count(PRGid)))
	{
		std::cerr << "! extendedReferenceGenome_levelTranslation.count(PRGid)" << "\n";
		std::cerr << "\t" << "BAMid" << ": " << BAMid << "\n";
		std::cerr << "\t" << "PRGid" << ": " << PRGid << "\n";
		std::cerr << std::flush;
	}
	assert(extendedReferenceGenome_levelTranslation.count(PRGid));
	
	reads::PRGContigBAMAlignment PRGcontigAlignment;
	transformBAMreadToInternalAlignment(
			(_currentBAM_isExtendedReferenceGenome ? extendedReferenceGenomeSequences.at(BAMid) : PRGonlyReferenceGenomeSequences.at(BAMid)),
			extendedReferenceGenome_levelTranslation.at(PRGid),
			std::get<1>(al),
			std::get<2>(al),
			complete_read_sequence,
			complete_read_qualities,
			PRGcontigAlignment
	);

	{
		assert(complete_read_sequence.find('_') == std::string::npos);
		std::string subSequence = complete_read_sequence.substr(PRGcontigAlignment.sequence_aligned_startInRaw, PRGcontigAlignment.sequence_aligned_stopInRaw - PRGcontigAlignment.sequence_aligned_startInRaw + 1);


		std::string sequence_aligned_noGaps = filter(PRGcontigAlignment.sequence_aligned, [](unsigned char c){return (c != '_');});

		if(!(sequence_aligned_noGaps == subSequence))
		{
			std::cerr << "Mismatch in sequence check!\n";
			std::cerr << "\t" << "complete_read_sequence" << ": " << complete_read_sequence << "\n";
			std::cerr << "\t" << "PRGcontigAlignment.sequence_aligned" << ": " << PRGcontigAlignment.sequence_aligned << "\n";
			std::cerr << "\t" << "PRGcontigAlignment.graph_aligned" << ": " << PRGcontigAlignment.graph_aligned << "\n";
			std::cerr << "\t" << "sequence_aligned_noGaps" << ": " << sequence_aligned_noGaps << "\n";
			std::cerr << "\t" << "subSequence" << ": " << subSequence << "\n";
			std::cerr << "\t" << "PRGcontigAlignment.sequence_aligned_startInRaw" << ": " << PRGcontigAlignment.sequence_aligned_startInRaw << "\n";
			std::cerr << "\t" << "PRGcontigAlignment.sequence_aligned_stopInRaw" << ": " << PRGcontigAlignment.sequence_aligned_stopInRaw << "\n";
			std::cerr << "\n" << std::flush;
			
			reads::protoSeeds::printAlignmentInfo(1, std::get<2>(al));
			
			assert(1 == 0);
		}		
		
	}
	
	if(verbose)
	{
		std::cout << std::flush;
		std::cerr << std::flush;

		std::string subSequence = complete_read_sequence.substr(PRGcontigAlignment.sequence_aligned_startInRaw, PRGcontigAlignment.sequence_aligned_stopInRaw - PRGcontigAlignment.sequence_aligned_startInRaw + 1);
		std::string sequence_aligned_noGaps = filter(PRGcontigAlignment.sequence_aligned, [](unsigned char c){return (c != '_');});

		std::cerr << "processBAM::alignment2Chain(..): Verbose mode!\n";
		std::cerr << "\t" << "BAM Reference ID" << ": " << BAMid << "\n";
		std::cerr << "\t" << "complete_read_sequence" << ": " << complete_read_sequence << "\n";
		std::cerr << "\t" << "PRGcontigAlignment.sequence_aligned" << ": " << PRGcontigAlignment.sequence_aligned << "\n";
		std::cerr << "\t" << "sequence_aligned_noGaps" << ": " << sequence_aligned_noGaps << "\n";
		std::cerr << "\t" << "subSequence" << ": " << subSequence << "\n";
		std::cerr << "\t" << "PRGcontigAlignment.sequence_aligned_startInRaw" << ": " << PRGcontigAlignment.sequence_aligned_startInRaw << "\n";
		std::cerr << "\t" << "PRGcontigAlignment.sequence_aligned_stopInRaw" << ": " << PRGcontigAlignment.sequence_aligned_stopInRaw << "\n";
		std::cerr << "\n" << std::flush;

		reads::protoSeeds::printAlignmentInfo(1, std::get<2>(al));

		std::cout << std::flush;
		std::cerr << std::flush;
	}

	assert(paranoid);
	if(paranoid)
	{
		PRGcontigAlignment.checkAlignmentConcordanceWithSequence(complete_read_sequence);
	}

	reads::verboseSeedChain forReturn_graphSeed;
	reads::verboseSeedChain forReturn_sequenceSeed;

	// std::cout << "BAM Reference ID" << ": " << BAMid << "\n" << std::flush;
	PRGContigAlignment2Seed(g, PRGcontigAlignment, paranoid, forReturn_graphSeed, forReturn_sequenceSeed, inGraphGapStretch);

	forReturn_graphSeed.is_from_BWAseed.resize(forReturn_graphSeed.graph_aligned_levels.size(), 1);
	forReturn_sequenceSeed.is_from_BWAseed.resize(forReturn_sequenceSeed.graph_aligned_levels.size(), 1);

	return forReturn_sequenceSeed;
}

reads::verboseSeedChainPair processBAM::alignOneReadPair(const reads::protoSeeds& protoSeed, const boost::math::normal& rnd_InsertSize, double max_insertsize_penalty_log, simulator::trueReadLevels* trueReadLevels, aligner::statistics* statisticsStore) const
{
	// bool verbose = ((protoSeed.read1_alignments.size() > 1) || (protoSeed.read2_alignments.size() > 1));
	bool verbose = false;

	std::set<std::string> printReadIDs;
	printReadIDs.insert("H781JADXX131217:1:2115:10484:82360");
	printReadIDs.insert("H781JADXX131217:2:1205:17348:100758");
	printReadIDs.insert("H781JADXX131217:2:1216:9390:54498");

	size_t r1_primary = protoSeed.read1_getPrimaryAlignmentI();
	size_t r2_primary = protoSeed.read2_getPrimaryAlignmentI();

	std::string r1_QueryBases_alignmentOrientation = std::get<2>(protoSeed.read1_alignments.at(r1_primary)).QueryBases;
	std::string r2_QueryBases_alignmentOrientation = std::get<2>(protoSeed.read2_alignments.at(r2_primary)).QueryBases;
	std::string r1_Qualities_alignmentOrientation = std::get<2>(protoSeed.read1_alignments.at(r1_primary)).Qualities;
	std::string r2_Qualities_alignmentOrientation = std::get<2>(protoSeed.read2_alignments.at(r2_primary)).Qualities;
	assert(r1_QueryBases_alignmentOrientation.length() == r1_Qualities_alignmentOrientation.length());
	assert(r2_QueryBases_alignmentOrientation.length() == r2_Qualities_alignmentOrientation.length());

	reads::oneRead r1(
		std::get<2>(protoSeed.read1_alignments.at(r1_primary)).Name,
		r1_QueryBases_alignmentOrientation,
		std::get<2>(protoSeed.read1_alignments.at(r1_primary)).Qualities
	);
	if(std::get<2>(protoSeed.read1_alignments.at(r1_primary)).IsReverseStrand())
	{
		r1.invert();
	}
	reads::oneRead r2(
		std::get<2>(protoSeed.read2_alignments.at(r2_primary)).Name,
		r2_QueryBases_alignmentOrientation,
		std::get<2>(protoSeed.read2_alignments.at(r2_primary)).Qualities
	);
	if(std::get<2>(protoSeed.read2_alignments.at(r2_primary)).IsReverseStrand())
	{
		r2.invert();
	}
	reads::oneReadPair rP(r1, r2, 0);

	if(printReadIDs.count(std::get<2>(protoSeed.read1_alignments.at(r1_primary)).Name) || printReadIDs.count(std::get<2>(protoSeed.read2_alignments.at(r2_primary)).Name))
	{
		verbose = true;
	}

	if(verbose)
	{
		std::cout << "Debug output processBAM::alignOneReadPair(..) for " << std::get<2>(protoSeed.read1_alignments.at(r1_primary)).Name << ":\n";
		std::cout << "==============================================================\n";
		std::cout << "Read 1: " <<  r1.sequence << "\n";
		std::cout << "Read 2: " <<  r2.sequence << "\n";
	}

	// now process chains
	if(verbose)
	{
		std::cout << "Read 1 chains:\n";
	}

	bool skipIdenticalCoordinates = true;

	std::vector<reads::verboseSeedChain> read1_extendedChains;
	std::vector<double> read1_extendedChains_log_likelihoods;
	// std::cout << "Read chains " << protoSeed.read1_alignments.size() << "\n";

	size_t considered_chains = 0;
	size_t n_calledChainExtension = 0;
	size_t n_calledChainExtension_primary = 0;
	size_t considered_chain_pairs = 0;

	std::map<std::string, int> read1_alignments_scores;
	for(size_t chainI = 0; chainI < protoSeed.read1_alignments.size(); chainI++)
	{
		considered_chains++;

		const std::tuple<std::string, int,  BamTools::BamAlignment, int>& al = protoSeed.read1_alignments.at(chainI);
		std::pair<int, int> al_startStop_PRG = alignment_get_startstop_PRGcoordinates(al);

		std::string id = Utilities::ItoStr(al_startStop_PRG.first) + "//" + Utilities::ItoStr(al_startStop_PRG.second);
		int score = getAlignmentScore(std::get<2>(al));

		if(verbose)
		{
			std::cout << "\tRead 1 chain " << chainI << " / " << protoSeed.read1_alignments.size() << " id: " << id << "\n";
			std::cout << "\t\tScore: " << score << "\n" << std::flush;						
		}		
		
		if(std::get<2>(al).IsReverseStrand() != std::get<2>(protoSeed.read1_alignments.at(r1_primary)).IsReverseStrand())
		{
			if(verbose)
			{
				std::cout << "\t\tSkip reverse!" << "\n" << std::flush;		
			}			
			continue;
		}
		
		reads::verboseSeedChain al_Chain = alignment2Chain(al, r1_QueryBases_alignmentOrientation, r1_Qualities_alignmentOrientation);
		
		if(verbose)
		{
			std::set<std::string> segments = al_Chain.getSegments(g_level_names);
			std::vector<std::string> segments_vector(segments.begin(), segments.end());
			std::cout << "\t\tSegments: " << Utilities::join(segments_vector, ", ") << "\n" << std::flush;
		}
		
		if(skipIdenticalCoordinates && read1_alignments_scores.count(id) && (read1_alignments_scores.at(id) >= score))
		{
			if(verbose)
			{
				std::cout << "\t\tSkip identical coordinates!" << "\n" << std::flush;		
			}					
			continue;
		}



		// std::cout << "\t" << al_startStop_PRG.first << " " << al_startStop_PRG.second << " " << score << "\n" << std::flush;

		assert(al_Chain.is_from_BWAseed.size() == al_Chain.graph_aligned_levels.size());

		if(verbose)
		{
			std::cout << "\tRead 1 chain " << chainI << " / " << protoSeed.read1_alignments.size() << "\n";
			al_Chain.print(g_level_names);
		}
 
		if(paranoid)
			al_Chain.checkChainConcordanceWithSequence(r1_QueryBases_alignmentOrientation);

		n_calledChainExtension++;
		if(chainI == r1_primary)
			n_calledChainExtension_primary++;

		reads::verboseSeedChain al_Chain_extended = eA->extendSeedChain(r1_QueryBases_alignmentOrientation, al_Chain);
		
		if(verbose)
		{
			std::cout << "\tRead 1 EXTENDED chain " << chainI << " / " << protoSeed.read1_alignments.size() << "\n";
			al_Chain_extended.print(g_level_names);
		}		
		
		assert(al_Chain_extended.is_from_BWAseed.size() == al_Chain_extended.graph_aligned_levels.size());

		read1_extendedChains.push_back(al_Chain_extended);
		read1_extendedChains_log_likelihoods.push_back(eA->scoreOneAlignment(al_Chain_extended, r1));

		if((read1_alignments_scores.count(id) == 0) || (read1_alignments_scores.at(id) < score))
		{
			read1_alignments_scores[id] = score;
		}
	}

	if(verbose)
	{
		std::cout << "Read 2 chains:\n";
	}

	std::vector<reads::verboseSeedChain> read2_extendedChains;
	std::vector<double> read2_extendedChains_log_likelihoods;
	std::map<std::string, int> read2_alignments_scores;
	for(size_t chainI = 0; chainI < protoSeed.read2_alignments.size(); chainI++)
	{
		considered_chains++;

		const std::tuple<std::string, int,  BamTools::BamAlignment, int>& al = protoSeed.read2_alignments.at(chainI);
		std::pair<int, int> al_startStop_PRG = alignment_get_startstop_PRGcoordinates(al);
		
		std::string id = Utilities::ItoStr(al_startStop_PRG.first) + "//" + Utilities::ItoStr(al_startStop_PRG.second);
		int score = getAlignmentScore(std::get<2>(al));
		
		if(verbose)
		{
			std::cout << "\tRead 2 chain " << chainI << " / " << protoSeed.read1_alignments.size() << " id: " << id << "\n";		
			std::cout << "\t\tScore: " << score << "\n" << std::flush;
		}
		
		if(std::get<2>(al).IsReverseStrand() != std::get<2>(protoSeed.read2_alignments.at(r2_primary)).IsReverseStrand())
		{
			if(verbose)
			{
				std::cout << "\t\tSkip reverse!" << "\n" << std::flush;				
			}
			continue;
		}		
		
	
		reads::verboseSeedChain al_Chain = alignment2Chain(al, r2_QueryBases_alignmentOrientation, r2_Qualities_alignmentOrientation);
		
		if(verbose)
		{
			std::set<std::string> segments = al_Chain.getSegments(g_level_names);
			std::vector<std::string> segments_vector(segments.begin(), segments.end());
			std::cout << "\t\tSegments: " << Utilities::join(segments_vector, ", ") << "\n" << std::flush;
		}
				
		if(skipIdenticalCoordinates && read2_alignments_scores.count(id) && (read2_alignments_scores.at(id) >= score))
		{
			if(verbose)
			{
				std::cout << "\t\tSkip identical coordinates!" << "\n" << std::flush;		
			}								
			continue;
		}



		assert(al_Chain.is_from_BWAseed.size() == al_Chain.graph_aligned_levels.size());

		if(verbose)
		{
			std::cout << "\tRead 2 chain " << chainI << " / " << protoSeed.read2_alignments.size() << "\n";
			al_Chain.print(g_level_names);
		}

		if(paranoid)
			al_Chain.checkChainConcordanceWithSequence(r2_QueryBases_alignmentOrientation);

		n_calledChainExtension++;
		if(chainI == r2_primary)
			n_calledChainExtension_primary++;
		reads::verboseSeedChain al_Chain_extended = eA->extendSeedChain(r2_QueryBases_alignmentOrientation, al_Chain);

		assert(al_Chain_extended.is_from_BWAseed.size() == al_Chain_extended.graph_aligned_levels.size());

		read2_extendedChains.push_back(al_Chain_extended);
		read2_extendedChains_log_likelihoods.push_back(eA->scoreOneAlignment(al_Chain_extended, r2));

		if((read2_alignments_scores.count(id) == 0) || (read2_alignments_scores.at(id) < score))
		{
			read2_alignments_scores[id] = score;
		}
	}

	if(verbose)
	{
		std::vector<double> ll_r1 = read1_extendedChains_log_likelihoods;
		std::vector<double> ll_r2 = read2_extendedChains_log_likelihoods;
		ll_r1 = Utilities::normalize_log_vector(ll_r1);
		ll_r2 = Utilities::normalize_log_vector(ll_r2);

		std::cout << "R1 extended chains:\n";
		for(unsigned int cI = 0; cI < ll_r1.size(); cI++)
		{
			std::set<std::string> chain_segments = read1_extendedChains.at(cI).getSegments(g_level_names);
			std::cout << "\t" << cI << "\t"
					<< read1_extendedChains.at(cI).alignment_firstLevel() << " - " << read1_extendedChains.at(cI).alignment_lastLevel() << "\t"
					<< Utilities::join(std::vector<std::string>(chain_segments.begin(), chain_segments.end()), ", ") << "\t"
					<< ll_r1.at(cI)
					<< "\n";
		}

		std::cout << "R2 extended chains:\n";
		for(unsigned int cI = 0; cI < ll_r2.size(); cI++)
		{
			std::set<std::string> chain_segments = read2_extendedChains.at(cI).getSegments(g_level_names);
			std::cout << "\t" << cI << "\t"
					<< read2_extendedChains.at(cI).alignment_firstLevel() << " - " << read2_extendedChains.at(cI).alignment_lastLevel() << "\t"
					<< Utilities::join(std::vector<std::string>(chain_segments.begin(), chain_segments.end()), ", ") << "\t"
					<< ll_r2.at(cI)
					<< "\n";
		}

	}

	reads::verboseSeedChainPair forReturn;
	assert(read1_extendedChains.size() > 0);
	assert(read2_extendedChains.size() > 0);
	assert(std::get<2>(protoSeed.read1_alignments.at(r1_primary)).Name == std::get<2>(protoSeed.read2_alignments.at(r2_primary)).Name);
	forReturn.readID = std::get<2>(protoSeed.read1_alignments.at(r1_primary)).Name;
	std::vector<std::pair<unsigned int, unsigned int>> combinations_indices;
	std::vector<double> combinations_LL;
	std::vector<double> combinations_LL_insertSizeOnly;

	bool verboseReadPairing = verbose;
	if(verboseReadPairing)
	{
		std::cerr << "Read pairing process..\n";
	}

	for(unsigned int read1_alignmentI = 0; read1_alignmentI < read1_extendedChains.size(); read1_alignmentI++)
	{
		for(unsigned int read2_alignmentI = 0; read2_alignmentI < read2_extendedChains.size(); read2_alignmentI++)
		{
			considered_chain_pairs++;

			double combined_log_likelihood = read1_extendedChains_log_likelihoods.at(read1_alignmentI) +
												read2_extendedChains_log_likelihoods.at(read2_alignmentI);

			if(verboseReadPairing && 0)
			{
				std::cerr << "\tPairs " << read1_alignmentI << " / " << read2_alignmentI << "\n";
				std::cerr << "\t\tCombined likelihood " << combined_log_likelihood << "\n";
			}
			const reads::verboseSeedChain& chain_read1 = read1_extendedChains.at(read1_alignmentI);
			const reads::verboseSeedChain& chain_read2 = read2_extendedChains.at(read2_alignmentI);

			bool strandsValid = eA->alignedReadPair_strandsValid(chain_read1, chain_read2);

			if(verboseReadPairing && 0)
				std::cerr << "\t\t" << "strandsValid" << ": " << strandsValid << "\n";

			double log_likelihood_insertSize;
			if(strandsValid)
			{
				int graphDistance = eA->alignedReadPair_pairsDistanceInGraphLevels(chain_read1, chain_read2);
				std::set<int> underlyingSequencesDistances = eA->alignedReadPair_pairsDistancesUnderlyingSequences(chain_read1, chain_read2, graphLevel_2_underlyingSequencePositions);

				if(underlyingSequencesDistances.size())
				{
					std::vector<double> underlyingSequenceDistances_LLs;
					for(std::set<int>::iterator distanceIt = underlyingSequencesDistances.begin(); distanceIt != underlyingSequencesDistances.end(); distanceIt++)
					{
						int distance = *distanceIt;
						if(verboseReadPairing && 0)
						{
							std::cerr << "\t\t\tdistance: " << distance << "\n";
						}
						double distance_P = boost::math::pdf(rnd_InsertSize, distance);
						if(distance_P <= 0)
						{					
							/*
							std::cerr << "distance" << ": " << distance << "\n";
							std::cerr << "distance_P" << ": " << distance_P << "\n";
							std::cerr << "insertSize_mean" << ": " << rnd_InsertSize.mean() << "\n";
							std::cerr << "insertSize_sd" << ": " << rnd_InsertSize.standard_deviation() << "\n";
							std::cerr << "\n" << std::flush;
							assert(2 == 3);
							*/
							
							underlyingSequenceDistances_LLs.push_back(max_insertsize_penalty_log);
						}
						else
						{
							assert(distance_P > 0);
							assert(distance_P <= 1);
							underlyingSequenceDistances_LLs.push_back(log(distance_P));
						}
					}
					std::pair<double, unsigned int> best_underlyingSequenceDistances_LLs = Utilities::findVectorMax(underlyingSequenceDistances_LLs);
					log_likelihood_insertSize = best_underlyingSequenceDistances_LLs.first;
				}
				else
				{
					log_likelihood_insertSize = max_insertsize_penalty_log;
				}

//				std::set<int> underlyingSequencesDistances_properGraph = eA->alignedReadPair_properGraphDistance(chain_read1, chain_read2);
//				if(verboseReadPairing || (underlyingSequencesDistances.count(graphDistance) == 0))
//				{
//					std::cerr << "!!!!" << "\n";
//					std::cerr << "\t\t" << "graphDistance" << ": " << graphDistance << "\n";
//					std::cerr << "\t\t" << "Underlying sequence distances: " << "\n";
//					for(std::set<int>::iterator distanceIt = underlyingSequencesDistances.begin(); distanceIt != underlyingSequencesDistances.end(); distanceIt++)
//					{
//						std::cerr << "\t\t\t" << *distanceIt << "\n";
//					}
//					std::cerr << "\t\t" << "Proper graph distances: " << "\n";
//					for(std::set<int>::iterator distanceIt = underlyingSequencesDistances_properGraph.begin(); distanceIt != underlyingSequencesDistances_properGraph.end(); distanceIt++)
//					{
//						std::cerr << "\t\t\t" << *distanceIt << "\n";
//					}
//				}
			}
			else
			{
				log_likelihood_insertSize = max_insertsize_penalty_log;
			}

			combined_log_likelihood += log_likelihood_insertSize;

			combinations_LL.push_back(combined_log_likelihood);
			combinations_indices.push_back(make_pair(read1_alignmentI, read2_alignmentI));

			combinations_LL_insertSizeOnly.push_back(log_likelihood_insertSize);

			std::cerr << std::flush;
		}
	}

	if(verbose)
	{
		std::vector<double> combinations_LL_normalized = combinations_LL;
		combinations_LL_normalized = Utilities::normalize_log_vector(combinations_LL_normalized);

		std::cout << "Combinations:\n";

		unsigned int combinationI = 0;
		for(unsigned int read1_alignmentI = 0; read1_alignmentI < read1_extendedChains.size(); read1_alignmentI++)
		{
			for(unsigned int read2_alignmentI = 0; read2_alignmentI < read2_extendedChains.size(); read2_alignmentI++)
			{
				const reads::verboseSeedChain& chain_read1 = read1_extendedChains.at(read1_alignmentI);
				const reads::verboseSeedChain& chain_read2 = read2_extendedChains.at(read2_alignmentI);
				bool strandsValid = eA->alignedReadPair_strandsValid(chain_read1, chain_read2);

				// assert(combinations_LL_insertSizeOnly.at(combinationI) > combinations_LL.at(combinationI));
				// double prop_insertSize = exp(combinations_LL_insertSizeOnly.at(combinationI) - combinations_LL.at(combinationI));

				std::cout << "\t" << combinationI << " " << read1_alignmentI << "/" << read2_alignmentI
							<< "\t" << "strandsValid: " << strandsValid
							<< "\t" << "P: " << combinations_LL_normalized.at(combinationI)
							<< "\t" << "IS score (log): " << combinations_LL_insertSizeOnly.at(combinationI)
							<< "\n";

				combinationI++;
			}
		}
	}

	std::pair<double, unsigned int> combinations_max = Utilities::findVectorMax(combinations_LL);
	unsigned int maxCombination_i1 = combinations_indices.at(combinations_max.second).first;
	unsigned int maxCombination_i2 = combinations_indices.at(combinations_max.second).second;

	forReturn.chains.first = read1_extendedChains.at(maxCombination_i1);
	forReturn.chains.second = read2_extendedChains.at(maxCombination_i2);

	assignMappingQualities(forReturn, combinations_indices, combinations_LL, combinations_max, read1_extendedChains, read2_extendedChains);
	
	assert(forReturn.chains.first.mapQ_perPosition.size() == forReturn.chains.first.sequence_aligned.size());
	assert(forReturn.chains.second.mapQ_perPosition.size() == forReturn.chains.second.sequence_aligned.size());
	
	if(trueReadLevels != 0)
	{
		#pragma omp critical
		{
			trueReadLevels->evaluateAlignment(forReturn, rP, eA);
		}
	}

	if(verbose)
	{
		std::cout << "Resulting alignment:\n";
		forReturn.print(g_level_names);
		std::cout << "end debug\n\n";
	}


	if(statisticsStore != 0)
	{
		statisticsStore->n_calls_alignOneReadPair++;
		
		size_t selectedChains_totalColumns = 0;
		size_t selectedChains_totalColumns_fromBWASeed = 0;
		for(bool fromBWA : forReturn.chains.first.is_from_BWAseed)
		{
			selectedChains_totalColumns++;
			if(fromBWA)
				selectedChains_totalColumns_fromBWASeed++;
		}
		for(bool fromBWA : forReturn.chains.second.is_from_BWAseed)
		{
			selectedChains_totalColumns++;
			if(fromBWA)
				selectedChains_totalColumns_fromBWASeed++;
		}

		int selectedPair_removed_columns_noGap_restriction =  forReturn.chains.first.removed_columns_noGap_restriction + forReturn.chains.second.removed_columns_noGap_restriction;
		double selectedPair_improvement_through_bt =  forReturn.chains.first.improvement_through_bt + forReturn.chains.second.improvement_through_bt;
		if(!(selectedPair_removed_columns_noGap_restriction >= 0))
		{
			std::cerr << "forReturn.chains.first.removed_columns_noGap_restriction" << ": " << forReturn.chains.first.removed_columns_noGap_restriction  << "\n";
			std::cerr << "forReturn.chains.second.removed_columns_noGap_restriction" << ": " << forReturn.chains.second.removed_columns_noGap_restriction  << "\n";
			std::cerr << std::flush;
		}
		assert(selectedPair_removed_columns_noGap_restriction >= 0);

		statisticsStore->alignOneReadPair_considered_chains += considered_chains;
		statisticsStore->alignOneReadPair_n_calledChainExtension += n_calledChainExtension;
		statisticsStore->alignOneReadPair_n_calledChainExtension_primary += n_calledChainExtension_primary;

		statisticsStore->alignOneReadPair_return_selectedChains_totalColumns += selectedChains_totalColumns;
		statisticsStore->alignOneReadPair_return_selectedChains_totalColumns_fromBWASeed += selectedChains_totalColumns_fromBWASeed;

		statisticsStore->alignOneReadPair_selectedPair_removed_columns_noGap_restriction += selectedPair_removed_columns_noGap_restriction;

		statisticsStore->alignOneReadPair_selectedPair_improvement_through_bt += selectedPair_improvement_through_bt;


	}


	return forReturn;
}

std::pair<int, int> processBAM::alignment_get_startstop_PRGcoordinates(const std::tuple<std::string, int,  BamTools::BamAlignment, int>& al) const
{
	//const std::string& BAMid = (_currentBAM != "") ? R.GetReferenceData().at(al.second.RefID).RefName : R1.GetReferenceData().at(al.second.RefID).RefName;

	const std::string& BAMid = std::get<0>(al);

	if(_currentBAM == "")
	{
		assert(_currentBAM1 != "");
		assert(R2.GetReferenceData().at(std::get<2>(al).RefID).RefName == BAMid);
	}
	assert(BAMid_2_PRGid.count(BAMid));
	int PRGid = BAMid_2_PRGid.at(BAMid);
	assert(extendedReferenceGenome_levelTranslation.count(PRGid));

	int readStart = std::get<2>(al).Position;
	int readStop = std::get<2>(al).GetEndPosition(false, true);

	auto transform_sequenceCoordinate_2_PRGcoordinate = [&](int levelAlongReferenceSequence) -> int {
		int reference2level_offset = std::get<1>(al);
		levelAlongReferenceSequence -= reference2level_offset;
		if(!(levelAlongReferenceSequence < (int)extendedReferenceGenome_levelTranslation.at(PRGid).size()))
		{
			std::cerr << "Error!" << "\n";
			std::cerr << "\tal.RefID: " << std::get<2>(al).RefID << "\n";
			std::cerr << "\t" << "reference2level_offset: " <<  reference2level_offset << "\n";
			std::cerr << "\t" << "levelAlongReferenceSequence" << ": " << levelAlongReferenceSequence << "\n";
			std::cerr << "\t" << "extendedReferenceGenome_levelTranslation.at(PRGid).size()" << ": " << extendedReferenceGenome_levelTranslation.at(PRGid).size() << "\n" << std::flush;

			reads::protoSeeds::printAlignmentInfo(1, std::get<2>(al));
		}
		assert(levelAlongReferenceSequence >= 0);
		assert(levelAlongReferenceSequence < (int)extendedReferenceGenome_levelTranslation.at(PRGid).size());
		int levelAlongPRG = extendedReferenceGenome_levelTranslation.at(PRGid).at(levelAlongReferenceSequence);
		return levelAlongPRG;
	};

	int readStart_PRG = transform_sequenceCoordinate_2_PRGcoordinate(readStart);
	int readStop_PRG = transform_sequenceCoordinate_2_PRGcoordinate(readStop);

	return make_pair(readStart_PRG, readStop_PRG);

}

void processBAM::assignMappingQualities(reads::verboseSeedChainPair& forReturn, const std::vector<std::pair<unsigned int, unsigned int>> combinations_indices, const std::vector<double>& combinations_LL, const std::pair<double, unsigned int>& combinations_max, const std::vector<reads::verboseSeedChain>& read1_extendedChains, const std::vector<reads::verboseSeedChain>& read2_extendedChains)
{
	if(combinations_indices.size() > 1)
	{
		double max_LL = combinations_max.first;
		unsigned int maxCombination_i1 = combinations_indices.at(combinations_max.second).first;
		unsigned int maxCombination_i2 = combinations_indices.at(combinations_max.second).second;

		std::vector<double> combinations_PP = combinations_LL;
		for(unsigned int i = 0; i < combinations_PP.size(); i++)
		{
			double combiantion_ll = combinations_PP.at(i);
			combinations_PP.at(i) = exp(combiantion_ll - max_LL);
			if(!(combinations_PP.at(i) >= 0))
			{
				std::cerr << "! assert(combinations_PP.at(i) >= 0);";
				std::cerr << "\t" << "max_LL" << ":" << max_LL<< "\n";
				std::cerr << "\t" << "combiantion_ll" << ":" << combiantion_ll<< "\n";
				std::cerr << "\n" << std::flush;
			}
			assert(combinations_PP.at(i) >= 0);
			assert(combinations_PP.at(i) <= 1);
		}
		Utilities::normalize_vector(combinations_PP);

		forReturn.chains.first = read1_extendedChains.at(maxCombination_i1);
		forReturn.chains.second = read2_extendedChains.at(maxCombination_i2);
		double mapQ = combinations_PP.at(combinations_max.second);
		forReturn.mapQ = mapQ;

		double mapQ_chain1 = 0;
		double mapQ_chain2 = 0;
		for(unsigned int i = 0; i < combinations_PP.size(); i++)
		{
			double p = combinations_PP.at(i);
			if(combinations_indices.at(i).first == maxCombination_i1)
			{
				mapQ_chain1 += p;
			}
			if(combinations_indices.at(i).second == maxCombination_i2)
			{
				mapQ_chain2 += p;
			}
		}
		assert(mapQ_chain1 >= mapQ);
		assert(mapQ_chain2 >= mapQ);

		// std::cerr << mapQ << "\t" << mapQ_chain1 << "\t" << mapQ_chain2 << "\n" << std::flush;

		assert(mapQ_chain1 >= 0);
		assert(mapQ_chain2 >= 0);
		assert(mapQ_chain1 <= (1+1e-3));
		assert(mapQ_chain2 <= (1+1e-3));

		if(mapQ_chain1 > 1)
		{
			mapQ_chain1 = 1;
		}
		if(mapQ_chain2 > 1)
		{
			mapQ_chain2 = 1;
		}

		forReturn.chains.first.mapQ = mapQ_chain1;
		forReturn.chains.second.mapQ = mapQ_chain2;

		std::map<std::string, double> alignmentPositionConfidences;

		auto alignedSequence2SequenceIndex = [](std::string sequence_aligned, bool reverse) -> std::vector<int> {
			std::vector<int> forReturn;
			forReturn.reserve(sequence_aligned.length());
			int i_noGap = -1;
			int sequenceLength_noGap = 0;
			for(unsigned int i = 0; i < sequence_aligned.length(); i++)
			{
				char sequenceAlignmentChar = sequence_aligned.at(i);
				if(sequenceAlignmentChar != '_')
				{
					sequenceLength_noGap++;
				}
			}

			for(unsigned int i = 0; i < sequence_aligned.length(); i++)
			{
				char sequenceAlignmentChar = sequence_aligned.at(i);
				int sequenceIndex;
				if(sequenceAlignmentChar == '_')
				{
					sequenceIndex = -1;
				}
				else
				{
					i_noGap++;
					if(reverse)
					{
						sequenceIndex = sequenceLength_noGap - i_noGap - 1;
					}
					else
					{
						sequenceIndex = i_noGap;
					}
					assert(sequenceIndex >= 0);
					assert(sequenceIndex < sequenceLength_noGap);
				}

				forReturn.push_back(sequenceIndex);
			}

			return forReturn;
		};

		for(unsigned int i = 0; i < combinations_indices.size(); i++)
		{
			unsigned int i1 = combinations_indices.at(i).first;
			unsigned int i2 = combinations_indices.at(i).second;

			std::vector<int> first_sequence_index = alignedSequence2SequenceIndex(read1_extendedChains.at(i1).sequence_aligned, read1_extendedChains.at(i1).reverse);
			std::vector<int> second_sequence_index = alignedSequence2SequenceIndex(read2_extendedChains.at(i2).sequence_aligned, read2_extendedChains.at(i2).reverse);

			std::string r1_reverse = (read1_extendedChains.at(i1).reverse) ? "minus" : "plus";
			std::string r2_reverse = (read2_extendedChains.at(i2).reverse) ? "minus" : "plus";

			for(unsigned int j = 0; j < read1_extendedChains.at(i1).graph_aligned.length(); j++)
			{
				std::string c_graph = read1_extendedChains.at(i1).graph_aligned.substr(j, 1);
				int l_graph = read1_extendedChains.at(i1).graph_aligned_levels.at(j);
				int idx_sequence = first_sequence_index.at(j);
				std::string positionID = c_graph + ":" + Utilities::ItoStr(l_graph) + ":" +"r1" + ":" + r1_reverse +":" + Utilities::ItoStr(idx_sequence);
				if(alignmentPositionConfidences.count(positionID) == 0)
				{
					alignmentPositionConfidences[positionID] = 0;
				}
				alignmentPositionConfidences.at(positionID) += combinations_PP.at(i);
			}

			for(unsigned int j = 0; j < read2_extendedChains.at(i2).graph_aligned.length(); j++)
			{
				std::string c_graph = read2_extendedChains.at(i2).graph_aligned.substr(j, 1);
				int l_graph = read2_extendedChains.at(i2).graph_aligned_levels.at(j);
				int idx_sequence = second_sequence_index.at(j);
				std::string positionID = c_graph + ":" + Utilities::ItoStr(l_graph) + ":" +"r2" + ":" + r2_reverse + ":" + Utilities::ItoStr(idx_sequence);
				if(alignmentPositionConfidences.count(positionID) == 0)
				{
					alignmentPositionConfidences[positionID] = 0;
				}
				alignmentPositionConfidences.at(positionID) += combinations_PP.at(i);
			}
		}

		auto test_phredConversion = [](double p, double tolerance)
		{
			assert(p >= 0);
			assert(p <= 1);
			char phred = Utilities::PCorrectToPhred(p);
			double pBack = Utilities::PhredToPCorrect(phred);
			if(!(abs(p - pBack) <= pBack))
			{
				std::cerr << "p = " << p << "\n";
				std::cerr << "phred = " << (int)phred << "\n";
				std::cerr << "pBack = " << pBack << "\n";
				std::cerr << "abs(p - pBack) = " << abs(p - pBack) << "\n";
				std::cerr << "tolerance: " << tolerance << "\n";
				std::cerr << std::flush;
			}
			assert(abs(p - pBack) <= pBack);
		};

		test_phredConversion(0.0, 1e-5);
		test_phredConversion(0.5, 1e-1);
		test_phredConversion(0.6, 1e-1);
		test_phredConversion(0.7, 1e-1);
		test_phredConversion(0.8, 1e-2);
		test_phredConversion(0.9, 1e-2);
		test_phredConversion(0.9, 1e-2);
		test_phredConversion(0.99, 1e-3);
		test_phredConversion(0.999, 1e-3);
		test_phredConversion(0.9999, 1e-4);
		test_phredConversion(1, 1e-5);

		// actual alignment
		{
			std::vector<int> first_sequence_index = alignedSequence2SequenceIndex(read1_extendedChains.at(maxCombination_i1).sequence_aligned, read1_extendedChains.at(maxCombination_i1).reverse);
			std::vector<int> second_sequence_index = alignedSequence2SequenceIndex(read2_extendedChains.at(maxCombination_i2).sequence_aligned, read2_extendedChains.at(maxCombination_i2).reverse);

			std::string r1_reverse = (read1_extendedChains.at(maxCombination_i1).reverse) ? "minus" : "plus";
			std::string r2_reverse = (read2_extendedChains.at(maxCombination_i2).reverse) ? "minus" : "plus";

			forReturn.chains.first.mapQ_perPosition.reserve(read1_extendedChains.at(maxCombination_i1).graph_aligned.length());
			for(unsigned int j = 0; j < read1_extendedChains.at(maxCombination_i1).graph_aligned.length(); j++)
			{
				std::string c_graph = read1_extendedChains.at(maxCombination_i1).graph_aligned.substr(j, 1);
				int l_graph = read1_extendedChains.at(maxCombination_i1).graph_aligned_levels.at(j);
				int idx_sequence = first_sequence_index.at(j);
				std::string positionID = c_graph + ":" + Utilities::ItoStr(l_graph) + ":" +"r1" + ":" + r1_reverse +":" + Utilities::ItoStr(idx_sequence);

				double Q = alignmentPositionConfidences.at(positionID);
				if(!((Q - 1) <= 1e-5))
				{
					std::cerr << "Q: " << Q << "\n";
					std::cerr << "positionID: " << positionID << "\n";
					std::cerr << "j: " << j << "\n";
					std::cerr << "forReturn.at(i).first.graph_aligned.length(): " << read1_extendedChains.at(maxCombination_i1).graph_aligned.length() << "\n";
					std::cerr << std::flush;
				}
				assert((Q - 1) <= 1e-5);
				if(Q > 1)
				{
					Q = 1;
				}
				assert(Q >= mapQ);
				forReturn.chains.first.mapQ_perPosition.push_back(Utilities::PCorrectToPhred(Q));
			}


			forReturn.chains.second.mapQ_perPosition.reserve(read2_extendedChains.at(maxCombination_i2).graph_aligned.length());
			for(unsigned int j = 0; j < read2_extendedChains.at(maxCombination_i2).graph_aligned.length(); j++)
			{
				std::string c_graph = read2_extendedChains.at(maxCombination_i2).graph_aligned.substr(j, 1);
				int l_graph = read2_extendedChains.at(maxCombination_i2).graph_aligned_levels.at(j);
				int idx_sequence = second_sequence_index.at(j);
				std::string positionID = c_graph + ":" + Utilities::ItoStr(l_graph) + ":" +"r2" + ":" + r2_reverse + ":" + Utilities::ItoStr(idx_sequence);
				double Q = alignmentPositionConfidences.at(positionID);
				if(!((Q - 1) <= 1e-5))
				{
					std::cerr << "Q: " << Q << "\n";
					std::cerr << "positionID: " << positionID << "\n";
					std::cerr << "j: " << j << "\n";
					std::cerr << "read2_extendedChains.at(maxCombination_i2).graph_aligned.length(): " << read2_extendedChains.at(maxCombination_i2).graph_aligned.length() << "\n";
					std::cerr << std::flush;
				}
				assert((Q - 1) <= 1e-5);
				if(Q > 1)
				{
					Q = 1;
				}
				assert(Q >= mapQ);
				forReturn.chains.second.mapQ_perPosition.push_back(Utilities::PCorrectToPhred(Q));
			}
		}
	}
	else
	{
		forReturn.mapQ = 1;
		forReturn.chains.first.mapQ = 1;
		forReturn.chains.second.mapQ = 1;

		char Phred1 = Utilities::PCorrectToPhred(1);
		forReturn.chains.first.mapQ_perPosition.resize(forReturn.chains.first.graph_aligned.size(), Phred1);
		forReturn.chains.second.mapQ_perPosition.resize(forReturn.chains.second.graph_aligned.size(), Phred1);
	}
}

int processBAM::getAlignmentScore(const BamTools::BamAlignment& al)
{
	
	std::string tagName = "AS";
	u_int32_t v_u;
	int32_t v;
	if(al.GetTag(tagName, v))
	{
		return v;
	}
	else if(al.GetTag(tagName, v_u))
	{
		return v_u;
	}
	else
	{
		std::cerr << "Can't get AS tag!\n" << std::flush;
		reads::protoSeeds::printAlignmentInfo(1, al);
		assert(1 == 0);
	}

	return v;
}

void processBAM::updateBAMSelector()
{
	std::string regionID =(*interestingIntervals_iterator).first;
	int refIDidx = R.GetReferenceID(regionID);
	assert(refIDidx != -1);

	int left = (*interestingIntervals_iterator).second.at(interestingIntervals_i).start_1based - 1;
	int right = (*interestingIntervals_iterator).second.at(interestingIntervals_i).stop_1based - 1;
	assert(left >= 0);
	assert(left < right);

	BamTools::BamRegion stretch;
	stretch.LeftRefID = refIDidx;
	stretch.LeftPosition = left;
	stretch.RightRefID = refIDidx;
	stretch.RightPosition =  right + 1;

	bool success = R.SetRegion(stretch);
	assert(success);

	// std::cout << "\tprocessBAM::updateBAMSelector(): " << regionID << ":" << left << "-" << right << "\n" << std::flush;
}

void processBAM::updateBAMSelectors()
{
	std::string regionID =(*interestingIntervals_iterator).first;
	int refIDidx = R1.GetReferenceID(regionID);
	assert(refIDidx == R2.GetReferenceID(regionID));
	
	assert(refIDidx != -1);

	int left = (*interestingIntervals_iterator).second.at(interestingIntervals_i).start_1based - 1;
	int right = (*interestingIntervals_iterator).second.at(interestingIntervals_i).stop_1based - 1;
	assert(left >= 0);
	assert(left < right);

	BamTools::BamRegion stretch;
	stretch.LeftRefID = refIDidx;
	stretch.LeftPosition = left;
	stretch.RightRefID = refIDidx;
	stretch.RightPosition =  right + 1;

	bool success1 = R1.SetRegion(stretch);
	assert(success1);

	bool success2 = R2.SetRegion(stretch);
	assert(success2);

	// std::cout << "\tprocessBAM::updateBAMSelector(): " << regionID << ":" << left << "-" << right << "\n" << std::flush;
}

void processBAM::_loadMapping(int ID)
{
	if(extendedReferenceGenome_levelTranslation.count(ID) == 0)
	{
		std::string fileName = graphDir + "/translation/" + Utilities::ItoStr(ID) + ".txt";
		if(! Utilities::fileExists(fileName))
		{
			throw std::runtime_error("Expected coordinate translation file not found: "+fileName);
		}

		std::vector<int> translated_levels;

		std::ifstream translationStream;
		translationStream.open(fileName.c_str());
		std::string line;
		assert(translationStream.is_open());

		while(translationStream.good())
		{
			std::getline(translationStream, line);
			Utilities::eraseNL(line);
			translated_levels.push_back(Utilities::StrtoI(line));
		}

		extendedReferenceGenome_levelTranslation[ID] = translated_levels;

//		int firstLevel = translated_levels.front();
//		int lastLevel = translated_levels.back();
//		assert(lastLevel > firstLevel);
//		unsigned int parallelToGraph_length = lastLevel - firstLevel + 1;
//		std::vector<int> translated_levels_parallelToGraph;
//		translated_levels_parallelToGraph.reserve(parallelToGraph_length);
//		int lastInsertedLevel = -1;
//		for(unsigned int levelII = 0; levelII < translated_levels.size(); levelII++)
//		{
//			int thisLevel = translated_levels.at(levelII);
//			if(levelII > 0)
//			{
//				while((lastInsertedLevel + 1) != thisLevel)
//				{
//					lastInsertedLevel++;
//					translated_levels_parallelToGraph.push_back(lastInsertedLevel);
//				}
//			}
//			translated_levels_parallelToGraph.push_back(thisLevel);
//			lastInsertedLevel = thisLevel;
//		}
//		assert(translated_levels_parallelToGraph.size() == parallelToGraph_length);

		for(unsigned int positionAlongRawSequence = 0; positionAlongRawSequence < translated_levels.size(); positionAlongRawSequence++)
		{
			int thisLevel = translated_levels.at(positionAlongRawSequence);
			assert(thisLevel >= 0);
			if(!(thisLevel < graphLevel_2_underlyingSequencePositions.size()))
			{ 
				std::cerr << "ID: " << ID << "\n";
				std::cerr << "thisLevel: " << thisLevel << "\n";
				std::cerr << "graphLevel_2_underlyingSequencePositions.size(: " << graphLevel_2_underlyingSequencePositions.size() << "\n";
				std::cerr << std::flush;
				 
			}	
			assert(thisLevel < graphLevel_2_underlyingSequencePositions.size());
			
			graphLevel_2_underlyingSequencePositions.at((unsigned int)thisLevel)[ID] = positionAlongRawSequence;
		}
		
	}
}

void processBAM::restrictInitialAlignmentToNoGapAreas(std::vector<int>& graph_aligned_levels, std::string& graph_aligned, std::string& sequence_aligned, int& graphSpace_sequence_aligned_startInRaw, int& graphSpace_sequence_aligned_stopInRaw, const std::vector<bool>& inGraphGapStretch)
{
	int runningStretch_nonGap_begin = -1;
	int sequenceCharacters = 0;
	std::vector<std::pair<int, int>> possibleStretches_preFilter;
	for(unsigned int lI = 0; lI < graph_aligned_levels.size(); lI++)
	{
		if(sequence_aligned.at(lI) != '_')
		{
			sequenceCharacters++;
		}

		int graph_level = graph_aligned_levels.at(lI);
		if((graph_level != -1) && (inGraphGapStretch.at(graph_level)))
		{
			if(runningStretch_nonGap_begin != -1)
			{
				int runningStretch_nonGap_end = lI - 1;
				assert(runningStretch_nonGap_end >= runningStretch_nonGap_begin);
				possibleStretches_preFilter.push_back(make_pair(runningStretch_nonGap_begin, runningStretch_nonGap_end));
				runningStretch_nonGap_begin = -1;
			}
		}
		else
		{
			if(runningStretch_nonGap_begin == -1)
			{
				runningStretch_nonGap_begin = lI;
			}
		}
	}
	if(runningStretch_nonGap_begin != -1)
	{
		if(runningStretch_nonGap_begin != 0)
		{
			int runningStretch_nonGap_end = (graph_aligned_levels.size() - 1);
			if(runningStretch_nonGap_end >= runningStretch_nonGap_begin)
			{
				possibleStretches_preFilter.push_back(make_pair(runningStretch_nonGap_begin, runningStretch_nonGap_end));
			}
		}
	}

	std::vector<std::pair<int, int>> possibleStretches;
	possibleStretches.reserve(possibleStretches_preFilter.size());

	for(unsigned int sI = 0; sI < possibleStretches_preFilter.size(); sI++)
	{
		std::pair<int, int> thisStretch = possibleStretches_preFilter.at(sI);

		while(graph_aligned_levels.at(thisStretch.first) == -1)
		{
			thisStretch.first++;
			if((thisStretch.first > thisStretch.second) || (thisStretch.first > ((int)graph_aligned_levels.size() - 1)))
			{
				break;
			}
		}

		while(graph_aligned_levels.at(thisStretch.second) == -1)
		{
			thisStretch.second--;
			if((thisStretch.second < thisStretch.first) || (thisStretch.second < 0))
			{
				break;
			}
		}

		if(thisStretch.second >= thisStretch.first)
		{
			possibleStretches.push_back(thisStretch);
		}
	}


//	std::cerr << "Initial alignment of " << sequenceCharacters << " characters, found " << possibleStretches.size() << " possible stretches.\n" << std::flush;

	std::sort(possibleStretches.begin(), possibleStretches.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b){
		int aL = a.second - a.first + 1;
		int bL = b.second - b.first + 1;
		assert(aL >= 1);
		assert(bL >= 1);
		return (aL < bL);
	});

	if(possibleStretches.size() > 1)
	{
		int lastL = possibleStretches.back().second - possibleStretches.back().first + 1;
		assert(lastL >= 1);
		int secondLastL = possibleStretches.at(possibleStretches.size() - 1).second - possibleStretches.at(possibleStretches.size() - 1).first + 1;
		assert(secondLastL >= 1);
		assert(lastL >= secondLastL);
	}

	if(possibleStretches.size() > 0)
	{
		std::pair<int, int> selectedStretch = possibleStretches.back();
		std::vector<int> new_graph_aligned_levels;
		std::string  new_graph_aligned;
		std::string new_sequence_aligned;

		int stretchL = selectedStretch.second - selectedStretch.first + 1;
		new_graph_aligned_levels.reserve(stretchL);
		new_graph_aligned.reserve(stretchL);
		new_graph_aligned.reserve(stretchL);

		int new_graphSpace_sequence_aligned_startInRaw = graphSpace_sequence_aligned_startInRaw;
		int new_graphSpace_sequence_aligned_stopInRaw = graphSpace_sequence_aligned_stopInRaw;

		if(selectedStretch.first != 0)
		{
			for(int lI = 0; lI < selectedStretch.first; lI++)
			{
				if(sequence_aligned.at(lI) !='_')
				{
					new_graphSpace_sequence_aligned_startInRaw++;
				}
			}
		}
		if(selectedStretch.second != (int)(graph_aligned_levels.size() - 1))
		{
			for(int lI = (selectedStretch.second + 1); lI < (int)graph_aligned_levels.size(); lI++)
			{
				if(sequence_aligned.at(lI) !='_')
				{
					new_graphSpace_sequence_aligned_stopInRaw--;
				}
			}
		}

		int stretch_sequenceCharacters = 0;
		for(int lI = selectedStretch.first; lI <= selectedStretch.second; lI++)
		{

			new_graph_aligned_levels.push_back(graph_aligned_levels.at(lI));
			new_graph_aligned.push_back(graph_aligned.at(lI));
			new_sequence_aligned.push_back(sequence_aligned.at(lI));
			if(sequence_aligned.at(lI) != '_')
			{
				stretch_sequenceCharacters++;
			}
		}

//		std::cerr << "\tSelected stretch of length " << stretchL << " with " << stretch_sequenceCharacters << " characters.\n" << std::flush;

		if(((double)stretch_sequenceCharacters/(double)sequenceCharacters) > 0.3)
		{
			graph_aligned_levels = new_graph_aligned_levels;
			graph_aligned = new_graph_aligned;
			sequence_aligned = new_sequence_aligned;
			graphSpace_sequence_aligned_startInRaw = new_graphSpace_sequence_aligned_startInRaw;
			graphSpace_sequence_aligned_stopInRaw = new_graphSpace_sequence_aligned_stopInRaw;
		}
		else
		{
//			std::cerr << "\tDISCARD, not enough characters!" << std::flush;
		}
	}
}

void processBAM::cleanInitialAlignment(std::vector<int>& graph_aligned_levels, std::string& graph_aligned, std::string& sequence_aligned)
{
	bool paranoid = true;

	// we want to remove sequences of sequence insertions, graph gaps
	std::vector<int> potentialBreakPoints; // if value x is in here, the breakpoint is between string coordinates x - 1 and x
	bool inInterestingStretch = false;
	int interestingStretchStart = -1;
	bool sawSequenceInsertion = false;
	bool sawGraphInternalGap = false;
	int balance = 0;
	bool cleaningNecessary = false;

	std::vector<int> levelsBeforeCleaning;
	std::string graphBeforeCleaning;
	std::string sequenceBeforeCleaning;

	if(paranoid)
	{
		levelsBeforeCleaning = graph_aligned_levels;
		graphBeforeCleaning = graph_aligned;
		sequenceBeforeCleaning = sequence_aligned;
	}

	for(unsigned int pI = 0; pI < graph_aligned_levels.size(); pI++)
	{
		if((graph_aligned_levels.at(pI) == -1) || ((graph_aligned.at(pI) == '_') && (sequence_aligned.at(pI) == '_')))
		{
			if(inInterestingStretch == false)
			{
				interestingStretchStart = pI;
				inInterestingStretch = true;
			}

			if((graph_aligned_levels.at(pI) == -1))
			{
				sawSequenceInsertion = true;
				balance++;
			}

			if(((graph_aligned.at(pI) == '_')) && (sequence_aligned.at(pI) == '_'))
			{
				assert(graph_aligned_levels.at(pI) != -1);
				sawGraphInternalGap = true;
				balance--;
			}

		}
		else
		{
			if(inInterestingStretch)
			{
				int interestingStretchStop = pI - 1;

				// std::cout << "Stretch from " << interestingStretchStart << " to " << interestingStretchStop << ", balance " << balance << "\n" << std::flush;

				if(balance == 0)
				{
					int L = interestingStretchStop - interestingStretchStart + 1;
					assert(L >= 2);
					assert((L % 2) == 0);
					std::string insertedSequenceCharacters;
					std::vector<int> graphInternalGapLevels;
					for(int pII = interestingStretchStart; pII <= interestingStretchStop; pII++)
					{
						if((graph_aligned_levels.at(pII) == -1))
						{
							assert(sequence_aligned.at(pII) != '_');
							insertedSequenceCharacters.push_back(sequence_aligned.at(pII));
						}
						else
						{
							assert(((graph_aligned.at(pII) == '_')) && (sequence_aligned.at(pII) == '_'));
							assert(graph_aligned_levels.at(pII) != -1);
							graphInternalGapLevels.push_back(graph_aligned_levels.at(pII));
						}
					}

					assert(insertedSequenceCharacters.size() == graphInternalGapLevels.size());
					assert((int)insertedSequenceCharacters.size() == (L / 2));

					cleaningNecessary = true;

					for(int pII = interestingStretchStart; pII <= interestingStretchStop; pII++)
					{
						int iCharacter = pII - interestingStretchStart;
						if(iCharacter < (L/2))
						{
							graph_aligned_levels.at(pII) = graphInternalGapLevels.at(iCharacter);
							graph_aligned.at(pII) = '_';
							sequence_aligned.at(pII) = insertedSequenceCharacters.at(iCharacter);
						}
						else
						{
							graph_aligned_levels.at(pII) = -1;
							graph_aligned.at(pII) = '_';
							sequence_aligned.at(pII) = '_';
						}
					}
				}
				else
				{
					potentialBreakPoints.push_back(interestingStretchStart);
					potentialBreakPoints.push_back(interestingStretchStop);
				}

				inInterestingStretch = false;
				sawSequenceInsertion = false;
				sawGraphInternalGap = false;
				interestingStretchStart = -1;
				balance = 0;
			}
		}
	}

	if(cleaningNecessary)
	{
		std::vector<int> new_graph_aligned_levels;
		std::string new_graph_aligned;
		std::string new_sequence_aligned;

		new_graph_aligned_levels.reserve(graph_aligned_levels.size());
		new_graph_aligned.reserve(graph_aligned_levels.size());
		new_sequence_aligned.reserve(graph_aligned_levels.size());

		int deletedPositions = 0;
		for(unsigned int pI = 0; pI < graph_aligned_levels.size(); pI++)
		{
			if( (graph_aligned_levels.at(pI) == -1) &&
				(graph_aligned.at(pI) == '_') &&
				(sequence_aligned.at(pI) == '_')
			)
			{
				deletedPositions++;
			}
			else
			{
				new_graph_aligned_levels.push_back(graph_aligned_levels.at(pI));
				new_graph_aligned.push_back(graph_aligned.at(pI));
				new_sequence_aligned.push_back(sequence_aligned.at(pI));
			}
		}

		assert(deletedPositions > 0);
		graph_aligned_levels = new_graph_aligned_levels;
		graph_aligned = new_graph_aligned;
		sequence_aligned = new_sequence_aligned;
	}


	if(cleaningNecessary && 0)
	{

		std::cout << "Extracted alignment BEFORE modification:\n";
		std::cout << "\tlevelsBeforeCleaning                : " << Utilities::join(Utilities::ItoStr(levelsBeforeCleaning), ", ") << "\n";
		std::cout << "\tgraphBeforeCleaning                : " << graphBeforeCleaning << "\n";
		std::cout << "\tsequenceBeforeCleaning                : " << sequenceBeforeCleaning << "\n";

		std::cout << "Extracted alignment AFTER modification:\n";
		std::cout << "\tLevels                : " << Utilities::join(Utilities::ItoStr(graph_aligned_levels), ", ") << "\n";
		std::cout << "\tGraph_aligned         : " << graph_aligned << "\n";
		std::cout << "\tSequence_aligned      : " << sequence_aligned << "\n";

		std::cout << std::flush;
	}

	if(paranoid)
	{
		assert(Utilities::removeGaps(graph_aligned) == Utilities::removeGaps(graphBeforeCleaning));
		assert(Utilities::removeGaps(sequence_aligned) == Utilities::removeGaps(sequenceBeforeCleaning));
	}
}

bool processBAM::transformBAMreadToInternalAlignment(const std::string& referenceSequence, const std::vector<int>& reference2level, int reference2level_offset, const BamTools::BamAlignment& al, const std::string& queryBases, const std::string& qualitiesString, reads::PRGContigBAMAlignment& forReturn)
{
	if(! al.IsMapped())
	{
		return false;
	}

	forReturn.graph_aligned_levels.clear();
	forReturn.graph_aligned.clear();
	forReturn.sequence_aligned.clear();

	forReturn.sequence_aligned_startInRaw = -1;
	forReturn.sequence_aligned_stopInRaw = -1;

	int readStart = al.Position;

	std::string CIGARstring;
	std::vector< BamTools::CigarOp > CIGAR_Compressed = al.CigarData;
	std::vector<char> CIGAR;
	CIGAR.reserve(al.AlignedBases.length()*1.1);
	for(unsigned int cigarI = 0; cigarI < CIGAR_Compressed.size(); cigarI++)
	{
		BamTools::CigarOp thisOp = CIGAR_Compressed.at(cigarI);
		if(thisOp.Type != 'P')
		{
			for(unsigned int thisOpChar = 0; thisOpChar < thisOp.Length; thisOpChar++)
			{
				CIGAR.push_back(thisOp.Type);
			}
		}

		std::string thisCigarOP = Utilities::ItoStr(thisOp.Length);
		thisCigarOP.push_back(thisOp.Type);
		CIGARstring.append(thisCigarOP);
	}

	// We modify these indexes as we walk along the read
	int index_along_genome_fromReadStart = 0; // how far have we got along the genome?
	int index_along_read = 0; // how far have we got along the read as reported by BAMTools' AlignedBases?
							  // (which ignores hard clipping and displays deletions as "-")
	int index_along_unaligned_read = 0; // how far have we got along the unaligned (but potentially hard-clipped) read?
	int index_along_unclipped_read = 0; // how far have we got along the unclipped, unaligned (ie raw) read?

	const std::string& alignedBases = al.AlignedBases;
	// std::string alignedBases = "";
	// const std::string& queryBases = al.QueryBases;
	// const std::string& qualitiesString = al.Qualities;

	// populate read_forLL
	/*
	read_forLL.name = al.Name;
	if(softClipped)
	{
		read_forLL.name = "SOFTCLIPPED";
	}
	read_forLL.sequence = queryBases;
	read_forLL.quality = qualitiesString;
	if ( al.IsReverseStrand() ) {
		std::reverse(read_forLL.quality.begin(), read_forLL.quality.end());
		read_forLL.sequence = Utilities::seq_reverse_complement(read_forLL.sequence);
	}
	*/

	// population alignment_forLL
	forReturn.reverse = al.IsReverseStrand();

	forReturn.graph_aligned_levels.reserve(CIGAR.size());
	forReturn.graph_aligned.reserve(CIGAR.size());
	forReturn.sequence_aligned.reserve(CIGAR.size());

	int readLength = queryBases.length(); // should NOT include hard clipped bases!

	// If we begin with hard clipping, we need to increment both readLength and index_along_unclipped_read...
	assert(CIGAR_Compressed.size() >= 1);
	if(CIGAR_Compressed.at(0).Type == 'H')
	{
		assert(CIGAR_Compressed.size() >= 2);
		assert(CIGAR_Compressed.at(1).Type != 'H');
		readLength += CIGAR_Compressed.at(0).Length;
		index_along_unclipped_read += CIGAR_Compressed.at(0).Length;
	}

	// .. and if hard clipping at the end, we need to account for that by increasing read length...
	if(CIGAR_Compressed.at(CIGAR_Compressed.size() - 1).Type == 'H')
	{
		assert(((int)CIGAR_Compressed.size() - 2) >= 0);
		assert(CIGAR_Compressed.at(CIGAR_Compressed.size() - 2).Type != 'H');
		readLength += CIGAR_Compressed.at(CIGAR_Compressed.size() - 1).Length;
	}

	// as many quality characters as bases in the raw output!
	assert(qualitiesString.length() == queryBases.length());

	std::string alleleForPosition_allele;
	std::vector<char> alleleForPosition_qualities;
	std::string alleleForPosition_genome;
	std::vector<int> alleleForPosition_genome_graphLevel;

	int alleleForPosition_genomicPosition;

	std::string alleleForPositionM1_allele;
	std::vector<char> alleleForPositionM1_qualities;
	std::string alleleForPositionM1_genome;
	std::vector<int> alleleForPositionM1_genome_graphLevel;

	int alleleForPositionM1_genomicPosition;

	bool debug = false;
	std::vector<std::string> debug_sequence_aligned;
	if(debug)
		debug_sequence_aligned.resize(CIGAR.size()+1);

	// we move along the CIGAR string and reconstruct our columns...
	for(unsigned int cigarI = 0; cigarI < CIGAR.size(); cigarI++)
	{
		int index_into_genome = readStart + index_along_genome_fromReadStart;
		int index_into_results_vector = index_into_genome - readStart;
		int index_into_results_vectorM1;

		assert(index_into_results_vector >= 0);

		if(debug)
		{
			if(index_into_results_vector > ((int)debug_sequence_aligned.size() - 1))
			{
				debug_sequence_aligned.resize(index_into_results_vector+2);
			}
		}

		alleleForPosition_allele.clear();
		alleleForPositionM1_allele.clear();
		int alleleForPosition_allele_start_index_along_unclipped_read = -1;

		alleleForPosition_qualities.clear();
		alleleForPositionM1_qualities.clear();

		alleleForPosition_genome.clear();
		alleleForPositionM1_genome.clear();

		alleleForPosition_genome_graphLevel.clear();
		alleleForPositionM1_genome_graphLevel.clear();
 
		alleleForPosition_genomicPosition = index_into_genome;
		alleleForPositionM1_genomicPosition = index_into_genome - 1;

		// CIGAR operation
		char CIGARoperation = CIGAR.at(cigarI);

		auto printDebug = [&]() {
			std::cerr << "Debug information:\n";
			//std::cerr << "\tal.BuildCharData(): " << al.BuildCharData() << "\n";
			std::cerr << "\tal.Name: " << al.Name << "\n";
			std::cerr << "\tal.Position: " << al.Position << "\n";
			std::cerr << "\tal.IsMapped: " << al.IsMapped() << "\n";
			std::cerr << "\tal.IsFailedQC: " << al.IsFailedQC() << "\n";
			std::cerr << "\tindex_along_read: " << index_along_read << "\n";
			std::cerr << "\tindex_along_genome_fromReadStart: " << index_along_genome_fromReadStart << "\n";
			std::cerr << "\talignedBases: " << alignedBases << "\n";
			std::cerr << "\tqueryBases: " << queryBases << "\n";
			std::cerr << "\tcigarI: " << cigarI << "\n";

			std::cerr << "\tCIGAR: ";
			for(unsigned int cigarI = 0; cigarI < CIGAR_Compressed.size(); cigarI++)
			{
				std::cerr <<  CIGAR_Compressed.at(cigarI).Type << CIGAR_Compressed.at(cigarI).Length << " ";

			}
			std::cerr  << "\n" << std::flush;
		};

		switch(CIGARoperation)
		{
		// all of these operations take one character out of alignedBases
		case 'M':
		case '=':
		case 'X':
		case 'D':
			if(CIGARoperation == 'D')
			{
				if(alignedBases.length())
				{
					assert(alignedBases.substr(index_along_read, 1) == "-");
				}
				alleleForPosition_allele = "_";
				alleleForPosition_qualities.push_back(0);
				alleleForPosition_genome = referenceSequence.substr(index_into_genome, 1);
				alleleForPosition_genome_graphLevel.push_back(index_into_genome);
			}
			else
			{
				// alleleForPosition_allele = alignedBases.substr(index_along_read, 1);

				if(!(
						((index_along_unaligned_read >= 0) && (index_along_unaligned_read < (int)qualitiesString.size())) &&
						((index_along_unaligned_read >= 0) && (index_along_unaligned_read < (int)queryBases.size()))
				))
				{
					reads::protoSeeds::printAlignmentInfo(1, al);
				}
					assert((index_along_unaligned_read >= 0) && (index_along_unaligned_read < (int)qualitiesString.size()));
					assert((index_along_unaligned_read >= 0) && (index_along_unaligned_read < (int)queryBases.size()));

					if(!(
							((index_along_unclipped_read >= 0) && (index_along_unclipped_read < (int)qualitiesString.size())) &&
							((index_along_unclipped_read >= 0) && (index_along_unclipped_read < (int)queryBases.size()))
					))
					{
						std::cout << "index_along_unclipped_read: " << index_along_unclipped_read << "\n" << std::flush;
						reads::protoSeeds::printAlignmentInfo(1, al);
					}

					assert((index_along_unclipped_read >= 0) && (index_along_unclipped_read < (int)qualitiesString.size()));
					assert((index_along_unclipped_read >= 0) && (index_along_unclipped_read < (int)queryBases.size()));


				alleleForPosition_allele = queryBases.substr(index_along_unclipped_read, 1);

				alleleForPosition_qualities.push_back(qualitiesString.at(index_along_unclipped_read));
				alleleForPosition_genome = referenceSequence.substr(index_into_genome, 1);
				alleleForPosition_genome_graphLevel.push_back(index_into_genome);

				assert(alleleForPosition_allele.size() > 0);
				if(alignedBases.length())
				{
					if(!(alleleForPosition_allele.at(0) == alignedBases.at(index_along_read)))
					{
						std::cout << "alleleForPosition_allele.at(0): " << alleleForPosition_allele.at(0) << "\n";
						std::cout << "alignedBases.at(index_along_read): " << alignedBases.at(index_along_read) << "\n";
						std::cout << "index_along_read" << ": " << index_along_read << "\n";
						std::cout << "index_along_unclipped_read" << ": " << index_along_unclipped_read << "\n";
						std::cout << "queryBases.substr(index_along_unclipped_read, 1)" << ": " << queryBases.substr(index_along_unclipped_read, 1) << "\n";
						std::cout << "alignedBases.at(index_along_read)" << ": " << alignedBases.at(index_along_read) << "\n";
						std::cout << "queryBases" << ": " << queryBases << "\n";
						std::cout << "alignedBases" << ": " << alignedBases << "\n";
						
						reads::protoSeeds::printAlignmentInfo(1, al);

						std::cout << std::flush;
					}

					//std::cout << alleleForPosition.allele.at(0) << " " << queryBases.at(index_along_unaligned_read) << "\n";
					assert(alleleForPosition_allele.at(0) == alignedBases.at(index_along_read));
				}
			}

			alleleForPosition_allele_start_index_along_unclipped_read = index_along_unclipped_read;

			// if the next position begins an insertion, all associated nucleotides end up in this column!
			while(((cigarI+1) < CIGAR.size()) && (CIGAR.at(cigarI+1) == 'I'))
			{
				int nextPositionUnaligned = (CIGARoperation == 'D') ? index_along_unaligned_read : (index_along_unaligned_read + 1);
				int nextPositionUnclipped = (CIGARoperation == 'D') ? index_along_unclipped_read : (index_along_unclipped_read + 1);

				std::string forAppend = queryBases.substr(nextPositionUnclipped, 1);
				if(alignedBases.length())
				{
					assert(forAppend == alignedBases.substr(index_along_read+1, 1));
				}

				alleleForPosition_allele.append(forAppend);
				alleleForPosition_qualities.push_back(qualitiesString.at(nextPositionUnclipped));

				assert(alleleForPosition_allele.at(alleleForPosition_allele.size()-1) == queryBases.at(nextPositionUnclipped));

				alleleForPosition_genome.append("_");
				alleleForPosition_genome_graphLevel.push_back(-1);

				index_along_read++;
				index_along_unaligned_read++;
				index_along_unclipped_read++;
				cigarI++;
			}

			if(CIGARoperation == 'D')
			{
				index_along_genome_fromReadStart++;
				index_along_read++;
			}
			else
			{
				index_along_genome_fromReadStart++;
				index_along_read++;
				index_along_unaligned_read++;
				index_along_unclipped_read++;
			}

			break;

		// if we begin with an insertion, we collect all inserted nucleotides and add them to the previous column
		// note that we do not increment index_along_genome_fromReadStart - the reported position for the read in this
		// case should refer to the first (matched) position after the insertion
		case 'I':
			assert(index_along_read == 0);
			// alleleForPositionM1_allele = alignedBases.substr(index_along_read, 1);
			alleleForPositionM1_allele = queryBases.substr(index_along_unclipped_read, 1);
			alleleForPositionM1_qualities.push_back(qualitiesString.at(index_along_unclipped_read));

			assert(forReturn.sequence_aligned_startInRaw == -1);

			forReturn.sequence_aligned_startInRaw = index_along_unclipped_read;

			if(alignedBases.length())
			{
				assert(alleleForPositionM1_allele == alignedBases.substr(index_along_read, 1));
			}

			alleleForPositionM1_genome.append("_");
			alleleForPositionM1_genome_graphLevel.push_back(-1);

			while(((cigarI+1) < CIGAR.size()) && (CIGAR.at(cigarI+1) == 'I'))
			{
				// alleleForPositionM1_allele.append(alignedBases.substr(index_along_read+1, 1));
				alleleForPositionM1_allele.append(queryBases.substr(index_along_unclipped_read+1, 1));

				alleleForPositionM1_qualities.push_back(qualitiesString.at(index_along_unclipped_read+1));

				if(alignedBases.length() > 0)
				{	
					bool b1 = (alignedBases.length() > (index_along_read+1));
					bool b2 = (queryBases.length() > (index_along_unclipped_read+1));
					
					if(!(b1 && b2))
					{
						printDebug();
					}
					assert(b1);
					assert(b2);
					
					assert(alignedBases.substr(index_along_read+1, 1).at(0) == queryBases.at(index_along_unclipped_read+1));			
					// printDebug();
				}
				
				/*
				if(alignedBases.length())
				{
					assert(alleleForPositionM1_allele == alignedBases.substr(index_along_read+1, 1));
				}
				*/
				
				alleleForPositionM1_genome.append("_");
				alleleForPositionM1_genome_graphLevel.push_back(-1);

				index_along_read++; 
				index_along_unaligned_read++;
				index_along_unclipped_read++;
				cigarI++;
			}

			index_along_read++;
			index_along_unaligned_read++;
			index_along_unclipped_read++;

			index_into_results_vectorM1 = index_into_results_vector - 1;

			// We ignore the case in which a read starts with an insertion

			assert(alleleForPositionM1_allele.length() == alleleForPositionM1_genome.length());
			assert(alleleForPositionM1_allele.length() == alleleForPositionM1_genome_graphLevel.size());

			forReturn.sequence_aligned.append(alleleForPositionM1_allele);
			forReturn.graph_aligned.append(alleleForPositionM1_genome);
			forReturn.graph_aligned_levels.insert(forReturn.graph_aligned_levels.end(), alleleForPositionM1_genome_graphLevel.begin(), alleleForPositionM1_genome_graphLevel.end());

			forReturn.sequence_aligned_stopInRaw = index_along_unclipped_read - 1;

			if(debug)
			{
				assert(index_into_results_vectorM1 >= -1);
				assert((index_into_results_vectorM1 + 1) < (int)debug_sequence_aligned.size());
				debug_sequence_aligned.at(index_into_results_vectorM1+1).append(alleleForPositionM1_allele);
			}

			break;
		case 'N':
			throw std::runtime_error("N character in CIGAR - should only be the case for RNASeq data!\n");
			break;
		// soft-clipped reads are not assumed to be aligned, but they appear in SEQ (and in BAMTools' alignedBases,
		// hopefully. (this latter sentence does not seem to hold!)
		case 'S':
			//index_along_genome_fromReadStart++;
			//index_along_read++;
			index_along_unaligned_read++;
			index_along_unclipped_read++;
			break;
		case 'H':
			break;
		case 'P':
			throw std::runtime_error("P character in CIGAR - should have been removed in advance!\n");
			break;
		default:
			throw std::runtime_error("Unknown element of CIGAR string");
		}

		if(alleleForPosition_allele != "")
		{
			assert(alleleForPosition_allele.length() == alleleForPosition_genome.length());
			assert(alleleForPosition_allele.length() == alleleForPosition_genome_graphLevel.size());


			forReturn.sequence_aligned.append(alleleForPosition_allele);
			forReturn.graph_aligned.append(alleleForPosition_genome);
			forReturn.graph_aligned_levels.insert(forReturn.graph_aligned_levels.end(), alleleForPosition_genome_graphLevel.begin(), alleleForPosition_genome_graphLevel.end());

			if(forReturn.sequence_aligned_startInRaw == -1)
			{
				assert(alleleForPosition_allele_start_index_along_unclipped_read != -1);
				// assert(index_along_unclipped_read > 0);
				forReturn.sequence_aligned_startInRaw = alleleForPosition_allele_start_index_along_unclipped_read;
			}
			forReturn.sequence_aligned_stopInRaw = index_along_unclipped_read - 1;

			if(debug)
			{
				assert(index_into_results_vector >= -1);
				assert((index_into_results_vector + 1) < (int)debug_sequence_aligned.size());
				debug_sequence_aligned.at(index_into_results_vector+1).append(alleleForPosition_allele);
			}
		}
	}

	// If everything went well, index_along_read should be equal to length of the alignedBases string --
	// if not, complain and print debug info!
	if(alignedBases.length() && (!(index_along_read == (int)alignedBases.length())))
	{
		std::cerr << "Debug information:\n";
		std::cerr << "\tindex_along_read: " << index_along_read << "\n";
		std::cerr << "\tindex_along_genome_fromReadStart: " << index_along_genome_fromReadStart << "\n";
		std::cerr << "\talignedBases: " << alignedBases << "\n";
		std::cerr << "\tqueryBases: " << queryBases << "\n";

		std::cerr << "\tCIGAR: ";
		for(unsigned int cigarI = 0; cigarI < CIGAR_Compressed.size(); cigarI++)
		{
			std::cerr <<  CIGAR_Compressed.at(cigarI).Type << CIGAR_Compressed.at(cigarI).Length << " ";

		}
		std::cerr  << "\n" << std::flush;
	}
	if(alignedBases.length())
	{
		assert(index_along_read == (int)alignedBases.length());
	}

	if(forReturn.reverse)
	{
//		forReturn.sequence_aligned_startInRaw = readLength - forReturn.sequence_aligned_startInRaw - 1;
//		forReturn.sequence_aligned_stopInRaw = readLength - forReturn.sequence_aligned_stopInRaw - 1;
//		int third = forReturn.sequence_aligned_startInRaw;
//		forReturn.sequence_aligned_startInRaw = forReturn.sequence_aligned_stopInRaw;
//		forReturn.sequence_aligned_stopInRaw = third;
	}
	assert(forReturn.sequence_aligned_startInRaw < forReturn.sequence_aligned_stopInRaw); // could also be <=

	if(debug)
	{
		std::string debug_sequence_aligned_concatenated = Utilities::join(debug_sequence_aligned, "");
		if(!(forReturn.sequence_aligned == debug_sequence_aligned_concatenated))
		{
			std::cerr << "!(alignment_forLL.sequence_aligned == debug_sequence_aligned_concatenated)" << "\n";
			std::cerr << "alignment_forLL.sequence_aligned" << ": " << forReturn.sequence_aligned << "\n";
			std::cerr << "debug_sequence_aligned_concatenated" << ": " << debug_sequence_aligned_concatenated << "\n";
			std::cerr << std::flush;
		}
		assert(forReturn.sequence_aligned == debug_sequence_aligned_concatenated);
	}


	bool haveNonMinusOne = false;
	for(unsigned int i = 0; i <forReturn. graph_aligned_levels.size(); i++)
	{
		if(forReturn.graph_aligned_levels.at(i) != -1)
		{
			haveNonMinusOne = true;
			break;
		}
	}

	if(! haveNonMinusOne)
	{
		std::cerr << "Unexpected problem: apparently all insertions!\n";
		std::cerr << "\tal.RefID: " << al.RefID << "\n";
		std::cerr << "\tPosition: " << al.Position << "\n";
		std::cerr << "\talignedBases: " << alignedBases << "\n";
		std::cerr << "\tqueryBases: " << queryBases << "\n";
		std::cerr << "\tCIGAR: ";
		for(unsigned int cigarI = 0; cigarI < CIGAR_Compressed.size(); cigarI++)
		{
			std::cerr <<  CIGAR_Compressed.at(cigarI).Type << CIGAR_Compressed.at(cigarI).Length << " ";
		}
		std::cerr  << "\n" << std::flush;
		std::cerr << "\n" << std::flush;

		return false;
	}

	size_t positions_nonGap = 0;
	for(size_t columnI = 0; columnI < forReturn.graph_aligned_levels.size(); columnI++)
	{
		int levelAlongReferenceSequence = forReturn.graph_aligned_levels.at(columnI);
		if(levelAlongReferenceSequence != -1)
		{
			levelAlongReferenceSequence -= reference2level_offset;
			if(!(levelAlongReferenceSequence < (int)reference2level.size()))
			{
				std::cerr << "Error!" << "\n";
				std::cerr << "\tforReturn.graph_aligned_levels.at(columnI): " << forReturn.graph_aligned_levels.at(columnI) << "\n";
				std::cerr << "\tal.RefID: " << al.RefID << "\n";
				std::cerr << "\treference2level_offset: " << reference2level_offset << "\n";
				std::cerr << "\t" << "levelAlongReferenceSequence" << ": " << levelAlongReferenceSequence << "\n";
				std::cerr << "\t" << "reference2level.size()" << ": " << reference2level.size() << "\n" << std::flush;
				
				reads::protoSeeds::printAlignmentInfo(1, al);					
			}
			assert(levelAlongReferenceSequence >= 0);
			assert(levelAlongReferenceSequence < (int)reference2level.size());
			int levelAlongPRG = reference2level.at(levelAlongReferenceSequence);
			forReturn.graph_aligned_levels.at(columnI) = levelAlongPRG;
		}

		if(forReturn.sequence_aligned.at(columnI) != '_')
		{
			positions_nonGap++;
		}
	}

	if(! ((long long)positions_nonGap == (forReturn.sequence_aligned_stopInRaw - forReturn.sequence_aligned_startInRaw + 1)))
	{
		std::cerr << "Lenght mismatch!\n";
		std::cerr << "\tpositions_nonGap: " << positions_nonGap << "\n";
		std::cerr << "\tCIGAR: " << CIGARstring << "\n";
		std::cerr << "\tsequence_aligned_stopInRaw: " << forReturn.sequence_aligned_stopInRaw << "\n";
		std::cerr << "\tsequence_aligned_startInRaw: " << forReturn.sequence_aligned_startInRaw << "\n";
		std::cerr << "\t(sequence_aligned_stopInRaw - sequence_aligned_startInRaw + 1): " << (forReturn.sequence_aligned_stopInRaw - forReturn.sequence_aligned_startInRaw + 1) << "\n";
		std::cerr << "\tforReturn.reverse: " << forReturn.reverse << "\n\n";
		std::cerr << Utilities::join(Utilities::ItoStr(forReturn.graph_aligned_levels), ", ") << "\n";
		std::cerr << forReturn.graph_aligned << "\n";
		std::cerr << forReturn.sequence_aligned << "\n";

		std::cerr << std::flush;
	}
	assert((long long)positions_nonGap == (forReturn.sequence_aligned_stopInRaw - forReturn.sequence_aligned_startInRaw + 1));

	return true;
}

Graph* processBAM::getGraph()
{
	return g;
}



} /* namespace mapper */
