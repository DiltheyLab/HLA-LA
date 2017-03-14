/*
< * simpleGraphSimulator.cpp
 *
 *  Created on: 22.09.2015
 *      Author: AlexanderDilthey
 */

#include <fstream>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <iostream>

#include "../../mapper/bwa/BWAmapper.h"
#include "../graphSimulator/simpleGraphSimulator.h"

#include "../../simulator/readSimulator.h"
#include "../../mapper/reads/oneReadPair.h"

simpleGraphSimulator::simpleGraphSimulator(bool noGaps_switchContigs_, const pathFinder& pF) : pF(pF) {
	noGaps_switchContigs = noGaps_switchContigs_;
	graphLength = 25000;
	forPRG_mutations = 2;
	forPRG_largeGaps = 1;
	additionalForReads_mutations = 2;
	mutationDensity = 0.02;
	gapStartFrequency = (noGaps_switchContigs) ? 0 : 0.01;
	gapExpectedLength = 10;
	g = 0;

	init();
}

Graph* simpleGraphSimulator::getGraph()
{
	return g;
}



void simpleGraphSimulator::storeLikeRealPRG(std::string directory)
{
	Utilities::make_or_clearDirectory(directory + "/PRG");
	g->writeToFile(directory + "/PRG/graph.txt");

	Utilities::make_or_clearDirectory(directory + "/mapping");
	Utilities::make_or_clearDirectory(directory + "/mapping_PRGonly");

	Utilities::make_or_clearDirectory(directory + "/translation");
	Utilities::make_or_clearDirectory(directory + "/referenceGenome");

	std::string sequences_output_fn = directory + "/sequences.txt";
	std::ofstream sequencesStream;
	sequencesStream.open(sequences_output_fn.c_str());
	assert(sequencesStream.is_open());
	sequencesStream << Utilities::join({"SequenceID", "Name", "FASTAID", "Chr", "Start_1based", "Stop_1based"}, "\t") << "\n";

	std::map<std::string, std::string> PRG_underlyingSequences = hP.getAllHaplotypes();
	std::map<std::string, std::string> PRG_underlyingSequences_noGaps;
	std::vector<std::string> ordered_sequenceIDs;
	for(std::map<std::string, std::string>::iterator seqIt = PRG_underlyingSequences.begin(); seqIt != PRG_underlyingSequences.end(); seqIt++)
	{
		std::string originalS = seqIt->second;
		std::string S = Utilities::removeGaps(originalS);
		PRG_underlyingSequences_noGaps[seqIt->first] = S;

		std::map<std::string, std::string> oneSeq;
		oneSeq[seqIt->first] = S;

		ordered_sequenceIDs.push_back(seqIt->first);

		Utilities::writeFASTA(directory + "/mapping/" + Utilities::ItoStr(ordered_sequenceIDs.size()) + ".fa", oneSeq);

		std::vector<int> translatedPositions;
		for(size_t posI = 0; posI < originalS.length(); posI++)
		{
			if(originalS.at(posI) != '_')
			{
				translatedPositions.push_back(posI);
			}
		}

		assert(translatedPositions.size() == S.length());
		assert(PRGID_2_INTid.count(seqIt->first));
		std::string translationFile = directory + "/translation/" + Utilities::ItoStr(PRGID_2_INTid.at(seqIt->first)) + ".txt";

		std::ofstream translationStream;
		translationStream.open(translationFile.c_str());
		assert(translationStream.is_open());
		translationStream << Utilities::join(Utilities::ItoStr(translatedPositions), "\n") << "\n";
		translationStream.close();

		std::vector<std::string> sequences_fields;
		sequences_fields.resize(6);
		sequences_fields.at(0) = Utilities::ItoStr(PRGID_2_INTid.at(seqIt->first));
		sequences_fields.at(1) = seqIt->first;
		sequences_fields.at(2) = seqIt->first;

//		sequences_fields.at(3) = seqIt->first;
//		sequences_fields.at(4) = Utilities::ItoStr(1);
//		sequences_fields.at(5) = Utilities::ItoStr(S.length());

		sequencesStream << Utilities::join(sequences_fields, "\t") << "\n";
	}

	sequencesStream.close();

	std::string referenceGenomeFASTAPath = directory + "/referenceGenome/ref.fa";
	Utilities::writeFASTA(referenceGenomeFASTAPath, PRG_underlyingSequences_noGaps);
	Utilities::printToFile(directory + "/extendedReferenceGenomePath.txt", referenceGenomeFASTAPath);

	std::string PRGonly_referenceGenomeFASTAPath = directory + "/mapping_PRGonly/referenceGenome.fa";
	Utilities::writeFASTA(PRGonly_referenceGenomeFASTAPath, PRG_underlyingSequences_noGaps);


	mapper::bwa::BWAmapper bwaMapper(pF);
	bwaMapper.index(referenceGenomeFASTAPath);
}


std::vector<Edge*> simpleGraphSimulator::string_2_edgePath(std::string in)
{
	std::vector<Edge*> forReturn;
	for(unsigned int i = 0; i < in.length(); i++)
	{
		Edge* e = new Edge();
		e->emission = in.at(i);
		forReturn.push_back(e);
		generated_edges.insert(e);
	}
	return forReturn;
}

simpleGraphSimulator::~simpleGraphSimulator() {
	for(std::set<Edge*>::iterator edgeIt = generated_edges.begin(); edgeIt != generated_edges.end(); edgeIt++)
	{
		delete(*edgeIt);
	}
	delete(g);
}

void simpleGraphSimulator::init()
{
	assert(g == 0);

	std::string basicContig = Utilities::generateRandomSequence(graphLength);
	size_t PRG_mutations = 0;
	size_t PRG_gaps = 0;

	contigs_by_category["PRGScaffold"].push_back(string_2_edgePath(basicContig));
	for(unsigned int i = 0; i < forPRG_mutations; i++)
	{
		std::string modifiedContigSequence = basicContig;
		for(unsigned int cI = 0; cI < modifiedContigSequence.size(); cI++)
		{
			if(Utilities::randomDouble() <= mutationDensity)
			{
				if(Utilities::randomDouble() < 0.3)
				{
					if(!noGaps_switchContigs)
					{
						modifiedContigSequence.at(cI) = '_';
						PRG_gaps++;
					}
				}
				else
				{
					modifiedContigSequence.at(cI) = Utilities::randomNucleotide();
					PRG_mutations++;
				}
			}
		}
		contigs_by_category["PRGMutated"].push_back(string_2_edgePath(modifiedContigSequence));
	}

	/*
		for(unsigned int i = 0; i < additionalForReads_mutations; i++)
		{
			std::string modifiedContigSequence = basicContig;
			for(unsigned int cI = 0; cI < modifiedContigSequence.size(); cI++)
			{
				if(Utilities::randomDouble() <= mutationDensity)
				{
					if(Utilities::randomDouble() < 0.3)
					{
						if(!noGaps_additionalContigs)
						{
							modifiedContigSequence.at(cI) = '_';
							PRG_gaps++;
						}
					}
					else
					{
						modifiedContigSequence.at(cI) = Utilities::randomNucleotide();
						PRG_mutations++;
					}
				}
			}
			contigs_by_category["nonPRGMutated"].push_back(string_2_edgePath(modifiedContigSequence));
		}
	*/

	for(unsigned int i = 0; i < forPRG_largeGaps; i++)
	{
		std::string modifiedContigSequence = basicContig;
		for(unsigned int cI = 0; cI < modifiedContigSequence.size(); cI++)
		{
			if(Utilities::randomDouble() <= gapStartFrequency)
			{
				int gapLength = Utilities::randomPoisson(gapExpectedLength);
				if(gapLength > 0)
				{
					int lastPosition = cI + gapLength - 1;
					if(lastPosition > (modifiedContigSequence.length() - 1))
					{
						lastPosition = modifiedContigSequence.length() - 1;
					}
					for(unsigned int cII = cI; cII <= lastPosition; cII++)
					{
						if(!noGaps_switchContigs)
						{
							modifiedContigSequence.at(cII) = '_';
							PRG_gaps++;
						}
					}
					cI = lastPosition;
				}
			}
		}
		contigs_by_category["PRGLargeGaps"].push_back(string_2_edgePath(modifiedContigSequence));
	}

//	std::cout << "Construct PRG from sequences:\n";
	size_t haplotypeI = 0;
	std::vector<std::string> contigCategoriesForPRG = {"PRGScaffold", "PRGMutated", "PRGLargeGaps"};
	std::vector<std::string> loci;
	for(unsigned int categoryI = 0; categoryI < contigCategoriesForPRG.size(); categoryI++)
	{
		std::string category = contigCategoriesForPRG.at(categoryI);
		for(unsigned int contigI = 0; contigI < contigs_by_category.at(category).size(); contigI++)
		{
			std::vector<Edge*> contig_asEdges = contigs_by_category.at(category).at(contigI);
			std::string contigAsString;
			for(unsigned int posI = 0; posI < contig_asEdges.size(); posI++)
			{
				contigAsString.append(contig_asEdges.at(posI)->getEmission());
			}
			assert(contigAsString.length() == contig_asEdges.size());
			haplotypeI++;
			if(loci.size() == 0)
			{
				for(unsigned int locusI = 0; locusI < contigAsString.length(); locusI++)
				{
					std::string locusID = "L" + Utilities::ItoStr(locusI);
					loci.push_back(locusID);
				}
			}
			std::string PRGid = "PRG_" + Utilities::ItoStr(haplotypeI);
			hP.addString(loci, PRGid, contigAsString);
			std::cout << PRGid << "\t" << contigAsString << "\n";

			PRGID_2_internalID[PRGid] = make_pair(category, contigI);
			PRGID_2_INTid[PRGid] = haplotypeI;
		}
	}

	g = new Graph();
	g->buildFromHaplotypes(hP, false, 10);

	// std::cout << "Constructed graph from " << haplotypeI << " sequences.\n" << std::flush;
}

std::map<std::string, std::vector<mapper::reads::PRGContigBAMAlignment>> simpleGraphSimulator::simulateBAMAlignments(double contigCoverage, std::string qualityMatrixFile_, int read_length_, double insertSize_mean_, double insertSize_sd_, bool error)
{
	std::map<std::string, std::vector<mapper::reads::PRGContigBAMAlignment>> simulatedAlignments;

	auto oneSimulatedRead2PRGContigBAMAlignment = [error](std::vector<Edge*> contig_edgePath, simulator::oneRead R) -> mapper::reads::PRGContigBAMAlignment {
		mapper::reads::PRGContigBAMAlignment forReturn;
		forReturn.graph_aligned_levels = R.coordinates_edgePath;
		unsigned int sequenceCharacters = 0;

		for(unsigned int i = 0; i < R.coordinates_edgePath.size(); i++)
		{
			if(R.coordinates_edgePath.at(i) == -1)
			{
				forReturn.graph_aligned.push_back('_');
			}
			else
			{
				assert(R.coordinates_edgePath.at(i) < (int)contig_edgePath.size());
				std::string edgeEmission = contig_edgePath.at(R.coordinates_edgePath.at(i))->getEmission();
				assert(edgeEmission.length() == 1);
				forReturn.graph_aligned.push_back(edgeEmission.at(0));
			}

			if(R.sequence.at(i) != '_')
			{
				sequenceCharacters++;
			}
		}
		forReturn.sequence_aligned = R.sequence;
		forReturn.sequence_aligned_startInRaw = 0;
		forReturn.sequence_aligned_stopInRaw = sequenceCharacters - 1;
		forReturn.reverse = false;

		assert(forReturn.sequence_aligned.length() == forReturn.graph_aligned.length());

		if(! error)
		{
			//assert(forReturn.graph_aligned == forReturn.sequence_aligned);
		}

		/*
		std::cout << "oneSimulatedRead2PRGContigBAMAlignment" << "\n";
		std::cout << "\t" << "forReturn.graph_aligned_levels" << ": " << Utilities::join(Utilities::ItoStr(forReturn.graph_aligned_levels), ", ") << "\n";
		std::cout << "\t" << "forReturn.graph_aligned       " << ": " << forReturn.graph_aligned << "\n";
		std::cout << "\t" << "forReturn.sequence_aligned    " << ": " << forReturn.sequence_aligned << "\n";
		std::cout << "\n" << std::flush;
		*/
		return forReturn;
	};

	simulator::readSimulator rS(qualityMatrixFile_, read_length_);
	for(std::map<std::string, std::vector<std::vector<Edge*>>>::iterator categoryIt = contigs_by_category.begin(); categoryIt != contigs_by_category.end(); categoryIt++)
	{
		for(unsigned int contigI = 0; contigI < categoryIt->second.size(); contigI++)
		{
//			std::cout << "Simulate " << contigI << " " << categoryIt->first << "\n" << std::flush;
			std::vector<Edge*> contig_edgepath = categoryIt->second.at(contigI);
			int otherContigI = contigI + 1;
			if(otherContigI >= categoryIt->second.size())
			{
				otherContigI = contigI - 1;
			}
			if(otherContigI < 0)
			{
				otherContigI = contigI;
			}
			assert(otherContigI >= 0);
			assert(otherContigI < categoryIt->second.size());
			std::vector<Edge*> otherContig = categoryIt->second.at(otherContigI);

			std::vector<simulator::oneReadPair> simulatedReadPairs = rS.simulate_paired_reads_from_edgePath(contig_edgepath, contigCoverage, insertSize_mean_, insertSize_sd_, ! error, std::string(""), true);
//			std::cout << "a\n" << std::flush;
			for(unsigned int pairI = 0; pairI < simulatedReadPairs.size(); pairI++)
			{
				simulator::oneReadPair rP = simulatedReadPairs.at(pairI);
				if(rP.firstRead_minusStrand)
				{
					rP.reads.first.invert();
				}
				else
				{
					rP.reads.second.invert();
				}
				mapper::reads::PRGContigBAMAlignment BA1 = oneSimulatedRead2PRGContigBAMAlignment(noGaps_switchContigs ? otherContig : contig_edgepath, rP.reads.first);
				mapper::reads::PRGContigBAMAlignment BA2 = oneSimulatedRead2PRGContigBAMAlignment(noGaps_switchContigs ? otherContig : contig_edgepath, rP.reads.second);
				simulatedAlignments[categoryIt->first].push_back(BA1);
				simulatedAlignments[categoryIt->first].push_back(BA2);
			}
//			std::cout << "b\n" << std::flush;
		}
	}


	return simulatedAlignments;
}

