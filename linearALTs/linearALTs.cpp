/*
 * linearALTs.cpp
 *
 *  Created on: 14.01.2016
 *      Author: AlexanderDilthey
 */

#include "linearALTs.h"

#include "../mapper/bwa/BWAmapper.h"

#include "../Utilities.h"
#include <map>
#include <locale>
#include <vector>
#include <set>
#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include <stdexcept>
#include <exception>
#include "../Graph/Graph.h"

#include "../mapper/reads/protoSeeds.h"

double r_indels = 0.005;
double logP_matches = log((1-r_indels) * 0.99);
double logP_mismatches = log((1-r_indels) * 0.01 * (1.0/3.0));
double logP_insertions = log(r_indels*0.5);
double logP_deletions = log(r_indels*0.5);
double penalty_noPairing_P = 1e-6;

namespace linearALTs {

linearALTs::linearALTs(std::string directory, const pathFinder& pF) : pF(pF) {

	fastaFile_equalLengthHaplotypes = directory + "/equalLengthHaplotypesBlock/haplotypes.fa";
	assert(Utilities::fileExists(fastaFile_equalLengthHaplotypes));

	geneGraphDirectory = directory + "/genePRG";
	assert(Utilities::directoryExists(geneGraphDirectory));

	geneGraphFile = geneGraphDirectory + "/PRG/graph.txt";
	assert(Utilities::fileExists(geneGraphFile));

	std::string annotationsFile = directory + "/equalLengthHaplotypesBlock/haplotypes.annotation";
	assert(Utilities::fileExists(annotationsFile));

	std::string extendedReferenceGenomeFile = directory + "/extendedReferenceGenomePath.txt";
	assert(Utilities::fileExists(extendedReferenceGenomeFile));

	std::string extendedGenomeCoveredRegionsFile = directory + "/extendedGenome_coveredRegions.txt";
	assert(Utilities::fileExists(extendedGenomeCoveredRegionsFile));

	std::string equalLengthHaplotypesInfoFile = directory + "/equalLengthHaplotypesBlock/haplotypes_information.txt";
	assert(Utilities::fileExists(equalLengthHaplotypesInfoFile));

	fastaFile_explicitGenes = directory + "/regionalHaplotypesWithExplicitGenes/mappable.fa";
	assert(Utilities::fileExists(fastaFile_explicitGenes));

	std::string explicitGenesInformationFile = directory + "/regionalHaplotypesWithExplicitGenes/sequenceIDs.txt";
	assert(Utilities::fileExists(explicitGenesInformationFile));

	geneSequencesFile = directory + "/geneSequences/allSequences.mfa";
	assert(Utilities::fileExists(geneSequencesFile)); 

	simulateGeneSequences_alignedSequences = geneGraphDirectory + "/alignedTestSequences/sequences.mfa";
	assert(Utilities::fileExists(simulateGeneSequences_alignedSequences));


	equalLengthHaplotypes = Utilities::readFASTA(fastaFile_equalLengthHaplotypes);

	assert(equalLengthHaplotypes.size() > 1); // Sanity check
	for(std::map<std::string, std::string>::iterator ALTit = equalLengthHaplotypes.begin(); ALTit != equalLengthHaplotypes.end(); ALTit++)
	{
		std::string ALTid = ALTit->first;
		double propN = Utilities::proportionN(equalLengthHaplotypes.at(ALTid));
		if(propN > 0.2)
		{
			std::cerr << "Exclude because " << ALTid << " has " << propN << " Ns.\n";
			continue;
		}
		equalLengthHaplotypeIDs.push_back(ALTid);
	}

	mapper::bwa::BWAmapper bwaMapper(pF);
	bwaMapper.make_sure_ref_is_indexed(fastaFile_equalLengthHaplotypes);
	bwaMapper.make_sure_ref_is_indexed(fastaFile_explicitGenes);

	std::string line;

	// covered regions information
	std::ifstream extendedReferenceGenomeCoveredRegionsStream;
	extendedReferenceGenomeCoveredRegionsStream.open(extendedGenomeCoveredRegionsFile.c_str());
	assert(extendedReferenceGenomeCoveredRegionsStream.is_open());
	assert(extendedReferenceGenomeCoveredRegionsStream.good());
	std::getline(extendedReferenceGenomeCoveredRegionsStream, line);
	std::string headerLine = line;
	std::vector<std::string> headerFields = Utilities::split(headerLine, "\t");
	assert(headerFields.at(0) == "ContigID");
	assert(headerFields.at(1) == "AlternativeID");
	assert(headerFields.at(2) == "Start");
	assert(headerFields.at(3) == "Stop");
	while(extendedReferenceGenomeCoveredRegionsStream.good())
	{
		std::getline(extendedReferenceGenomeCoveredRegionsStream, line);
		Utilities::eraseNL(line);

		if(line.length() == 0)
		{
			continue;
		}

		std::vector<std::string> fields = Utilities::split(line, "\t");
		assert(fields.size() == headerFields.size());

		std::string extendedReferenceGenomeID = fields.at(0);
		std::string haplotypeID = fields.at(1);
		int begin_0based = fields.at(2).length() ? Utilities::StrtoI(fields.at(2)) : -1;
		int end_0based = fields.at(3).length() ? Utilities::StrtoI(fields.at(3)) : -1;
		assert(((begin_0based == -1) && (end_0based == -1)) || ((begin_0based != -1) && (end_0based != -1)));

		assert(haplotypeID_2_extendedReferenceGenomeID.count(haplotypeID) == 0);
		haplotypeID_2_extendedReferenceGenomeID[haplotypeID] = extendedReferenceGenomeID;

		extendedReferenceGenome_coveredRegions[extendedReferenceGenomeID] = std::make_pair(begin_0based, end_0based);
	}

	// equal length haplotypes information

	std::ifstream equalLengthHaplotypesInformationStream;
	equalLengthHaplotypesInformationStream.open(equalLengthHaplotypesInfoFile.c_str());
	assert(equalLengthHaplotypesInformationStream.is_open());
	assert(equalLengthHaplotypesInformationStream.good());
	std::getline(equalLengthHaplotypesInformationStream, line);
	headerLine = line;
	headerFields = Utilities::split(headerLine, "\t");
	assert(headerFields.at(0) == "HaplotypeID");
	assert(headerFields.at(1) == "ExtendedReferenceGenomeID");
	assert(headerFields.at(3) == "ExtendedReferenceGenomeID_Start0Based");
	assert(headerFields.at(4) == "ExtendedReferenceGenomeID_Stop0Based");

	while(equalLengthHaplotypesInformationStream.good())
	{
		std::getline(equalLengthHaplotypesInformationStream, line);
		Utilities::eraseNL(line);

		if(line.length() == 0)
		{
			continue;
		}

		std::vector<std::string> fields = Utilities::split(line, "\t");

		std::string haplotypeID = fields.at(0);
		std::string extendedReferenceGenomeID = fields.at(1);
		int begin_0based = Utilities::StrtoI(fields.at(3));
		int end_0based = Utilities::StrtoI(fields.at(4));
		assert(haplotypeID_2_extendedReferenceGenomeID.count(haplotypeID));
		equalLengthHaplotypes_inExtendedRef_0based[extendedReferenceGenomeID] = std::make_pair(begin_0based, end_0based);
	}

	// equal length haplotypes annotation
	std::ifstream annotationStream;
	annotationStream.open(annotationsFile.c_str());
	assert(annotationStream.is_open());
	assert(annotationStream.good());
	std::getline(annotationStream, line);
	Utilities::eraseNL(line);
	headerLine = line;
	headerFields = Utilities::split(headerLine, "\t");
	assert(headerFields.at(0) == "HaplotypeID");
	assert(headerFields.at(1) == "GeneID");
	assert(headerFields.at(3) == "Start_0based");
	assert(headerFields.at(4) == "Stop_0based");

	while(annotationStream.good())
	{
		std::getline(annotationStream, line);
		Utilities::eraseNL(line);

		if(line.length() == 0)
		{
			continue;
		}

		std::vector<std::string> fields = Utilities::split(line, "\t");

		std::string haplotypeID = fields.at(0);
		std::string gene = fields.at(1);
		int start = Utilities::StrtoI(fields.at(3));
		int stop = Utilities::StrtoI(fields.at(4));
		assert(start < stop);

		for(int i = start; i <= stop; i++)
		{
			equalLengthHaplotypes_annotations_perPosition[haplotypeID][i].insert(gene);
		}
	}

	// explicit genes sequences information
	std::ifstream explicitGenesInformationStream;
	explicitGenesInformationStream.open(explicitGenesInformationFile.c_str());
	assert(explicitGenesInformationStream.is_open());
	assert(explicitGenesInformationStream.good());
	std::getline(explicitGenesInformationStream, line);
	Utilities::eraseNL(line);
	headerLine = line;
	headerFields = Utilities::split(headerLine, "\t");
	assert(headerFields.at(0) == "OriginalID");
	assert(headerFields.at(1) == "FastaID");
	assert(headerFields.at(2) == "IsGene");

	while(explicitGenesInformationStream.good())
	{
		std::getline(explicitGenesInformationStream, line);
		Utilities::eraseNL(line);

		if(line.length() == 0)
		{
			continue;
		}

		std::vector<std::string> fields = Utilities::split(line, "\t");

		std::string haplotypeID = fields.at(0);
		std::string fastaID = fields.at(1);
		bool isGene = Utilities::StrtoI(fields.at(2));

		assert(explicitGenes_haplotypeID_2_fastaID.count(haplotypeID) == 0);
		explicitGenes_haplotypeID_2_fastaID[haplotypeID] = fastaID;

		if(isGene)
		{
			explicitGenes_fastaIDs_genes.insert(fastaID);
		}
	}

	for(std::map<std::string, std::map<int, std::set<std::string>>>::iterator annotationIt = equalLengthHaplotypes_annotations_perPosition.begin(); annotationIt != equalLengthHaplotypes_annotations_perPosition.end(); annotationIt++)
	{
		std::string haplotypeID = annotationIt->first;
		for(std::map<int, std::set<std::string>>::iterator positionIt = annotationIt->second.begin(); positionIt != annotationIt->second.end(); positionIt++)
		{
			int position = positionIt->first;
			int count = positionIt->second.size();
			if(count > 1)
			{
				std::cout << "Overlap warning: " << haplotypeID << ":" << position << " " << count << " annotations.\n";
			}
		}
	}

	assert(equalLengthHaplotypes.size() == haplotypeID_2_extendedReferenceGenomeID.size());
	for(std::map<std::string, std::string>::iterator ALTit = equalLengthHaplotypes.begin(); ALTit != equalLengthHaplotypes.end(); ALTit++)
	{
		std::string ALTid = ALTit->first;
		if(equalLengthHaplotypes_annotations_perPosition.count(ALTid) == 0)
		{
			std::cerr << "No annotations for " << ALTid << " -- correct?\n" << std::flush;
			equalLengthHaplotypes_annotations_perPosition[ALTid] = std::map<int, std::set<std::string>>();
		}
		// assert(haplotypeID_2_extendedReferenceGenomeID.count(ALTid));
		assert(explicitGenes_haplotypeID_2_fastaID.count(ALTid));
	}
}

double linearALTs::scoreAlignment(const BamTools::BamAlignment& alignment)
{
	int insertions = 0;
	int deletions = 0;
	int matches = 0;
	int mismatches = 0;
	int hard_padding = 0;
	int soft_padding = 0;
	for(BamTools::CigarOp op : alignment.CigarData)
	{
		switch(op.Type) {
		case 'H':
			insertions += op.Length;
			hard_padding += op.Length;
			break;
		case 'S':
			insertions += op.Length;
			soft_padding += op.Length;
			break;
		case 'I':
			insertions += op.Length;
			break;
		case 'D':
			deletions += op.Length;
			break;
		case 'M':
			matches += op.Length;
			break;
		default:
			std::cerr << "CIGAR type " << op.Type << " not implemented!\n" << std::flush;
			assert( 1 == 0 );
		};
	}


	if(alignment.Length)
	{
		int expected_read_length = (insertions + matches - hard_padding);
		if(!(expected_read_length == alignment.Length))
		{
			mapper::reads::protoSeeds::printAlignmentInfo(0, alignment);
			std::cerr << "expected_read_length" << ": " << expected_read_length << "\n";
			std::cerr << "alignment.Length" << ": " << alignment.Length << "\n";
			std::cerr << std::flush;
		}
		assert(expected_read_length == alignment.Length);
	}

	u_int32_t edit_distance;
	assert(alignment.GetTag("NM", edit_distance));

	if(!((int)edit_distance >= ((insertions - hard_padding - soft_padding) + deletions)))
	{
		mapper::reads::protoSeeds::printAlignmentInfo(0, alignment);
		std::cerr << "edit_distance" << ": " << edit_distance << "\n";
		std::cerr << "insertions" << ": " << insertions << "\n";
		std::cerr << "deletions" << ": " << deletions << "\n";
		std::cerr << std::flush;
	}
	assert((int)edit_distance >= ((insertions - hard_padding - soft_padding) + deletions));

	mismatches = edit_distance - ((insertions - hard_padding - soft_padding) + deletions);
	if(!(mismatches >= 0))
	{
		mapper::reads::protoSeeds::printAlignmentInfo(0, alignment);
		std::cerr << "mismatches" << ": " << mismatches << "\n";
		std::cerr << "edit_distance" << ": " << edit_distance << "\n";
		std::cerr << "insertions" << ": " << insertions << "\n";
		std::cerr << "deletions" << ": " << deletions << "\n";
		std::cerr << std::flush;
	}
	assert(mismatches >= 0);
	assert(mismatches <= matches);

	matches = matches - mismatches;
	assert(matches >= 0);

	double log_return = (
		logP_matches * matches +
		logP_mismatches * mismatches +
		logP_insertions * insertions +
		logP_deletions * deletions
	);

	return log_return;
}

std::string linearALTs::CIGARasString(const BamTools::BamAlignment& al)
{
	std::string forReturn;
	for(auto CIGARop : al.CigarData)
	{
		forReturn.push_back(CIGARop.Type);
		forReturn += Utilities::ItoStr(CIGARop.Length);
	}
	return forReturn;
}

std::pair<std::string, std::string> linearALTs::simulateEqualLengthDiploidHaplotypes(std::string outputDirectory, const simulator::simulator& S, std::map<std::string, std::set<std::string>>& read2gene_R1, std::map<std::string, std::set<std::string>>& read2gene_R2, double haploidCoverage, bool withError)
{

	read2gene_R1.clear();
	read2gene_R2.clear();

	std::string h1_id = Utilities::choose_uniformly_from_vector(equalLengthHaplotypeIDs);
	std::string h2_id = Utilities::choose_uniformly_from_vector(equalLengthHaplotypeIDs);

	std::cout << "Underlying haplotypes:\n";
	std::cout << h1_id << "\n";
	std::cout << h2_id << "\n";
	std::cout << "\n";


	S.simulateFromDiploidHaplotypes(equalLengthHaplotypes.at(h1_id), equalLengthHaplotypes.at(h2_id), outputDirectory, haploidCoverage, withError);

	std::string levels_1 = outputDirectory + "/R_1.levels";
	std::string levels_2 = outputDirectory + "/R_2.levels";

	auto readLevels = [&](std::string file, std::map<std::string, std::set<std::string>>& read2gene)
	{
		std::ifstream levelsStream;
		levelsStream.open(file.c_str());
		assert(levelsStream.is_open());

		std::string line;
		while(levelsStream.good())
		{
			std::getline(levelsStream, line);
			Utilities::eraseNL(line);

			if(line.length() == 0)
			{
				continue;
			}

			std::string readID = line;
			if(!(readID.substr(0, 1) == "@"))
			{
				std::cerr << "This line is supposed to be a read ID, but it is not!\n" << line << "\n" << std::flush;
			}

			assert(readID.substr(0, 1) == "@");
			readID = readID.substr(1);
			assert(levelsStream.good());

			if(!((readID.substr(0, 3) == "H2r") || (readID.substr(0, 3) == "H1r")))
			{
				std::cerr << "Weird read ID: " << readID << "\n" << std::flush;
			}
			assert((readID.substr(0, 3) == "H2r") || (readID.substr(0, 3) == "H1r"));
			std::string haplotypeID = ( (readID.substr(0, 3) == "H2r") ? h2_id : h1_id );
			assert(equalLengthHaplotypes_annotations_perPosition.count(haplotypeID));

			std::getline(levelsStream, line);
			std::string levelsLine = line;
			std::vector<std::string> levelsFields = Utilities::split(levelsLine, " ");

			std::vector<int> levels = Utilities::StrtoI(levelsFields);

			read2gene[readID] = std::set<std::string>();
			for(auto level : levels)
			{
				if(equalLengthHaplotypes_annotations_perPosition.at(haplotypeID).count(level))
				{
					read2gene[readID].insert(equalLengthHaplotypes_annotations_perPosition.at(haplotypeID).at(level).begin(), equalLengthHaplotypes_annotations_perPosition.at(haplotypeID).at(level).end());
					// std::cerr << "Have gene " << *(ALTs_annotations_perPosition.at(haplotypeID).at(level).begin()) << " for " << readID << "\n";
				}
			}
		}
	};

	readLevels(levels_1, read2gene_R1);
	readLevels(levels_2, read2gene_R2);

	std::string genes_1 = outputDirectory + "/R_1.genes";
	std::string genes_2 = outputDirectory + "/R_2.genes";

	auto printGenes = [](std::string output, std::map<std::string, std::set<std::string>> read2gene) {
		std::ofstream outputStream;
		outputStream.open(output.c_str());
		assert(outputStream.is_open());
		for(auto readData : read2gene)
		{
			outputStream << readData.first << "\n";
			outputStream << Utilities::join(std::vector<std::string>(readData.second.begin(), readData.second.end()), "") << "\n";
		}
		outputStream.close();
	};

	printGenes(genes_1, read2gene_R1);
	printGenes(genes_2, read2gene_R2);

	std::string fastq_1 = outputDirectory + "/R_1.fq";
	std::string fastq_2 = outputDirectory + "/R_2.fq";

	std::string outputBAM = outputDirectory + "/reads.bam";
	mapper::bwa::BWAmapper bwaMapper(pF);
	bwaMapper.map(fastaFile_equalLengthHaplotypes, fastq_1, fastq_2, outputBAM);
	std::cout << "Generated " << outputBAM << "\n";

	std::string outputBAM_1 = outputDirectory + "/reads_1.bam";
	std::string outputBAM_2 = outputDirectory + "/reads_2.bam";

	bwaMapper.map_all_unpaired_unsorted(fastaFile_equalLengthHaplotypes, fastq_1,  outputBAM_1);
	std::cout << "Generated " << outputBAM_1 << "\n";

	bwaMapper.map_all_unpaired_unsorted(fastaFile_equalLengthHaplotypes, fastq_2,  outputBAM_2);
	std::cout << "Generated " << outputBAM_2 << "\n";

	return make_pair(h1_id, h2_id);

}

std::string linearALTs::getGeneGraphPRGDirectory()
{
	return geneGraphDirectory;
}

void linearALTs::simulateGeneSequences(std::string outputDirectory, const simulator::simulator& S, bool simulateFromGraph, double haploidCoverage, bool withError)
{
	Utilities::make_or_clearDirectory(outputDirectory);

	if(simulateFromGraph)
	{
		Graph* g = new Graph();
		g->readFromFile(geneGraphFile);
		S.simulateFromGraph(g, 1, outputDirectory, haploidCoverage, withError);
		delete(g);
	}
	else
	{
		std::map<std::string, std::string> alignedTestSequences = Utilities::readFASTA(simulateGeneSequences_alignedSequences, true);
		assert(alignedTestSequences.size());
		std::vector<std::string> available_kir_gene_sequences;
		for(auto seq : alignedTestSequences)
		{
			std::string alleleSequence = seq.first;
			// for(size_t i = 0; i < alleleSequence.length(); i++)
			// {
				// alleleSequence.at(i) = std::toupper(alleleSequence.at(i));
				// if(alleleSequence.at(i) == '-')
				// {
					// alleleSequence.at(i) = '_';
				// }
			// }
			available_kir_gene_sequences.push_back(alleleSequence);
		}

		std::string seqID = Utilities::choose_uniformly_from_vector(available_kir_gene_sequences);

		std::cout << "Simulating one indivdual from graph one KIR gene reads with sequence ID: " << seqID << "\n";

		S.simulateFromHaploidAlignedSequence(alignedTestSequences.at(seqID), outputDirectory, haploidCoverage, withError);
	}
}

std::pair<std::string, std::string> linearALTs::haplotypeLikelihoods(std::string outputDirectory, double insertSize_mean, double insertSize_sd, int fallBackReadLength, std::pair<std::string, std::string> printPair)
{
//	std::cout << "read2gene_R1.size()" << ": " << read2gene_R1.size() << "\n";
//	std::cout << "read2gene_R2.size()" << ": " << read2gene_R2.size() << "\n";

	bool verbose = false;
	if(verbose)
	{
		std::cout << Utilities::timestamp() << "Enter linearALTs::haplotypeLikelihoods(..)." << std::flush;
	}

	std::string outputBAM_1 = outputDirectory + "/reads_1.bam";
	std::string outputBAM_2 = outputDirectory + "/reads_2.bam";
	assert(Utilities::fileExists(outputBAM_1));
	assert(Utilities::fileExists(outputBAM_2));

	std::map<std::string, std::vector<double>> likelihoods_perRead_perHaplotype;
	boost::math::normal rnd_InsertSize(insertSize_mean, insertSize_sd);

	BamTools::BamReader R1;
	R1.Open(outputBAM_1);

	BamTools::BamReader R2;
	R2.Open(outputBAM_2);

	const BamTools::RefVector references_R1 = R1.GetReferenceData();
	const BamTools::RefVector references_R2 = R2.GetReferenceData();
	assert(references_R1.size() == references_R2.size());
	std::vector<std::string> references_IDs;
	std::vector<size_t> references_numericalIndices;
	std::vector<size_t> reference_lengths;
	std::map<std::string, size_t> referenceID_2_numericalIdx;
	for(size_t i = 0; i < references_R1.size(); i++)
	{
		assert(references_R1.at(i).RefLength == references_R2.at(i).RefLength);
		assert(references_R1.at(i).RefName == references_R2.at(i).RefName);
		references_IDs.push_back(references_R1.at(i).RefName);
		reference_lengths.push_back(references_R1.at(i).RefLength);
		referenceID_2_numericalIdx[references_R1.at(i).RefName] = i;
	}


	BamTools::BamAlignment currentAlignment1;
	BamTools::BamAlignment currentAlignment2;

	std::map<std::string, std::vector<BamTools::BamAlignment>> collectedAlignments1;
	std::map<std::string, std::vector<BamTools::BamAlignment>> collectedAlignments2;

	std::map<std::string, bool> finishedCollection1;
	std::map<std::string, bool> finishedCollection2;


	std::map<std::string, std::vector<double>> loglikelihoods_perRead_perReference;
	std::map<std::string, std::pair<std::map<std::string, double>, std::map<std::string, double>>> hitGenes_perRead;

	std::string currentReadID1;
	std::string currentReadID2;
	size_t total_r1 = 0;
	size_t total_r2 = 0;
	size_t processed_r1 = 0;
	size_t processed_r2 = 0;
	size_t read_at_once = 1000;
	while(1)
	{
		size_t read_from_1 = 0;
		while(R1.GetNextAlignment(currentAlignment1))
		{
			total_r1++;
			read_from_1++;
			if(currentReadID1.length())
			{
				assert(collectedAlignments1.count(currentReadID1));
			}
			collectedAlignments1[currentAlignment1.Name].push_back(currentAlignment1);

			if(currentReadID1 != currentAlignment1.Name)
			{
				if(currentReadID1.length())
				{
					finishedCollection1.at(currentReadID1) = true;
				}
				finishedCollection1[currentAlignment1.Name] = false;
				currentReadID1 = currentAlignment1.Name;
			}

			if(read_from_1 >= read_at_once)
			{
				break;
			}
		}

		size_t read_from_2 = 0;
		while(R2.GetNextAlignment(currentAlignment2))
		{
			total_r2++;
			read_from_2++;
			if(currentReadID2.length())
			{
				assert(collectedAlignments2.count(currentReadID2));
			}
			collectedAlignments2[currentAlignment2.Name].push_back(currentAlignment2);

			if(currentReadID2 != currentAlignment2.Name)
			{
				if(currentReadID2.length())
				{
					finishedCollection2.at(currentReadID2) = true;
				}
				finishedCollection2[currentAlignment2.Name] = false;
				currentReadID2 = currentAlignment2.Name;
			}

			if(read_from_2 >= read_at_once)
			{
				break;
			}
		}

		std::set<std::string> delete_readID;
		for(std::map<std::string, std::vector<BamTools::BamAlignment>>::iterator readIDit = collectedAlignments1.begin(); readIDit != collectedAlignments1.end(); readIDit++)
		{
			std::string readID = readIDit->first;
			if(collectedAlignments2.count(readID))
			{
				if(finishedCollection1.at(readID) && finishedCollection2.at(readID))
				{
					assert(loglikelihoods_perRead_perReference.count(readID) == 0);
					std::pair<std::map<std::string, double>, std::map<std::string, double>> hitGenes;

					loglikelihoods_perRead_perReference[readID] = processCollectedAlignments(collectedAlignments1.at(readID), collectedAlignments2.at(readID), hitGenes, processed_r1, processed_r2, references_IDs, {}, rnd_InsertSize, fallBackReadLength, verbose);
					hitGenes_perRead[readID] = hitGenes;

					delete_readID.insert(readID);
				}
			}
		}
		for(auto readID : delete_readID)
		{
			collectedAlignments1.erase(readID);
			collectedAlignments2.erase(readID);
		}

		if((read_from_1 == 0) && (read_from_2 == 0))
		{
			break;
		}
	}

	std::set<std::string> processed_readID;
	for(std::map<std::string, std::vector<BamTools::BamAlignment>>::iterator readIDit = collectedAlignments1.begin(); readIDit != collectedAlignments1.end(); readIDit++)
	{
		std::string readID = readIDit->first;
		if(collectedAlignments2.count(readID))
		{
			std::pair<std::map<std::string, double>, std::map<std::string, double>> hitGenes;

			assert(loglikelihoods_perRead_perReference.count(readID) == 0);
			loglikelihoods_perRead_perReference[readID] = processCollectedAlignments(collectedAlignments1.at(readID), collectedAlignments2.at(readID), hitGenes, processed_r1, processed_r2, references_IDs, {}, rnd_InsertSize, fallBackReadLength, verbose);
			hitGenes_perRead[readID] = hitGenes;

			processed_readID.insert(readID);
		}
		else
		{
			std::pair<std::map<std::string, double>, std::map<std::string, double>> hitGenes;

			assert(loglikelihoods_perRead_perReference.count(readID) == 0);
			loglikelihoods_perRead_perReference[readID] = processCollectedAlignments(collectedAlignments1.at(readID), {}, hitGenes, processed_r1, processed_r2, references_IDs, {}, rnd_InsertSize, fallBackReadLength, verbose);
			hitGenes_perRead[readID] = hitGenes;

			processed_readID.insert(readID);
		}
	}

	for(std::map<std::string, std::vector<BamTools::BamAlignment>>::iterator readIDit = collectedAlignments2.begin(); readIDit != collectedAlignments2.end(); readIDit++)
	{
		std::string readID = readIDit->first;
		if(processed_readID.count(readID) == 0)
		{
			std::pair<std::map<std::string, double>, std::map<std::string, double>> hitGenes;

			assert(collectedAlignments1.count(readID) == 0);
			assert(loglikelihoods_perRead_perReference.count(readID) == 0);
			loglikelihoods_perRead_perReference[readID] = processCollectedAlignments({}, collectedAlignments2.at(readID), hitGenes, processed_r1, processed_r2, references_IDs, {}, rnd_InsertSize, fallBackReadLength, verbose);
			hitGenes_perRead[readID] = hitGenes;

		}
	}

	assert(total_r1 == processed_r1);
	assert(total_r2 == processed_r2);

	if(verbose)
	{
		std::cout << Utilities::timestamp() << "linearALTs::haplotypeLikelihoods(..): Determine most likely pair of haplotypes...\n";
	}

	assert(references_IDs.size() > 1);
	std::vector<double> haplotype_pair_scores;
	std::vector<std::pair<size_t, size_t>> haplotype_pair_indices;

	for(size_t referenceI1 = 0; referenceI1 < references_IDs.size(); referenceI1++)
	{
		for(size_t referenceI2 = referenceI1; referenceI2 < references_IDs.size(); referenceI2++)
		{
			std::string referenceID1 = references_IDs.at(referenceI1);
			std::string referenceID2 = references_IDs.at(referenceI2);

			double pair_log_likelihood = 0;
			for(auto readLoglikelihoods : loglikelihoods_perRead_perReference)
			{
				std::string readID = readLoglikelihoods.first;

				double read_loglikelihood_h1 = readLoglikelihoods.second.at(referenceI1);
				double read_loglikelihood_h2 = readLoglikelihoods.second.at(referenceI2);

				if(readID == "H1p|||18919|||120|||69")
				{
					// std::cout << "Pair calculation for " << readID << " on " << referenceID1 << " and " << referenceID2 << ", likelihoods " << read_loglikelihood_h1 << "\t" << read_loglikelihood_h2 << "\n";
				}

				// debug
				double p_h1 = exp(read_loglikelihood_h1);
				double p_h2 = exp(read_loglikelihood_h2);
				assert((p_h1 >= 0) && (p_h1 <= 1));
				assert((p_h2 >= 0) && (p_h2 <= 1));

				double read_loglikelihood_avg = Utilities::logAvg(read_loglikelihood_h1, read_loglikelihood_h2);

				pair_log_likelihood += read_loglikelihood_avg;
			}

			if(verbose)
			{
				std::cout << "\tPair " << referenceI1 << " / " << referenceI2 << ": " << pair_log_likelihood << " (" << loglikelihoods_perRead_perReference.size() << " reads).\n";
			}

			haplotype_pair_scores.push_back(pair_log_likelihood);
			haplotype_pair_indices.push_back(make_pair(referenceI1, referenceI2));
		}
	}

	std::pair<double, unsigned int> maxHaplotypePairI = Utilities::findVectorMax(haplotype_pair_scores);
	std::cout << Utilities::timestamp() << "\tdone. (One) maximum pair is " << maxHaplotypePairI.second << " with LL = " << maxHaplotypePairI.first << "\n" << std::flush;

	std::vector<double> haplotype_pair_scores_normalized;
	double LL_max = maxHaplotypePairI.first;
	double P_sum = 0;
	for(unsigned int cI = 0; cI < haplotype_pair_scores.size(); cI++)
	{
		double LL = haplotype_pair_scores.at(cI);
		double P = exp(LL - LL_max);
		P_sum += P;
	}
	if(P_sum > 0)
	{
		for(unsigned int cI = 0; cI < haplotype_pair_scores.size(); cI++)
		{
			double LL = haplotype_pair_scores.at(cI);
			double P = exp(LL - LL_max);
			double P_normalized = P / P_sum;
			if(!((P_normalized >= 0) && (P_normalized <= 1)))
			{
				std::cerr << "P_normalized: " << P_normalized << "\n";
				std::cerr << "P: " << P << "\n";
				std::cerr << "LL: " << LL << "\n";
				std::cerr << "LL_max: " << LL_max << "\n";

				std::cerr << std::flush;
			}
			assert(P_normalized >= 0);
			assert(P_normalized <= 1);
			haplotype_pair_scores_normalized.push_back(P_normalized);
		}
	}
	else
	{
		for(unsigned int cI = 0; cI < haplotype_pair_scores.size(); cI++)
		{
			haplotype_pair_scores_normalized.push_back(1.0/(double)haplotype_pair_scores.size());
		}
	}

	std::set<std::string> calledPairs;
	std::cout << "NORMALIZED probability of maximum pair: " << haplotype_pair_scores_normalized.at(maxHaplotypePairI.second) << "\n";
	std::cout << "\t" << references_IDs.at(haplotype_pair_indices.at(maxHaplotypePairI.second).first) << "\n";
	std::cout << "\t" << references_IDs.at(haplotype_pair_indices.at(maxHaplotypePairI.second).second) << "\n";
	std::cout << std::flush;
	calledPairs.insert(references_IDs.at(haplotype_pair_indices.at(maxHaplotypePairI.second).first));
	calledPairs.insert(references_IDs.at(haplotype_pair_indices.at(maxHaplotypePairI.second).second));

	if(printPair.first.length())
	{
		assert(referenceID_2_numericalIdx.count(printPair.first));
		assert(referenceID_2_numericalIdx.count(printPair.second));

		for(size_t pairI = 0; pairI < haplotype_pair_indices.size(); pairI++)
		{
			std::pair<size_t, size_t> pair_indices = haplotype_pair_indices.at(pairI);
			std::string h1 = references_IDs.at(pair_indices.first);
			std::string h2 = references_IDs.at(pair_indices.second);

			if( ((h1 == printPair.first) && (h2 == printPair.second)) ||
				((h2 == printPair.first) && (h1 == printPair.second))	)
			{
				std::cout << "PRINTPAIR probability " << haplotype_pair_scores_normalized.at(pairI) << "\n";
				std::cout << "\t" << references_IDs.at(haplotype_pair_indices.at(pairI).first) << "\n";
				std::cout << "\t" << references_IDs.at(haplotype_pair_indices.at(pairI).second) << "\n";
				std::cout << std::flush;
			}
		}

		/*
		std::cout << "Print per-read probabilities for two print pairs.\n";
		std::cout << "\t\t" << printPair.first << "\t" << printPair.second << "\n";
		for(auto readLoglikelihoods : loglikelihoods_perRead_perReference)
		{
			std::vector<double> likelihoods_per_reference = readLoglikelihoods.second;
			std::vector<double> likelihoods_per_reference_normalized = Utilities::normalize_log_vector(likelihoods_per_reference);

			std::cout << "\t\t" << readLoglikelihoods.first << ": " << likelihoods_per_reference_normalized.at(referenceID_2_numericalIdx.at(printPair.first)) << "\t" << likelihoods_per_reference_normalized.at(referenceID_2_numericalIdx.a(printPair.second)) << "\n";
			// std::cout << "\t\t\t" << "#0: " << likelihoods_per_reference_normalized.at(0) << "\n";
		}
		*/
	}

	return make_pair(references_IDs.at(haplotype_pair_indices.at(maxHaplotypePairI.second).first), references_IDs.at(haplotype_pair_indices.at(maxHaplotypePairI.second).second));
}


void linearALTs::reads2Genes(std::string outputDirectory, double insertSize_mean, double insertSize_sd, int fallBackReadLength, std::map<std::string, std::set<std::string>>& true_read2gene_R1, std::map<std::string, std::set<std::string>>& true_read2gene_R2, std::map<std::string, std::map<std::string, int>>& ret_reads_trueGene_2_inferredGene, std::map<std::string, std::map<std::string, int>>& ret_reads_inferredGene_2_trueGene, const std::set<std::string>& restrictToHaplotypes, double T)
{
	bool verbose = false;
	if(verbose)
	{
		std::cout << Utilities::timestamp() << "Enter linearALTs::reads2Genes(..)." << std::flush;
	}

	std::string outputBAM_1 = outputDirectory + "/reads_1.bam";
	std::string outputBAM_2 = outputDirectory + "/reads_2.bam";
	assert(Utilities::fileExists(outputBAM_1));
	assert(Utilities::fileExists(outputBAM_2));

	std::map<std::string, std::vector<double>> likelihoods_perRead_perHaplotype;
	boost::math::normal rnd_InsertSize(insertSize_mean, insertSize_sd);

	BamTools::BamReader R1;
	R1.Open(outputBAM_1);

	BamTools::BamReader R2;
	R2.Open(outputBAM_2);

	const BamTools::RefVector references_R1 = R1.GetReferenceData();
	const BamTools::RefVector references_R2 = R2.GetReferenceData();
	assert(references_R1.size() == references_R2.size());
	std::vector<std::string> references_IDs;
	std::vector<size_t> references_numericalIndices;
	std::vector<size_t> reference_lengths;
	std::map<std::string, size_t> referenceID_2_numericalIdx;
	for(size_t i = 0; i < references_R1.size(); i++)
	{
		assert(references_R1.at(i).RefLength == references_R2.at(i).RefLength);
		assert(references_R1.at(i).RefName == references_R2.at(i).RefName);
		references_IDs.push_back(references_R1.at(i).RefName);
		reference_lengths.push_back(references_R1.at(i).RefLength);
		referenceID_2_numericalIdx[references_R1.at(i).RefName] = i;
	}


	BamTools::BamAlignment currentAlignment1;
	BamTools::BamAlignment currentAlignment2;

	std::map<std::string, std::vector<BamTools::BamAlignment>> collectedAlignments1;
	std::map<std::string, std::vector<BamTools::BamAlignment>> collectedAlignments2;

	std::map<std::string, bool> finishedCollection1;
	std::map<std::string, bool> finishedCollection2;

	size_t processed_r1 = 0;
	size_t processed_r2 = 0;
	std::map<std::string, std::vector<double>> loglikelihoods_perRead_perReference;
	std::map<std::string, std::pair<std::map<std::string, double>, std::map<std::string, double>>> hitGenes_perRead;

	std::string currentReadID1;
	std::string currentReadID2;
	size_t total_r1 = 0;
	size_t total_r2 = 0;
	size_t read_at_once = 1000;
	while(1)
	{
		size_t read_from_1 = 0;
		while(R1.GetNextAlignment(currentAlignment1))
		{
			total_r1++;
			read_from_1++;
			if(currentReadID1.length())
			{
				assert(collectedAlignments1.count(currentReadID1));
			}
			collectedAlignments1[currentAlignment1.Name].push_back(currentAlignment1);

			if(currentReadID1 != currentAlignment1.Name)
			{
				if(currentReadID1.length())
				{
					finishedCollection1.at(currentReadID1) = true;
				}
				finishedCollection1[currentAlignment1.Name] = false;
				currentReadID1 = currentAlignment1.Name;
			}

			if(read_from_1 >= read_at_once)
			{
				break;
			}
		}

		size_t read_from_2 = 0;
		while(R2.GetNextAlignment(currentAlignment2))
		{
			total_r2++;
			read_from_2++;
			if(currentReadID2.length())
			{
				assert(collectedAlignments2.count(currentReadID2));
			}
			collectedAlignments2[currentAlignment2.Name].push_back(currentAlignment2);

			if(currentReadID2 != currentAlignment2.Name)
			{
				if(currentReadID2.length())
				{
					finishedCollection2.at(currentReadID2) = true;
				}
				finishedCollection2[currentAlignment2.Name] = false;
				currentReadID2 = currentAlignment2.Name;
			}

			if(read_from_2 >= read_at_once)
			{
				break;
			}
		}

		std::set<std::string> delete_readID;
		for(std::map<std::string, std::vector<BamTools::BamAlignment>>::iterator readIDit = collectedAlignments1.begin(); readIDit != collectedAlignments1.end(); readIDit++)
		{
			std::string readID = readIDit->first;
			if(collectedAlignments2.count(readID))
			{
				if(finishedCollection1.at(readID) && finishedCollection2.at(readID))
				{
					assert(loglikelihoods_perRead_perReference.count(readID) == 0);
					std::pair<std::map<std::string, double>, std::map<std::string, double>> hitGenes;

					loglikelihoods_perRead_perReference[readID] = processCollectedAlignments(collectedAlignments1.at(readID), collectedAlignments2.at(readID), hitGenes, processed_r1, processed_r2, references_IDs, restrictToHaplotypes, rnd_InsertSize, fallBackReadLength, verbose);
					hitGenes_perRead[readID] = hitGenes;

					delete_readID.insert(readID);
				}
			}
		}
		for(auto readID : delete_readID)
		{
			collectedAlignments1.erase(readID);
			collectedAlignments2.erase(readID);
		}

		if((read_from_1 == 0) && (read_from_2 == 0))
		{
			break;
		}
	}

	std::set<std::string> processed_readID;
	for(std::map<std::string, std::vector<BamTools::BamAlignment>>::iterator readIDit = collectedAlignments1.begin(); readIDit != collectedAlignments1.end(); readIDit++)
	{
		std::string readID = readIDit->first;
		if(collectedAlignments2.count(readID))
		{
			std::pair<std::map<std::string, double>, std::map<std::string, double>> hitGenes;

			assert(loglikelihoods_perRead_perReference.count(readID) == 0);
			loglikelihoods_perRead_perReference[readID] = processCollectedAlignments(collectedAlignments1.at(readID), collectedAlignments2.at(readID), hitGenes, processed_r1, processed_r2, references_IDs, restrictToHaplotypes, rnd_InsertSize, fallBackReadLength, verbose);
			hitGenes_perRead[readID] = hitGenes;

			processed_readID.insert(readID);
		}
		else
		{
			std::pair<std::map<std::string, double>, std::map<std::string, double>> hitGenes;

			assert(loglikelihoods_perRead_perReference.count(readID) == 0);
			loglikelihoods_perRead_perReference[readID] = processCollectedAlignments(collectedAlignments1.at(readID), {}, hitGenes, processed_r1, processed_r2, references_IDs, restrictToHaplotypes, rnd_InsertSize, fallBackReadLength, verbose);
			hitGenes_perRead[readID] = hitGenes;

			processed_readID.insert(readID);
		}
	}

	for(std::map<std::string, std::vector<BamTools::BamAlignment>>::iterator readIDit = collectedAlignments2.begin(); readIDit != collectedAlignments2.end(); readIDit++)
	{
		std::string readID = readIDit->first;
		if(processed_readID.count(readID) == 0)
		{
			std::pair<std::map<std::string, double>, std::map<std::string, double>> hitGenes;

			assert(collectedAlignments1.count(readID) == 0);
			assert(loglikelihoods_perRead_perReference.count(readID) == 0);
			loglikelihoods_perRead_perReference[readID] = processCollectedAlignments({}, collectedAlignments2.at(readID), hitGenes, processed_r1, processed_r2, references_IDs, restrictToHaplotypes, rnd_InsertSize, fallBackReadLength, verbose);
			hitGenes_perRead[readID] = hitGenes;

		}
	}

	assert(total_r1 == processed_r1);
	assert(total_r2 == processed_r2);

	if(verbose)
	{
		std::cout << Utilities::timestamp() << "linearALTs::haplotypeLikelihoods(..): Analyse per-read gene hits...\n";
	}

	for(std::map<std::string, std::pair<std::map<std::string, double>, std::map<std::string, double>>>::iterator readIt = hitGenes_perRead.begin(); readIt != hitGenes_perRead.end(); readIt++)
	{
		bool haveInferredGenes = ((readIt->second.first.size() > 0) || (readIt->second.second.size() > 0));

		std::string readID = readIt->first;

		std::set<std::string> trueGenes_R1 = true_read2gene_R1[readID];
		std::set<std::string> trueGenes_R2 = true_read2gene_R2[readID];

		auto add_to_read_gene_statistics = [&](const std::map<std::string, double>& hitGenes, const std::set<std::string>& trueGenes)
		{
			std::string use_inferred_gene;
			std::string use_true_gene;

			if(hitGenes.size() > 0)
			{
				std::pair<double, std::string> maxGene_PP = Utilities::findStringMapMax(hitGenes);

				if(maxGene_PP.first >= T)
				{
					use_inferred_gene = maxGene_PP.second;
				}
				else
				{
					use_inferred_gene = "NONE";
				}
			}
			else
			{
				use_inferred_gene = "NONE";
			}

			if(trueGenes.size() == 0)
			{
				use_true_gene = "NONE";
			}
			else if(trueGenes.size() == 1)
			{
				use_true_gene = *(trueGenes.begin());
			}
			else
			{
				if(trueGenes.count(use_inferred_gene))
				{
					use_true_gene = use_inferred_gene;
				}
				else
				{
					use_true_gene = *(trueGenes.begin());
				}
			}


			if(ret_reads_trueGene_2_inferredGene[use_true_gene].count(use_inferred_gene) == 0)
			{
				ret_reads_trueGene_2_inferredGene[use_true_gene][use_inferred_gene] = 0;
			}
			ret_reads_trueGene_2_inferredGene[use_true_gene][use_inferred_gene]++;


			if(ret_reads_inferredGene_2_trueGene[use_inferred_gene].count(use_true_gene) == 0)
			{
				ret_reads_inferredGene_2_trueGene[use_inferred_gene][use_true_gene] = 0;
			}
			ret_reads_inferredGene_2_trueGene[use_inferred_gene][use_true_gene]++;
		};

		add_to_read_gene_statistics(readIt->second.first, trueGenes_R1);
		add_to_read_gene_statistics(readIt->second.second, trueGenes_R2);

		// if(haveInferredGenes || trueGenes_R1.size() || trueGenes_R2.size())
		// if(trueGenes_R1.size() || trueGenes_R2.size())
		if(0)
		{
			std::cout << "Analyze " << readID << "\n";
			std::cout << "\t" << "Read 1:" << "\n";
			std::cout << "\t\tTrue genes: " << Utilities::join(std::vector<std::string>(trueGenes_R1.begin(), trueGenes_R1.end()), ",") << "\n";
			for(std::map<std::string, double>::iterator geneIt = readIt->second.first.begin(); geneIt != readIt->second.first.end(); geneIt++)
			{
				std::cout << "\t\t\t" << geneIt->first << " PP: " << geneIt->second;
				if(trueGenes_R1.count(geneIt->first))
				{
					std::cout << "\t!";
				}
				std::cout << "\n";
			}
			std::cout << "\t" << "Read 2:" << "\n";
			std::cout << "\t\tTrue genes: " << Utilities::join(std::vector<std::string>(trueGenes_R2.begin(), trueGenes_R2.end()), ",") << "\n";
			for(std::map<std::string, double>::iterator geneIt = readIt->second.second.begin(); geneIt != readIt->second.second.end(); geneIt++)
			{
				std::cout << "\t\t\t" << geneIt->first << " PP: " << geneIt->second;
				if(trueGenes_R2.count(geneIt->first))
				{
					std::cout << "\t!";
				}
				std::cout << "\n";
			}
		}
	}
}


bool linearALTs::alignedReadPair_strandsValid(const BamTools::BamAlignment& alignment_read1, const BamTools::BamAlignment& alignment_read2)
{

	int r1_firstPos = alignment_read1.Position;
	int r1_lastPos = alignment_read1.GetEndPosition();

	int r2_firstPos = alignment_read2.Position;
	int r2_lastPos = alignment_read2.GetEndPosition();


	if (alignment_read1.IsReverseStrand() != alignment_read2.IsReverseStrand())
	{
		if(! alignment_read1.IsReverseStrand())
		{
			if(r1_firstPos < r2_firstPos)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			if(r1_lastPos > r2_lastPos)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}
	else
	{
		return false;
	}
}

int linearALTs::alignedReadPair_distance(const BamTools::BamAlignment& alignment_read1, const BamTools::BamAlignment& alignment_read2)
{
	int r1_firstPos = alignment_read1.Position;
	int r1_lastPos = alignment_read1.GetEndPosition();

	int r2_firstPos = alignment_read2.Position;
	int r2_lastPos = alignment_read2.GetEndPosition();

	if(r1_firstPos < r2_firstPos)
	{
		int D = (r2_firstPos - r1_lastPos - 1);
		return D;

	}
	else
	{
		assert(r1_firstPos >= r2_firstPos);
		int D = (r1_firstPos - r2_lastPos - 1);
		return D;
	}
}

std::vector<double> linearALTs::processCollectedAlignments(const std::vector<BamTools::BamAlignment>& Alignments_R1, const std::vector<BamTools::BamAlignment>& Alignments_R2, std::pair<std::map<std::string, double>, std::map<std::string, double>>& hitGenes, size_t& processed_r1, size_t& processed_r2, const std::vector<std::string>& references_IDs, const std::set<std::string>& restrictToHaplotypesForGenes, const boost::math::normal& rnd_InsertSize, int fallBackReadLength, bool verbose)
{
	if(verbose)
	{
		std::cout << "Call processCollectedAlignments with \n";
		std::cout << "\tR1: " << Alignments_R1.size() << " alignments.\n";
		std::cout << "\tR2: " << Alignments_R2.size() << " alignments.\n";
	}

	std::string readID;
	std::set<int32_t> references;

//		std::cout << "Read " << "\n";

	std::vector<int> read_lengths_R1;
	std::vector<int> read_lengths_R2;
	std::map<int32_t, std::vector<BamTools::BamAlignment>> Alignments_R1_byReferenceID;
	std::map<int32_t, std::vector<BamTools::BamAlignment>> Alignments_R2_byReferenceID;

//		std::cout << "\tR1" << "\n";
	if(Alignments_R1.size())
	{
		std::string rID = Alignments_R1.at(0).Name;
		for(size_t alignmentI = 0; alignmentI < Alignments_R1.size(); alignmentI++)
		{
			auto oneAlignment = Alignments_R1.at(alignmentI);
			assert(oneAlignment.Name == rID);
			references.insert(oneAlignment.RefID);
			Alignments_R1_byReferenceID[oneAlignment.RefID].push_back(oneAlignment);
			read_lengths_R1.push_back(oneAlignment.Length);

//				if(rID == "H1p|||18919|||120|||69")
//				{
//					std::cout << readID << "\t" << references_IDs.at(oneAlignment.RefID) << "\n" << std::flush;
//				}
//
//				if(oneAlignment.RefID == referenceID_2_numericalIdx.at(printPair.first))
//				{
//					std::cout << "\t\t" <<  printPair.first << "\t" << CIGARasString(oneAlignment) << "\n";
//				}
//				if(oneAlignment.RefID == referenceID_2_numericalIdx.at(printPair.second))
//				{
//					std::cout << "\t\t" <<  printPair.second << "\t" << CIGARasString(oneAlignment) << "\n";
//				}
		}
		readID = rID;
	}

//		std::cout << "\tR2" << "\n";
	if(Alignments_R2.size())
	{
		std::string rID = Alignments_R2.at(0).Name;
		for(size_t alignmentI = 0; alignmentI < Alignments_R2.size(); alignmentI++)
		{
			auto oneAlignment = Alignments_R2.at(alignmentI);
			assert(oneAlignment.Name == rID);
			references.insert(oneAlignment.RefID);
			Alignments_R2_byReferenceID[oneAlignment.RefID].push_back(oneAlignment);
			read_lengths_R2.push_back(oneAlignment.Length);

//				if(oneAlignment.RefID == referenceID_2_numericalIdx.at(printPair.first))
//				{
//					std::cout << "\t\t" <<  printPair.first << "\t" << CIGARasString(oneAlignment) << "\n";
//				}
//				if(oneAlignment.RefID == referenceID_2_numericalIdx.at(printPair.second))
//				{
//					std::cout << "\t\t" <<  printPair.second << "\t" << CIGARasString(oneAlignment) << "\n";
//				}
		}
		readID = rID;
	}

//		if(readID == "H1p|||18919|||120|||69")
//		{
////			assert(1 == 0);
//		}
//		std::cout << "\t" << readID << "\n";

	auto findReadLength = [&](const std::vector<int>& read_lengths) -> int
	{
		int inferredL = -1;
		for(int l : read_lengths)
		{
			if(l != 0)
			{
				if(inferredL == -1)
				{
					inferredL = l;
				}
				else
				{
					if(inferredL != l)
					{
						inferredL = -1;
						break;
					}
				}
			}
		}

		if(inferredL == -1)
		{
			inferredL = fallBackReadLength;
		}

		assert(inferredL > 50); // sanity check

		return inferredL;
	};

	int readLength_R1 = findReadLength(read_lengths_R1);
	int readLength_R2 = findReadLength(read_lengths_R2);

	processed_r1 += Alignments_R1.size();
	processed_r2 += Alignments_R2.size();

	if(Alignments_R1.size() && Alignments_R2.size())
	{
		assert(Alignments_R1.at(0).Name == Alignments_R2.at(0).Name);
	}

	if(verbose)
	{
		if(readID.length())
		{
			std::cout << "\t" << references.size() << " combined contigs.\n";
		}
	}

	if(verbose)
	{
		std::cout << "\t" << "Read 1 length" << ": " << readLength_R1 << "\n";
		std::cout << "\t" << "Read 2 length" << ": " << readLength_R2 << "\n";
		std::cout << "\tNow calculating per-contig likelihoods\n";
	}

	// verbose = (readID == "H1p|||18919|||120|||69");
	std::vector<std::pair<std::set<std::string>, std::set<std::string>>> alignmentCombination_hitGenes;
	std::vector<double> alignmentCombination_hitGenes_likelihoods;

	std::vector<double> likelihoods_per_reference;
	for(int referenceI = 0; referenceI < (int)references_IDs.size(); referenceI++)
	{
		if(verbose)
			std::cout << "\t\t" << "Reference ID" << ": " << references_IDs.at(referenceI) << "\n";

		assert(equalLengthHaplotypes_annotations_perPosition.count(references_IDs.at(referenceI)));

		std::vector<double> scores;
		std::vector<std::pair<size_t, size_t>> scores_pairs;

		std::vector<BamTools::BamAlignment> Alignments_R1 = Alignments_R1_byReferenceID[referenceI];
		std::vector<BamTools::BamAlignment> Alignments_R2 = Alignments_R2_byReferenceID[referenceI];

		for(size_t alignmentI_R1 = 0; alignmentI_R1 < (Alignments_R1.size() + 1); alignmentI_R1++)
		{
			double score_alignment_R1 = (alignmentI_R1 < Alignments_R1.size()) ? scoreAlignment(Alignments_R1.at(alignmentI_R1)) : (readLength_R1 * logP_insertions);

			std::set<std::string> genes_R1;
			if(alignmentI_R1 < Alignments_R1.size())
			{
				int alignment_start = Alignments_R1.at(alignmentI_R1).Position;
				int alignment_stop =  Alignments_R1.at(alignmentI_R1).GetEndPosition();
				assert(alignment_start < alignment_stop);
				for(int alignmentPos = alignment_start; alignmentPos <= alignment_stop; alignmentPos++)
				{
					// std::cout << "Check " << alignmentPos << " on " << references_IDs.at(referenceI) << "\n";
					if(equalLengthHaplotypes_annotations_perPosition.at(references_IDs.at(referenceI)).count(alignmentPos))
					{
						genes_R1.insert(equalLengthHaplotypes_annotations_perPosition.at(references_IDs.at(referenceI)).at(alignmentPos).begin(), equalLengthHaplotypes_annotations_perPosition.at(references_IDs.at(referenceI)).at(alignmentPos).end());
					}
				}
			}

			for(size_t alignmentI_R2 = 0; alignmentI_R2 < (Alignments_R2.size() + 1); alignmentI_R2++)
			{
				double score_alignment_R2 = (alignmentI_R2 < Alignments_R2.size()) ? scoreAlignment(Alignments_R2.at(alignmentI_R2)) : (readLength_R2 * logP_insertions);

				std::set<std::string> genes_R2;
				if(alignmentI_R2 < Alignments_R2.size())
				{
					int alignment_start = Alignments_R2.at(alignmentI_R2).Position;
					int alignment_stop =  Alignments_R2.at(alignmentI_R2).GetEndPosition();
					assert(alignment_start < alignment_stop);
					for(int alignmentPos = alignment_start; alignmentPos <= alignment_stop; alignmentPos++)
					{
						if(equalLengthHaplotypes_annotations_perPosition.at(references_IDs.at(referenceI)).count(alignmentPos))
						{
							genes_R2.insert(equalLengthHaplotypes_annotations_perPosition.at(references_IDs.at(referenceI)).at(alignmentPos).begin(), equalLengthHaplotypes_annotations_perPosition.at(references_IDs.at(referenceI)).at(alignmentPos).end());
						}
					}
				}
				double score_pairs = 0;

				if(verbose)
				{
					std::cout << "\t\t\t" << "Pair " << alignmentI_R1 << " / " << alignmentI_R2 << "\n";
				}

				if((alignmentI_R1 == Alignments_R1.size()) || (alignmentI_R2 == Alignments_R2.size()))
				{
					score_pairs = log(penalty_noPairing_P);
					if(verbose)
					{
						std::cout << "\t\t\t\tOnly one read aligned.\n";
					}
				}
				else
				{
					if(verbose)
					{
						std::cout << "\t\t\t\talignedReadPair_strandsValid: " << alignedReadPair_strandsValid(Alignments_R1.at(alignmentI_R1), Alignments_R2.at(alignmentI_R2)) << " alignedReadPair_distance: " << alignedReadPair_distance(Alignments_R1.at(alignmentI_R1), Alignments_R2.at(alignmentI_R2)) << "\n";
					}
					if(alignedReadPair_strandsValid(Alignments_R1.at(alignmentI_R1), Alignments_R2.at(alignmentI_R2)))
					{
						int distance = alignedReadPair_distance(Alignments_R1.at(alignmentI_R1), Alignments_R2.at(alignmentI_R2));
						double distance_P = boost::math::pdf(rnd_InsertSize, distance);

						if(!(distance_P > 0))
						{
							distance_P = penalty_noPairing_P;
						}

						if(!((distance_P > 0) && (distance_P <= 1)))
						{
							std::cerr << "penalty_noPairing_P: " << penalty_noPairing_P << "\n";
							std::cerr << "distance_P: " << distance_P << "\n" << std::flush;
						}
						assert((distance_P > 0) && (distance_P <= 1));

						score_pairs = log(distance_P);
					}
					else
					{
						score_pairs = log(penalty_noPairing_P);
					}
				}

				double combined_score = score_alignment_R1 + score_alignment_R2 + score_pairs;
				scores.push_back(combined_score);
				scores_pairs.push_back(make_pair(alignmentI_R1, alignmentI_R2));

				if(verbose)
				{
					std::cout << "\t\t\t\t" << combined_score << " (" << score_alignment_R1 << "\t" << score_alignment_R2 << "\t" << score_pairs << ")\n";
				}


				if((restrictToHaplotypesForGenes.size() == 0) || (restrictToHaplotypesForGenes.count(references_IDs.at(referenceI))))
				{
					alignmentCombination_hitGenes_likelihoods.push_back(combined_score);
					alignmentCombination_hitGenes.push_back(make_pair(genes_R1, genes_R2));
				}
			}
		}

		std::pair<double, unsigned int> maximumScore = Utilities::findVectorMax(scores);

		if(verbose)
		{
			std::cout << "\n" << "\t\t\t Maximum pair with log likelihood " << maximumScore.first << "\n\n";
		}
		likelihoods_per_reference.push_back(maximumScore.first);
	}

	assert(likelihoods_per_reference.size() == references_IDs.size());

	double P_sum = 0;
	alignmentCombination_hitGenes_likelihoods = Utilities::normalize_log_vector(alignmentCombination_hitGenes_likelihoods);
	for(size_t alignmentI = 0; alignmentI < alignmentCombination_hitGenes_likelihoods.size(); alignmentI++)
	{
		double P = alignmentCombination_hitGenes_likelihoods.at(alignmentI);
		assert(P >= 0);
		assert(P <= 1);

		const std::set<std::string>& hitGenes_R1 = alignmentCombination_hitGenes.at(alignmentI).first;
		const std::set<std::string>& hitGenes_R2 = alignmentCombination_hitGenes.at(alignmentI).second;

		for(auto geneName : hitGenes_R1)
		{
			if(hitGenes.first.count(geneName) == 0)
			{
				hitGenes.first[geneName] = 0;
			}
			hitGenes.first.at(geneName) += P;
		}

		for(auto geneName : hitGenes_R2)
		{
			if(hitGenes.second.count(geneName) == 0)
			{
				hitGenes.second[geneName] = 0;
			}
			hitGenes.second.at(geneName) += P;
		}

		P_sum += P;
	}
	assert(abs(1 - P_sum) <= 1e-3);

	verbose = false;

	return likelihoods_per_reference;
}

linearALTs::~linearALTs() {
	// TODO Auto-generated destructor stub
}

void linearALTs::extractReadsFromBAM(std::string outputDirectory, std::string BAM, const std::map<std::string, std::pair<int, int>>& regions, bool removePairNumbers)
{
	BamTools::BamReader reader;
	reader.Open(BAM);

	reader.LocateIndex();
    if ( ! reader.HasIndex() )
    {
		throw std::runtime_error("File "+BAM+" does not seem to be indexed - please specify indexed BAM!");
    }

    class rawRead {
    public:
    	std::string readID;
    	std::string sequence;
    	std::string qualities;
    };


    class fastq_readPair;
    class fastq_readPair {
    public:
    	bool have1;
    	bool have2;

    	rawRead a1;
    	rawRead a2;

    	fastq_readPair() : have1(false), have2(false)
    	{

    	}

    	bool takeAlignment(rawRead a, int which)
    	{
    		bool success = false;
    		assert((which == 1) || (which == 2));
    		if(which == 1)
    		{
    			if(have1 == false)
    			{
    				success = true;
    				a1 = a;
    				have1 = true;
    			}
    		}
    		else if(which == 2)
    		{
    			if(have2 == false)
    			{
    				success = true;
    				a2 = a;
    				have2 = true;
    			}
    		}
    		return success;
    	}

    	bool isComplete()
    	{
    		return (have1 && have2);
    	}

    	void printFASTQ(std::ofstream& f1, std::ofstream& f2, bool removePairNumbers)
    	{
    		assert(a1.sequence.length() == a1.qualities.length());
    		assert(a2.sequence.length() == a2.qualities.length());

    		f1 << "@" << a1.readID << (removePairNumbers ? "" : "/1") << "\n";
    		f1 << a1.sequence << "\n";
    		f1 << "+" << "\n";
    		f1 << a1.qualities << "\n";

    		f2 << "@" << a2.readID << (removePairNumbers ? "" : "/2") << "\n";
    		f2 << a2.sequence << "\n";
    		f2 << "+" << "\n";
    		f2 << a2.qualities << "\n";
    	}
    };

    std::map<std::string, fastq_readPair> readsStore;

	// std::cout << "Extraction " << regions.size() << " regions.\n" << std::flush;
	size_t examined_reads = 0;
    for(std::map<std::string, std::pair<int, int>>::const_iterator regionIt = regions.begin(); regionIt != regions.end(); regionIt++)
    {
    	std::string regionName = regionIt->first;
    	int regionStartSpec = regionIt->second.first;
    	int regionStopSpec = regionIt->second.second;

    	int id_for_name = reader.GetReferenceID(regionName);
    	if(id_for_name == -1)
    	{
    		throw std::runtime_error("File "+BAM+" does not seem to have contig "+ regionName);
    	}
    	assert(id_for_name != -1);
    	BamTools::BamRegion stretch_region_BAMTools;
    	stretch_region_BAMTools.LeftRefID = id_for_name;
    	stretch_region_BAMTools.RightRefID = id_for_name;
		int regionStartExtraction = regionStartSpec;
		int regionStopExtraction = regionStopSpec;
    	if(regionStartExtraction == -1)
    	{
			assert(regionStopExtraction == -1);
    		regionStartExtraction = 0;
    		regionStopExtraction = reader.GetReferenceData().at(id_for_name).RefLength;
    	}
		assert(regionStartExtraction >= 0);
		assert(regionStopExtraction >= 0);
		assert(regionStartExtraction < regionStopExtraction);
		
		stretch_region_BAMTools.LeftPosition = regionStartExtraction;
		stretch_region_BAMTools.RightPosition = regionStopExtraction;		

		
		// std::cout << "\t" << stretch_region_BAMTools.LeftRefID << " " << stretch_region_BAMTools.RightRefID << "\n";
		// std::cout << "\t\t" << stretch_region_BAMTools.LeftPosition << " - " << stretch_region_BAMTools.RightPosition << "\n" << std::flush;
    	
		assert(reader.SetRegion(stretch_region_BAMTools));

		BamTools::BamAlignment alignment;
		while(reader.GetNextAlignment(alignment))
		{
			examined_reads++;
			
			int startAlignment = alignment.Position;
			int stopAlignment = alignment.GetEndPosition();

			if( ( (startAlignment >= regionStartExtraction) && (startAlignment <= regionStopExtraction) ) &&
				( (stopAlignment >= regionStartExtraction) && (stopAlignment <= regionStopExtraction) )
			)
			{
				if(! alignment.IsPrimaryAlignment())
				{
					continue;
				}
				if(! alignment.IsPaired())
				{
					continue;
				}

				std::string readName = alignment.Name;
				std::string qualities = alignment.Qualities;
				std::string sequence  = alignment.QueryBases;
				if ( alignment.IsReverseStrand() ) {
					std::reverse(qualities.begin(), qualities.end());
					sequence = Utilities::seq_reverse_complement(sequence);
				}

				rawRead thisRead;
				thisRead.readID = readName;
				thisRead.sequence = sequence;
				thisRead.qualities = qualities;
				int whichMate =  (alignment.IsFirstMate()) ? 1 : 2;

				if(readsStore.count(readName) == 0)
				{
					fastq_readPair p;
					bool success = p.takeAlignment(thisRead, whichMate);
					assert(success);
					readsStore[readName] = p;
				}
				else
				{
					fastq_readPair& thisPair = readsStore.at(readName);
					bool success = thisPair.takeAlignment(thisRead, whichMate);
					if(! success)
					{
						std::cerr << "There is a problem with the read IDs in this BAM.\n";
						std::cerr << "Read ID: " << readName << "\n";
						std::cerr << "whichMate: " << whichMate << "\n";
						std::cerr << "thisPair.have1: " << thisPair.have1 << " with ID " << thisPair.a1.readID << "\n";
						std::cerr << "thisPair.have2: " << thisPair.have2 << " with ID " << thisPair.a2.readID << "\n" << std::flush;
					}
					assert(success);

				}
			}
		}
    }
	
	std::string fn_1 = outputDirectory  + "/R_1.fq";
	std::string fn_2 = outputDirectory  + "/R_2.fq";

	std::ofstream fastq_1_output;
	fastq_1_output.open(fn_1.c_str());
	if(! fastq_1_output.is_open())
	{
		throw std::runtime_error("readFilter::doFilter(): Cannot open file "+fn_1);
	}

	std::ofstream fastq_2_output;
	fastq_2_output.open(fn_2.c_str());
	if(! fastq_2_output.is_open())
	{
		throw std::runtime_error("readFilter::doFilter(): Cannot open file "+fn_2);
	}

	size_t printed_reads = 0;
    for(std::map<std::string, fastq_readPair>::iterator storedReadIt = readsStore.begin(); storedReadIt != readsStore.end(); storedReadIt++)
    {
    	if(storedReadIt->second.isComplete())
    	{
			printed_reads++;
    		storedReadIt->second.printFASTQ(fastq_1_output, fastq_2_output, removePairNumbers);
    	}
    }
}


void linearALTs::extractReads_equalLengthHaplotypes(std::string outputDirectory, std::string BAM)
{
	extractReadsFromBAM(outputDirectory, BAM, equalLengthHaplotypes_inExtendedRef_0based);

	std::string fastq_1 = outputDirectory + "/R_1.fq";
	std::string fastq_2 = outputDirectory + "/R_2.fq";
	
	mapper::bwa::BWAmapper bwaMapper(pF);
	
	// std::string outputBAM = outputDirectory + "/reads.bam";	
	// bwaMapper.map(fastaFile_equalLengthHaplotypes, fastq_1, fastq_2, outputBAM);
	// std::cout << "Generated " << outputBAM << "\n";

	std::string outputBAM_1 = outputDirectory + "/reads_1.bam";
	std::string outputBAM_2 = outputDirectory + "/reads_2.bam";

	bwaMapper.map_all_unpaired_unsorted(fastaFile_equalLengthHaplotypes, fastq_1,  outputBAM_1);

	bwaMapper.map_all_unpaired_unsorted(fastaFile_equalLengthHaplotypes, fastq_2,  outputBAM_2);
	
}

void linearALTs::extractReads_extendedReferenceGenome(std::string outputDirectory, std::string BAM)
{
	// std::cout << "extendedReferenceGenome_coveredRegions.size()" << extendedReferenceGenome_coveredRegions.size() << "\n";
	// for(auto oneS : extendedReferenceGenome_coveredRegions)
	// {
		// std::cout << "\t" << oneS.first << " " << oneS.second.first << ": " << (oneS.second.second - oneS.second.first) << "\n";
	// }
	// std::cout << "equalLengthHaplotypes_inExtendedRef_0based.size()" << equalLengthHaplotypes_inExtendedRef_0based.size() << "\n";
	// for(auto oneS : equalLengthHaplotypes_inExtendedRef_0based)
	// {
		// std::cout << "\t" << oneS.first << " " << oneS.second.first << ": " << (oneS.second.second - oneS.second.first) << "\n";
	// }	
	// std::cout << "\n" << std::flush;
	
	// assert(1 == 0);
	
	extractReadsFromBAM(outputDirectory, BAM, extendedReferenceGenome_coveredRegions);
}

void linearALTs::extractReads_geneGraph(std::string outputDirectory, std::string BAM)
{
	std::string outputDirectory_reads_overlapping_extendedReferenceGenome = outputDirectory + "/reads_coveringExtendedRegions";
	Utilities::make_or_clearDirectory(outputDirectory_reads_overlapping_extendedReferenceGenome);

	extractReadsFromBAM(outputDirectory_reads_overlapping_extendedReferenceGenome, BAM, extendedReferenceGenome_coveredRegions);
	std::string reads_coveringExtendedRegions_1 = outputDirectory_reads_overlapping_extendedReferenceGenome + "/R_1.fq";
	std::string reads_coveringExtendedRegions_2 = outputDirectory_reads_overlapping_extendedReferenceGenome + "/R_2.fq";
	assert(Utilities::fileExists(reads_coveringExtendedRegions_1));
	assert(Utilities::fileExists(reads_coveringExtendedRegions_2));

	std::string BAM_explicitGenes = outputDirectory_reads_overlapping_extendedReferenceGenome + "/mappedToExplicitGenes.bam";
	mapper::bwa::BWAmapper bwaMapper(pF);
	bwaMapper.map(fastaFile_explicitGenes, reads_coveringExtendedRegions_1, reads_coveringExtendedRegions_2, BAM_explicitGenes);
	assert(Utilities::fileExists(BAM_explicitGenes));

	std::map<std::string, std::pair<int, int>> map_genesExtraction;
	for(std::string fastaID : explicitGenes_fastaIDs_genes)
	{
		map_genesExtraction[fastaID] = make_pair(-1, -1);
	}

	extractReadsFromBAM(outputDirectory, BAM_explicitGenes, map_genesExtraction);
}

} /* namespace linearALTs */
