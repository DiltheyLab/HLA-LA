//============================================================================
// Name        : HLA-PRG-LA.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>

#include <map>
#include <assert.h>
#include <string>
#include <vector>
#include <exception>
#include <stdexcept>
#include <chrono>
#include <cstdio>
#include <fstream>

#include "mapper/processBAM.h"
#include "mapper/reads/PRGContigBAMAlignment.h"
#include "mapper/reads/verboseSeedChain.h"
#include "mapper/aligner/extensionAligner.h"
#include "mapper/bwa/BWAmapper.h"
#include "mapper/bowtie2/Bowtie2mapper.h"

#include "simulator/trueReadLevels.h"

#include "Graph/Graph.h"
#include "Graph/graphSimulator/simpleGraphSimulator.h"
#include "simulator/simulator.h"
#include "Graph/GraphAndEdgeIndex.h"

#include "fullLengthHMM/fullLengthHMM.h"

#include "Utilities.h"
#include "pathFinder.h"

#include "hla/HLATyper.h"

#include "linearALTs/linearALTs.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>

int main(int argc, char *argv[]) {

	std::vector<std::string> ARG (argv + 1, argv + argc + !argc);
	std::map<std::string, std::string> arguments;

	arguments["action"] = "readHMM";

	/*
	arguments["action"] = "PRGmapping";
	arguments["action"] = "testChainExtension";
	arguments["action"] = "testAlignments2Chains";
	arguments["action"] = "testPRGMapping";
	arguments["action"] = "testPRGMappingUnpaired";
	*/

	// arguments["action"] = "testRealBAM";

	// arguments["action"] = "prepareGraph";
	// arguments["action"] = "TestHLATyping";

	// arguments["action"] = "testCheckPresence";
	// arguments["action"] = "checkKIRgraph"; 

	// arguments["action"] = "KIR";

	// arguments["PRG_graph_dir"] = "/gpfs1/well/gsk_hla/HLA-PRG-LA/graphs/PRG_MHC_GRCh38_withIMGT";

	for(unsigned int i = 0; i < ARG.size(); i++)
	{
		if((ARG.at(i).length() > 2) && (ARG.at(i).substr(0, 2) == "--"))
		{
			std::string argname = ARG.at(i).substr(2);
			std::string argvalue = ARG.at(i+1);
			arguments[argname] = argvalue;
		}
	}


	/*
	
	std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>> trueHLA_test;
	std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>> inferredHLA_test;
	
	hla::HLATyper::read_inferred_types("NA12892", inferredHLA_test, "/Net/birch/data/dilthey/MHC-PRG/tmp/hla/C_Platinum_NA12892/R1_bestguess.txt");
	hla::HLATyper::read_true_types(trueHLA_test, "/Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended");
	hla::HLATyper::evaluate_HLA_types(trueHLA_test, inferredHLA_test);
	assert(2 == 5);
	*/

	// some overlap tests
	assert(! Utilities::intervalsOverlap(1, 10, 11, 20));
	assert(! Utilities::intervalsOverlap(5, 11, 1, 4));
	assert(Utilities::intervalsOverlap(5, 11, 8, 11));
	assert(Utilities::intervalsOverlap(8, 11, 1, 9));
	assert(Utilities::intervalsOverlap(8, 11, 1, 9));
	assert(Utilities::intervalsOverlap(8, 11, 9, 10));
	assert(Utilities::intervalsOverlap(9, 10, 8, 11));
	assert(Utilities::intervalsOverlap(1, 10, 2, 3));
	assert(Utilities::intervalsOverlap(2, 3, 1, 10));

	if(arguments.count("action") == 0)
	{
		std::cerr << "\n\nMissing --action parameter. Please don't try calling me directly; use HLA-LA.pl instead (see documentation on GitHub).\n" << std::endl;
		throw std::runtime_error("Missing arguments -- see above.");
	}
	
	std::set<std::string> noBinariesRequired = {"prepareGraph", "testBinary", "readHMM"};
	
	if((noBinariesRequired.count(arguments.at("action"))  == 0) && (! arguments.count("bwa_bin")))
	{
		throw std::runtime_error("Please specify arguments --bwa_bin");
	}
	if((noBinariesRequired.count(arguments.at("action"))  == 0) && (! arguments.count("samtools_bin")))
	{
		throw std::runtime_error("Please specify arguments --samtools_bin");
	}

	if(! std::system(NULL))
	{
		std::cerr << "\n\nMissing shell - std::system(NULL) has returned a 0 value.\n" << std::endl;
		throw std::runtime_error("Missing shell");		
	}

	pathFinder pF(arguments);
	assert(arguments.count("action"));
	if(arguments.at("action") == "testBinary")
	{
		std::cout << "\nHLA*LA binary functional!\n\n";
	}
	else if(arguments.at("action") == "readHMM")
	{
		if(arguments.count("inputPrefix") == 0)
		{	
			arguments["inputPrefix"] = Utilities::fileExists("../test/myTest.readAlleles") ? "../test/myTest" :  "C:\\Users\\Alexa\\Documents\\workspace\\HLA-PRG-LA\\test\\myTest";
		}
		
		std::map<std::string, unsigned int> gene_length;
		std::map<std::string, std::map<unsigned int, std::set<std::string>>> iteration_2_readIDs;
		std::map<std::string, unsigned int> readID_2_iteration;
		std::map<std::string, std::set<unsigned int>> gene2Iterations;

		std::map<std::string, std::string> reads_2_genes;
		std::map<std::string, std::map<std::string, std::pair<unsigned int, unsigned int>>> read_start_stop_positions;
		std::map<std::string, std::map<unsigned int, std::set<std::string>>> read_start_per_position;
		std::map<std::string, std::map<unsigned int, std::set<std::string>>> read_stop_per_position;
		std::map<std::string, std::map<std::string, std::map<unsigned int, std::string>>> read_genotypes_per_position;
		std::map<std::string, std::map<unsigned int, std::set<std::string>>> activeAlleles_per_position;
		std::map<std::string, std::map<std::string, std::string>> MSA_reference_sequences;
		std::map<std::string, std::map<std::string, std::string>> MSA_reference_sequences_whichHap;

		bool manualTest = false;
		if(manualTest)
		{
			std::cout << "Check0" << "\n" << std::flush;

			arguments.at("inputPrefix") = "../test/manualTest";

			std::string gene = "TEST";

			gene_length[gene] = 20;

			MSA_reference_sequences[gene][gene+"*01"] = "AAATAAAAATAAAATAAAAA";
			MSA_reference_sequences_whichHap[gene][gene+"*01"] = "1";

			MSA_reference_sequences[gene][gene+"*02"] = "AAACAAAAACAAAACAAAAA";
			MSA_reference_sequences_whichHap[gene][gene+"*02"] = "2";

			activeAlleles_per_position[gene][3].insert("T");
			activeAlleles_per_position[gene][3].insert("C");

			activeAlleles_per_position[gene][9].insert("T");
			activeAlleles_per_position[gene][9].insert("C");

			activeAlleles_per_position[gene][14].insert("T");
			activeAlleles_per_position[gene][14].insert("C");

			for(unsigned int pI = 0; pI < MSA_reference_sequences.at(gene).at(gene+"*01").length(); pI++)
			{
				if(activeAlleles_per_position[gene].count(pI) == 0)
				{
					activeAlleles_per_position[gene][pI].insert(MSA_reference_sequences.at(gene).at(gene+"*01").substr(pI, 1));
				}
			}

			//MSA_reference_sequences_whichHap[gene][gene+"*03"] = "AAAGAA--ACAAAAGAAAAA";

			auto addRead = [&](std::string readID, unsigned int start, unsigned int stop, std::vector<std::pair<unsigned int, std::string>> alleles) -> void {
				read_start_stop_positions[gene][readID] = std::make_pair(start, stop);
				read_start_per_position[gene][start].insert(readID);
				read_stop_per_position[gene][stop].insert(readID);
				reads_2_genes[readID] = gene;
				for(auto A : alleles)
				{
					read_genotypes_per_position[gene][readID][A.first] = A.second;
				}
				iteration_2_readIDs[gene][1].insert(readID);
			};

			std::cout << "Check1" << "\n" << std::flush;

			addRead("r1", 0, 19, {std::make_pair(3, "T"), std::make_pair(9, "T"), std::make_pair(14, "T")});
			addRead("r2", 0, 19, {std::make_pair(3, "T"), std::make_pair(9, "T"), std::make_pair(14, "T")});
			addRead("r5", 0, 19, {std::make_pair(3, "T"), std::make_pair(9, "T"), std::make_pair(14, "T")});

			addRead("r3", 0, 19, {std::make_pair(3, "C"), std::make_pair(9, "C"), std::make_pair(14, "C")});
			addRead("r4", 0, 19, {std::make_pair(3, "C"), std::make_pair(9, "C"), std::make_pair(14, "C")});
			addRead("r6", 0, 19, {std::make_pair(3, "C"), std::make_pair(9, "C"), std::make_pair(14, "C")});
		}
		else
		{
			{
				std::string inputFn_genes = arguments.at("inputPrefix") + ".genes";
				std::cout << "Now reading: " << inputFn_genes << "\n" << std::flush;

				assert(Utilities::fileExists(inputFn_genes));

				std::vector<std::string> lines_genes = Utilities::getAllLines(inputFn_genes);
				for(auto line : lines_genes)
				{
					Utilities::eraseNL(line);
					if(!line.length())
						continue;
					std::vector<std::string> line_fields = Utilities::split(line, "\t");
					assert(line_fields.size() == 2);
					std::string gene = line_fields.at(0);
					unsigned int this_gene_length = Utilities::StrtoI(line_fields.at(1));
					assert(gene_length.count(gene) == 0);
					gene_length[gene] = this_gene_length;
				}
			}

			{
				std::string inputFn_MSA = arguments.at("inputPrefix") + ".MSA";
				assert(Utilities::fileExists(inputFn_MSA));

				std::vector<std::string> lines_MSA = Utilities::getAllLines(inputFn_MSA);
				for(auto line : lines_MSA)
				{
					Utilities::eraseNL(line);
					if(!line.length())
						continue;
					std::vector<std::string> line_fields = Utilities::split(line, "\t");
					assert(line_fields.size() == 4);
					std::string gene = line_fields.at(0);
					std::string alleleID = line_fields.at(1);
					std::string whichHap = line_fields.at(2);
					std::string alleleSeq = line_fields.at(3);

					assert(MSA_reference_sequences[gene].count(alleleID) == 0);
					assert(alleleSeq.length() == gene_length.at(gene));

					MSA_reference_sequences[gene][alleleID] = alleleSeq;
					MSA_reference_sequences_whichHap[gene][alleleID] = whichHap;
				}
				for(auto gene : MSA_reference_sequences_whichHap)
				{
					std::set<std::string> existingHaps;
					for(auto allele : MSA_reference_sequences_whichHap.at(gene.first))
					{
						existingHaps.insert(allele.second);
					}					
					assert(existingHaps.size()>=1);
					if(existingHaps.size() == 1)
					{
						for(auto allele : MSA_reference_sequences_whichHap.at(gene.first))
						{
							MSA_reference_sequences_whichHap.at(gene.first).at(allele.first) = "?";
						}
					}
				}
			}

			{
				std::string inputFn_read_start_stop = arguments.at("inputPrefix") + ".readCoordinates";
				assert(Utilities::fileExists(inputFn_read_start_stop));

				std::vector<std::string> lines_read_start_stop = Utilities::getAllLines(inputFn_read_start_stop);
				for(auto line : lines_read_start_stop)
				{
					Utilities::eraseNL(line);
					if(!line.length())
						continue;
					std::vector<std::string> line_fields = Utilities::split(line, "\t");
					assert(line_fields.size() == 4);
					assert(read_start_stop_positions.count(line_fields.at(1)) == 0);
					std::string gene = line_fields.at(0);
					std::string readID = line_fields.at(1);
					unsigned int start = Utilities::StrtoI(line_fields.at(2));
					unsigned int stop = Utilities::StrtoI(line_fields.at(3));
					read_start_stop_positions[gene][readID] = std::make_pair(start, stop);
					read_start_per_position[gene][start].insert(readID);
					read_stop_per_position[gene][stop].insert(readID);
					assert(reads_2_genes.count(readID) == 0);
					reads_2_genes[readID] = gene;
				}
			}

			{
				std::string inputFn_read_alleles = arguments.at("inputPrefix") + ".readAlleles";
				assert(Utilities::fileExists(inputFn_read_alleles));

				std::vector<std::string> lines_read_alleles = Utilities::getAllLines(inputFn_read_alleles);
				for(auto line : lines_read_alleles)
				{
					Utilities::eraseNL(line);
					if(!line.length())
						continue;
					std::vector<std::string> line_fields = Utilities::split(line, "\t");
					assert(line_fields.size() == 3);
					unsigned int outerIteration = Utilities::StrtoI(line_fields.at(0));
					std::string readID = line_fields.at(1);
					std::string gene = reads_2_genes.at(readID);

					if(readID_2_iteration.count(readID))
					{
						assert(readID_2_iteration.at(readID) == outerIteration);
					}
					readID_2_iteration[readID] = outerIteration;
					iteration_2_readIDs[gene][outerIteration].insert(readID);
					gene2Iterations[gene].insert(outerIteration);

					assert(read_genotypes_per_position[gene].count(readID) == 0);
					std::vector<std::string> gt_fields = Utilities::split(line_fields.at(2), " ");
					for(unsigned int i = 0; i < gt_fields.size(); i++)
					{
						std::string oneGt = gt_fields.at(i);
						std::vector<std::string> oneGt_fields = Utilities::split(oneGt, ":");
						assert(oneGt_fields.size() == 2);
						unsigned int gt_pos = Utilities::StrtoI(oneGt_fields.at(0));
						std::string gt_value = oneGt_fields.at(1);
						assert(read_genotypes_per_position[gene][readID].count(gt_pos) == 0);
						read_genotypes_per_position[gene][readID][gt_pos] = gt_value;
					}
				}
				
				// assert(read_genotypes_per_position.count("A"));
				// assert(read_genotypes_per_position.at("A").count("HLAA_h0_A*02:90_43_578_0:0:0_1:0:0_68"));
				
				//std::cout << "Debug GT: " << read_genotypes_per_position.at("A").at("HLAA_h0_A*02:90_43_578_0:0:0_1:0:0_68").at(1375) << "\n";
				//assert( 2 == 4 );
			}

			/*
			{
				for(auto gene : read_start_stop_positions)
				{
					for(auto readStartStop : read_start_stop_positions.at(gene.first))
					{
						std::string readID = readStartStop.first;
						std::pair<unsigned int, unsigned int> start_stop = readStartStop.second;
						assert(start_stop.first < start_stop.second);
						for(unsigned int levelI = 0; levelI < )
					}
				}
			}
			*/

			{
				std::string inputFn_activeAlleles = arguments.at("inputPrefix") + ".activeAlleles";
				assert(Utilities::fileExists(inputFn_activeAlleles));

				std::vector<std::string> lines_activeAlleles = Utilities::getAllLines(inputFn_activeAlleles);
				for(auto line : lines_activeAlleles)
				{
					Utilities::eraseNL(line);
					if(!line.length())
						continue;
					std::vector<std::string> line_fields = Utilities::split(line, "\t");
					assert(line_fields.size() >= 2);

					std::string geneID = line_fields.at(0);
					unsigned int position = Utilities::StrtoI(line_fields.at(1));
					assert(position >= 0);
					assert(position < gene_length.at(geneID));

					for(unsigned int i = 2; i < line_fields.size(); i++)
					{
						std::string oneAllele = line_fields.at(i);
						std::vector<std::string> oneAllele_fields = Utilities::split(oneAllele, ";");
						assert(oneAllele_fields.size() == 3);
						activeAlleles_per_position[geneID][position].insert(oneAllele_fields.at(0));
					}
				}
			}
		}

		std::cout << "HMM-based full-gene inference: Have data for " << gene_length.size() << " genes\n" << std::flush;

		std::cout << "Check2" << "\n" << std::flush;

		fullLengthHMM myHMM(
				gene_length,
				reads_2_genes,
				read_start_stop_positions,
				read_start_per_position,
				read_stop_per_position,
				read_genotypes_per_position,
				activeAlleles_per_position,
				MSA_reference_sequences,
				MSA_reference_sequences_whichHap
		);

		std::cout << "Check3" << "\n" << std::flush;


		std::string outputFn = arguments.at("inputPrefix") + ".fullLengthInference.fasta";
		std::string outputFn_graphLevels = arguments.at("inputPrefix") + ".fullLengthInference.fasta.graphLevels";
		
		std::ofstream outputFastaStream;
		outputFastaStream.open(outputFn.c_str(), std::ios::out);
		assert(outputFastaStream.is_open());
		if(! outputFastaStream.is_open())
		{
			throw std::runtime_error("Cannot open file for writing");
		}
		
		std::ofstream outputGraphLevelsStream;
		outputGraphLevelsStream.open(outputFn_graphLevels.c_str(), std::ios::out);
		assert(outputGraphLevelsStream.is_open());

		if(! outputGraphLevelsStream.is_open())
		{
			throw std::runtime_error("Cannot open file for writing");
		}

		std::cout << "Check4" << "\n" << std::flush;


		for(auto gene : gene_length)
		{
			assert(iteration_2_readIDs.at(gene.first).count(1));
			std::cout << "Now making inference for " << gene.first << " -- " << gene2Iterations.at(gene.first).size() << " iterations\n" << std::flush;
			for(unsigned int iterationI = 1; iterationI <= gene2Iterations.at(gene.first).size(); iterationI++)
			{
				assert(gene2Iterations.at(gene.first).count(iterationI));
				std::cout << "\tIteration " << iterationI << "\n";
				myHMM.makeInference(gene.first, outputFastaStream, outputGraphLevelsStream, arguments.at("inputPrefix") + ".fullLengthInference.byGene.i" + gene.first + std::to_string(iterationI) + ".", iteration_2_readIDs.at(gene.first).at(iterationI));
			}
		}
		
		/*
		std::cout << "Now making inference for " << "A" << "\n" << std::flush;
		 myHMM.makeInference("A", outputFastaStream, arguments.at("inputPrefix") + ".fullLengthInference.byGene");

		 std::cout << "Now making inference for " << "DRB1" << "\n" << std::flush;
		 myHMM.makeInference("DRB1", outputFastaStream, arguments.at("inputPrefix") + ".fullLengthInference.byGene");	
		*/
		// std::cout << "Now making inference for " << "DRA" << "\n" << std::flush;
		// myHMM.makeInference("DRA", outputFastaStream, arguments.at("inputPrefix") + ".fullLengthInference.byGene");
		
	}
	else if(arguments.at("action") == "PRGmapping")
	{
		assert(1 == 0); // do we need this function? NB no re-mapping with -a!

		arguments["graph"] = "C:\\Users\\AlexanderDilthey\\OneDrive\\PRG-BWA\\testData\\testPRG";
		arguments["BAMs"] = "C:\\Users\\AlexanderDilthey\\OneDrive\\PRG-BWA\\testData\\NA12878_PRG\\sorted.bam";
		std::string outputDirectory = "";

		assert(arguments.count("graph"));
		assert(arguments.count("BAMs"));

		std::vector<std::string> BAMs = Utilities::split(arguments.at("BAMs"), ",");
		for(size_t BAMi = 0; BAMi < BAMs.size(); BAMi++)
		{
			mapper::processBAM BAMprocessor(arguments.at("graph"));
			std::pair<double, double> IS_estimate = BAMprocessor.estimateInsertSize(BAMs.at(BAMi), true);
			BAMprocessor.alignReads(BAMs.at(BAMi), 0, IS_estimate.first, IS_estimate.second, outputDirectory, true);
		}
	}
	else if(arguments.at("action") == "checkSeq")
	{
		assert(arguments.count("PRG_graph_dir"));
		std::string PRG_graph_dir = arguments.at("PRG_graph_dir");

		std::vector<std::pair<std::string, std::vector<int>>> check_sequences;
		std::string inputfile_sequences = "backbone_sequences.txt";
		std::vector<std::string> lines;

		std::ifstream inputSequences_stream;
		inputSequences_stream.open(inputfile_sequences.c_str());
		assert(inputSequences_stream.is_open());
		std::string line;
		while(inputSequences_stream.good())
		{
			std::getline(inputSequences_stream, line);
			Utilities::eraseNL(line);
			lines.push_back(line);
		}
		assert((lines.size() % 4) == 0);
		int n_seqs = lines.size() / 4;
		for(int i = 0; i < n_seqs; i++)
		{
			int start_index = i * 4; 
			std::string name = lines.at(start_index);
			std::string alignmentSequence = lines.at(start_index+1);
			std::vector<int> alignmentLevels = Utilities::StrtoI(Utilities::split(lines.at(start_index+3), " "));
			if(alignmentSequence.length() != alignmentLevels.size())
			{
				std::cerr << "alignmentSequence.length() != alignmentLevels.size()" << "\n";
				std::cerr << "alignmentSequence.length()" << ": " << alignmentSequence.length() << "\n";
				std::cerr << "alignmentLevels.size()" << ": " << alignmentLevels.size() << "\n";
				std::cerr << "i" << ": " << i << "\n";
				std::cerr << "read name" << ": " << name << "\n";
				std::cerr << std::flush;
				
			}
			assert(alignmentSequence.length() == alignmentLevels.size());
			check_sequences.push_back(make_pair(alignmentSequence, alignmentLevels));
		}
				
		mapper::processBAM BAMprocessor (PRG_graph_dir);
		Graph* g = BAMprocessor.getGraph();
		//g->checkAlignmentBackbonePresences(check_sequences);
		g->checkAlignmentBackbonePresences_ignoreGraphGaps(check_sequences);
	}
	else if(arguments.at("action") == "HLA2")
	{

		// ../bin/HLA-PRG-LA --action HLA --sampleID NA12878 --BAM /gpfs1/well/gsk_hla/temp_mapping_2/NA12878_PRG/merged.bam --outputDirectory /gpfs1/well/gsk_hla/HLA-PRG-LA/working/NA12878_PLATINUM --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended
		// ../bin/HLA-PRG-LA --action HLA --sampleID NA12892_REDUCED --BAM /gpfs1/well/gsk_hla/PRG_Remapping/BAMs/PLATINUM_NA12892/merged.bam --outputDirectory /gpfs1/well/gsk_hla/HLA-PRG-LA/working/NA12892_PLATINUM_RED --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended
	
		assert(arguments.count("sampleID"));
		assert(arguments.count("BAM"));
		assert(arguments.count("outputDirectory"));
		assert(arguments.count("PRG_graph_dir"));

		std::string sampleID = arguments.at("sampleID");
		std::string BAM = arguments.at("BAM");
		std::string outputDirectory = arguments.at("outputDirectory");
		std::string PRG_graph_dir = arguments.at("PRG_graph_dir");

		std::string file_true_HLA_types = (arguments.count("trueHLA") ? arguments.at("trueHLA") : "");

		if(! Utilities::directoryExists(outputDirectory))
		{
			Utilities::makeDir(outputDirectory);
		}

		std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>> inferredHLA;

		assert(Utilities::fileExists(BAM));

		mapper::processBAM BAMprocessor (PRG_graph_dir);
		Graph* g = BAMprocessor.getGraph();

		/*
		std::cout << Utilities::timestamp()  << "Now indexing graph!\n" << std::flush;
		GraphAndEdgeIndex gI(g, 25);
		std::cout << Utilities::timestamp()  << "Indexing done!\n" << std::flush;
		assert(1 == 0);
		*/
		
		std::pair<double, double> IS_estimate = BAMprocessor.estimateInsertSize(BAM, true);

		// debug
		//Utilities::make_or_clearDirectory(outputDirectory);

		std::string BAM1 = outputDirectory + "/BAM_1.bam";
		std::string BAM2 = outputDirectory + "/BAM_2.bam";
		// debug
		if(1 == 0)
		{
			std::cout << Utilities::timestamp() << "Carry out remapping step - input " << BAM << ", output " << BAM1 << ", " << BAM2 << "\n" << std::flush;

			std::string PRGonlyReferenceGenomePath = PRG_graph_dir + "/mapping_PRGonly/referenceGenome.fa";

		    // std::string extendedReferenceGenomePath = Utilities::getFirstLine(PRG_graph_dir + "/extendedReferenceGenomePath.txt");
		    std::string sequencesFile = PRG_graph_dir + "/sequences.txt";

		    std::string dir_extractedReads_for_remapping = outputDirectory + "/extractedReads_forRemapping";
			Utilities::make_or_clearDirectory(dir_extractedReads_for_remapping);

			std::map<std::string, std::pair<int, int>> regions_for_extraction;
			{
				// read regions for extraction

				std::ifstream PRG_covered_regions_stream;
				PRG_covered_regions_stream.open(sequencesFile.c_str());
				assert(PRG_covered_regions_stream.is_open());
				assert(PRG_covered_regions_stream.good());
				std::string line;
				std::getline(PRG_covered_regions_stream, line);
				std::string headerLine = line;
				Utilities::eraseNL(headerLine);				
				std::vector<std::string> headerFields = Utilities::split(headerLine, "\t");
				assert(headerFields.at(0) == "SequenceID");
				assert(headerFields.at(1) == "Name");
				assert(headerFields.at(2) == "FASTAID");
				assert(headerFields.at(3) == "Chr");
				assert(headerFields.at(4) == "Start_1based");
				assert(headerFields.at(5) == "Stop_1based");

				while(PRG_covered_regions_stream.good())
				{
					std::getline(PRG_covered_regions_stream, line);
					Utilities::eraseNL(line);

					if(line.length() == 0)
					{
						continue;
					}

					std::vector<std::string> fields = Utilities::split(line, "\t");
					assert(fields.size() == headerFields.size());

					std::string BAMid;
					int startIndex_0based;
					int stopIndex_0based;

					std::string FastaID = fields.at(2);
					std::string Chr = fields.at(3);
					std::string start_str = fields.at(4);
					std::string stop_str = fields.at(5);
					if(Chr != "")
					{
						BAMid = Chr;
						assert(start_str.length());
						assert(stop_str.length());

						startIndex_0based = Utilities::StrtoI(start_str);
						stopIndex_0based = Utilities::StrtoI(stop_str);
						assert(startIndex_0based >= 0);
						assert(stopIndex_0based >= 0);
						assert(startIndex_0based <= stopIndex_0based);
					}
					else
					{
						assert(FastaID.length());
						BAMid = FastaID;
						assert(start_str.length() == 0);
						assert(stop_str.length() == 0);
						startIndex_0based = -1;
						stopIndex_0based = -1;
					}


					assert(BAMid.length());
					assert(((startIndex_0based == -1) && (stopIndex_0based == -1)) || ((startIndex_0based != -1) && (stopIndex_0based != -1) && (startIndex_0based <= stopIndex_0based)));

					regions_for_extraction[BAMid] = std::make_pair(startIndex_0based, stopIndex_0based);
				}
			}

			linearALTs::linearALTs::extractReadsFromBAM(dir_extractedReads_for_remapping, BAM, regions_for_extraction, true);
			std::string fastq_1 = dir_extractedReads_for_remapping + "/R_1.fq";
			std::string fastq_2 = dir_extractedReads_for_remapping + "/R_2.fq";
			assert(Utilities::fileExists(fastq_1));
			assert(Utilities::fileExists(fastq_2));

			/*
			mapper::bwa::BWAmapper bwaMapper;
			bwaMapper.mapUnpaired(extendedReferenceGenomePath, fastq_1, BAM1, true);
			bwaMapper.mapUnpaired(extendedReferenceGenomePath, fastq_2, BAM2, true);
			*/


			mapper::bowtie2::Bowtie2mapper bowtie2Mapper(pF);
			std::string bowtie_idxname = PRGonlyReferenceGenomePath + "_bowtie2idx";
			bowtie2Mapper.make_sure_ref_is_indexed(PRGonlyReferenceGenomePath, bowtie_idxname);


			bowtie2Mapper.mapUnpaired(bowtie_idxname,
					fastq_1,
					BAM1
			);

			bowtie2Mapper.mapUnpaired(bowtie_idxname,
					fastq_2,
					BAM2
			);

			std::cout << Utilities::timestamp() << "Remapping done.\n" << std::flush;
			assert(Utilities::fileExists(BAM1));
			assert(Utilities::fileExists(BAM1+".bai"));
			assert(Utilities::fileExists(BAM2));
			assert(Utilities::fileExists(BAM2+".bai"));
		} 

		hla::HLATyper HLAtyper(g, PRG_graph_dir, "");

		// BAMprocessor.alignReads(BAM, 0, IS_estimate.first, IS_estimate.second, outputDirectory, &HLAtyper);
		BAMprocessor.alignReads2(BAM1, BAM2, 0, IS_estimate.first, IS_estimate.second, outputDirectory, false, &HLAtyper);

		std::string expected_HLA_type_inference_output = outputDirectory + "/hla/R1_bestguess.txt";
		assert(Utilities::fileExists(expected_HLA_type_inference_output));
		hla::HLATyper::read_inferred_types(sampleID, inferredHLA, expected_HLA_type_inference_output);

		if(file_true_HLA_types.length())
		{
			std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>> trueHLA;
			hla::HLATyper::read_true_types(trueHLA, file_true_HLA_types);
			hla::HLATyper::evaluate_HLA_types(trueHLA, inferredHLA);
		}
	}

	// what follows is the new multi-step mapping
	/*
	else if(arguments.at("action") == "HLA")
	{

		// ../bin/HLA-PRG-LA --action HLA --sampleID NA12878 --BAM /gpfs1/well/gsk_hla/temp_mapping_2/NA12878_PRG/merged.bam --outputDirectory /gpfs1/well/gsk_hla/HLA-PRG-LA/working/NA12878_PLATINUM --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended
		// ../bin/HLA-PRG-LA --action HLA --sampleID NA12892_REDUCED --BAM /gpfs1/well/gsk_hla/PRG_Remapping/BAMs/PLATINUM_NA12892/merged.bam --outputDirectory /gpfs1/well/gsk_hla/HLA-PRG-LA/working/NA12892_PLATINUM_RED --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended

		assert(arguments.count("sampleID"));
		assert(arguments.count("BAM") || (arguments.count("FASTQ1") && arguments.count("FASTQ2")));
		assert(arguments.count("outputDirectory"));
		assert(arguments.count("PRG_graph_dir"));

		std::string sampleID = arguments.at("sampleID");
		std::string BAM = arguments.at("BAM");
		std::string outputDirectory = arguments.at("outputDirectory");
		std::string PRG_graph_dir = arguments.at("PRG_graph_dir");
		if(arguments.count("FASTQ1"))
		{
			assert(arguments.count("mapAgainstCompleteGenome"));
			assert(! arguments.count("BAM"));
		}

		std::string file_true_HLA_types = (arguments.count("trueHLA") ? arguments.at("trueHLA") : "");

		if(! Utilities::directoryExists(outputDirectory))
		{
			Utilities::makeDir(outputDirectory);
		}
		//

		std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>> inferredHLA;

		mapper::processBAM BAMprocessor (PRG_graph_dir);
		std::pair<double, double> IS_estimate;

		if(arguments.count("BAM"))
		{
			assert(Utilities::fileExists(BAM));

			IS_estimate = BAMprocessor.estimateInsertSize(BAM, true);

			bool alwaysStartFromScratch = true;

			if(alwaysStartFromScratch)
			{
				Utilities::make_or_clearDirectory(outputDirectory);
			}
			else
			{
				if(!Utilities::directoryExists(outputDirectory))
				{
					Utilities::makeDir(outputDirectory);
				}
			}

			std::vector<std::string> sampledReferenceGenomes = Utilities::getAllLines(PRG_graph_dir + "/sampledReferenceGenomes.txt");

			// extract reads
			std::string fastq_extracted_forRemapping_1;
			std::string fastq_extracted_forRemapping_2;

			std::string dir_extractedReads_for_remapping = outputDirectory + "/extractedReads_forRemapping";
			if(alwaysStartFromScratch || (!Utilities::directoryExists(dir_extractedReads_for_remapping)))
			{
				Utilities::make_or_clearDirectory(dir_extractedReads_for_remapping);

				std::string PRG_sequences_file = PRG_graph_dir + "/sequences.txt";
				std::map<std::string, std::pair<int, int>> regions_for_extraction;
				{
					// read regions for extraction

					std::ifstream PRG_covered_regions_stream;
					PRG_covered_regions_stream.open(PRG_sequences_file.c_str());
					assert(PRG_covered_regions_stream.is_open());
					assert(PRG_covered_regions_stream.good());
					std::string line;
					std::getline(PRG_covered_regions_stream, line);
					std::string headerLine = line;
					Utilities::eraseNL(headerLine);
					std::vector<std::string> headerFields = Utilities::split(headerLine, "\t");

					assert(headerFields.at(0) == "SequenceID");
					assert(headerFields.at(1) == "Name");
					assert(headerFields.at(2) == "FASTAID");
					assert(headerFields.at(3) == "Chr");
					assert(headerFields.at(4) == "Start_1based");
					assert(headerFields.at(5) == "Stop_1based");

					while(PRG_covered_regions_stream.good())
					{
						std::getline(PRG_covered_regions_stream, line);
						Utilities::eraseNL(line);

						if(line.length() == 0)
						{
							continue;
						}

						std::vector<std::string> fields = Utilities::split(line, "\t");
						assert(fields.size() == headerFields.size());

						std::string BAMid;
						int startIndex_0based;
						int stopIndex_0based;

						std::string FastaID = fields.at(2);
						std::string Chr = fields.at(3);
						std::string start_str = fields.at(4);
						std::string stop_str = fields.at(5);
						if(Chr != "")
						{
							BAMid = Chr;
							assert(start_str.length());
							assert(stop_str.length());

							startIndex_0based = Utilities::StrtoI(start_str);
							stopIndex_0based = Utilities::StrtoI(stop_str);
							assert(startIndex_0based >= 0);
							assert(stopIndex_0based >= 0);
							assert(startIndex_0based <= stopIndex_0based);
						}
						else
						{
							assert(FastaID.length());
							BAMid = FastaID;
							assert(start_str.length() == 0);
							assert(stop_str.length() == 0);
							startIndex_0based = -1;
							stopIndex_0based = -1;
						}


						assert(BAMid.length());
						assert(((startIndex_0based == -1) && (stopIndex_0based == -1)) || ((startIndex_0based != -1) && (stopIndex_0based != -1) && (startIndex_0based <= stopIndex_0based)));

						regions_for_extraction[BAMid] = std::make_pair(startIndex_0based, stopIndex_0based);
					}

				}

				linearALTs::linearALTs::extractReadsFromBAM(dir_extractedReads_for_remapping, BAM, regions_for_extraction);
				fastq_extracted_forRemapping_1 = dir_extractedReads_for_remapping + "/R_1.fq";
				fastq_extracted_forRemapping_2 = dir_extractedReads_for_remapping + "/R_2.fq";
				assert(Utilities::fileExists(fastq_extracted_forRemapping_1));
				assert(Utilities::fileExists(fastq_extracted_forRemapping_2));
			}
			else
			{
				fastq_extracted_forRemapping_1 = dir_extractedReads_for_remapping + "/R_1.fq";
				fastq_extracted_forRemapping_2 = dir_extractedReads_for_remapping + "/R_2.fq";
				assert(Utilities::fileExists(fastq_extracted_forRemapping_1));
				assert(Utilities::fileExists(fastq_extracted_forRemapping_2));
			}


			std::vector<std::string> remapped_BAMs;
			for(unsigned int genomeI = 0; genomeI < sampledReferenceGenomes.size(); genomeI++)
			{
				std::string sampledReferenceGenome = sampledReferenceGenomes.at(genomeI);

				std::string BAM_remapped = outputDirectory + "/sampled_" + Utilities::ItoStr(genomeI) + ".bam";
				std::cout << Utilities::timestamp() << "Carry out remapping step for sampled genome " << genomeI << " - input " << BAM << ", output " << BAM_remapped << "\n" << std::flush;

				if(alwaysStartFromScratch || (! Utilities::fileExists(BAM_remapped)))
				{
					mapper::bwa::BWAmapper bwaMapper(pF);
					bwaMapper.map(sampledReferenceGenome, fastq_extracted_forRemapping_1, fastq_extracted_forRemapping_2, BAM_remapped, true);
				}
				else
				{
					assert(Utilities::fileExists(BAM_remapped));
				}

				std::cout << Utilities::timestamp() << "Remapping done.\n" << std::flush;
				assert(Utilities::fileExists(BAM_remapped+".bai"));
				remapped_BAMs.push_back(BAM_remapped);
			}
		}

		Graph* g = BAMprocessor.getGraph();
		hla::HLATyper HLAtyper(g, PRG_graph_dir, "");

		BAMprocessor.alignReadsMulti(remapped_BAMs, 0, IS_estimate.first, IS_estimate.second, outputDirectory, true, &HLAtyper);

		std::string expected_HLA_type_inference_output = outputDirectory + "/hla/R1_bestguess.txt";
		assert(Utilities::fileExists(expected_HLA_type_inference_output));
		hla::HLATyper::read_inferred_types(sampleID, inferredHLA, expected_HLA_type_inference_output);

		if(file_true_HLA_types.length())
		{
			std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>> trueHLA;
			hla::HLATyper::read_true_types(trueHLA, file_true_HLA_types);
			hla::HLATyper::evaluate_HLA_types(trueHLA, inferredHLA);
		}
	}
	*/
	// what follows is the old one-step remapping procedure
	else if(arguments.at("action") == "remapAndStop")
	{
		bool remap_with_a = true;
		assert(arguments.count("PRG_graph_dir"));
		assert(arguments.count("outputBAM"));
		assert(arguments.count("FASTQ1"));
		assert(arguments.count("FASTQ2"));
		assert(arguments.count("mapAgainstCompleteGenome"));
		assert(arguments.count("longReads"));

		std::string longReads = arguments.at("longReads");
		assert((longReads == "0") || (longReads == "ont2d") || (longReads == "pacbio"));

		if(longReads == "0")
		{
			longReads = "";
		}

		if(longReads.length())
		{
			assert(arguments.count("FASTQU"));
		}

		std::string FASTQU = (arguments.count("FASTQU") ? arguments.at("FASTQU") : "");

		bool mapAgainstCompleteGenome = Utilities::StrtoB(arguments.at("mapAgainstCompleteGenome"));

		unsigned int maxThreads = 1;
		if(arguments.count("maxThreads"))
		{
			maxThreads = Utilities::StrtoI(arguments.at("maxThreads"));
			std::cout << "Set maxThreads to " << maxThreads << "\n" << std::flush;
		}

		mapper::bwa::BWAmapper bwaMapper(pF, maxThreads);

		bwaMapper.createRemappedBAM_forGraph(
				arguments.at("PRG_graph_dir"),
				arguments.at("FASTQ1"),
				arguments.at("FASTQ2"),
				FASTQU,
				mapAgainstCompleteGenome,
				longReads,
				arguments.at("outputBAM"),
				remap_with_a
		);
	}
	else if(arguments.at("action") == "HLA")
	{

		// ../bin/HLA-PRG-LA --action HLA --sampleID NA12878 --BAM /gpfs1/well/gsk_hla/temp_mapping_2/NA12878_PRG/merged.bam --outputDirectory /gpfs1/well/gsk_hla/HLA-PRG-LA/working/NA12878_PLATINUM --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended
		// ../bin/HLA-PRG-LA --action HLA --sampleID NA12892_REDUCED --BAM /gpfs1/well/gsk_hla/PRG_Remapping/BAMs/PLATINUM_NA12892/merged.bam --outputDirectory /gpfs1/well/gsk_hla/HLA-PRG-LA/working/NA12892_PLATINUM_RED --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended

		unsigned int maxThreads = 1;

		assert(arguments.count("sampleID"));
		assert(arguments.count("BAM") || (arguments.count("FASTQ1") && arguments.count("FASTQ2")));
		assert(arguments.count("outputDirectory"));
		assert(arguments.count("PRG_graph_dir"));

		std::string sampleID = arguments.at("sampleID");
		std::string outputDirectory = arguments.at("outputDirectory");
		std::string PRG_graph_dir = arguments.at("PRG_graph_dir");
		if(arguments.count("FASTQ1"))
		{
			assert(arguments.count("mapAgainstCompleteGenome"));
			assert(! arguments.count("BAM"));
		}
		if(arguments.count("maxThreads"))
		{
			maxThreads = Utilities::StrtoI(arguments.at("maxThreads"));
			std::cout << "Set maxThreads to " << maxThreads << "\n" << std::flush;
		}

		std::string file_true_HLA_types = (arguments.count("trueHLA") ? arguments.at("trueHLA") : "");

		if(! Utilities::directoryExists(outputDirectory))
		{
			Utilities::makeDir(outputDirectory); 
		}
		// t

		std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>> inferredHLA;

		mapper::bwa::BWAmapper bwaMapper(pF, maxThreads);

		mapper::processBAM BAMprocessor (PRG_graph_dir, maxThreads);
		std::pair<double, double> IS_estimate;
		std::string BAM_remapped = outputDirectory + "/remapped_with_a.bam";
		std::string PRGonlyReferenceGenomePath = bwaMapper.getPRGonlyReferenceGenomePath_forGraph(PRG_graph_dir);
		/*
		std::string extendedReferenceGenomePath;
		if(Utilities::fileExists(PRG_graph_dir + "/extendedReferenceGenomePath.txt"))
		{
			extendedReferenceGenomePath  = Utilities::getFirstLine(PRG_graph_dir + "/extendedReferenceGenomePath.txt");
		}
		else
		{
			extendedReferenceGenomePath  = PRG_graph_dir + "/extendedReferenceGenome/extendedReferenceGenome.fa";
			assert(Utilities::fileExists(extendedReferenceGenomePath));
		}
		*/
		bool remap_with_a = true;
		if(arguments.count("remap_with_a"))
		{
			remap_with_a = Utilities::StrtoB(arguments.at("remap_with_a"));
		}
		
		bool remapped_against_extended_reference_genome;
		
		std::string longReads;
		if(arguments.count("BAM"))
		{
			assert((! arguments.count("longReads")) || (arguments.at("longReads") == "0"));

			std::string BAM = arguments.at("BAM");

			assert(Utilities::fileExists(BAM));

			IS_estimate = BAMprocessor.estimateInsertSize(BAM, true);
			
			if(1)
			{
				std::cout << Utilities::timestamp() << "Carry out remapping step - input " << BAM << ", output " << BAM_remapped << "\n" << std::flush;

				//if(1 == 1)
				if(1)
				{
				std::string PRG_sequences_file = PRG_graph_dir + "/sequences.txt";

				
				std::string dir_extractedReads_for_remapping = outputDirectory + "/extractedReads_forRemapping";
				Utilities::make_or_clearDirectory(dir_extractedReads_for_remapping);

				std::map<std::string, std::pair<int, int>> regions_for_extraction;
				{
					// read regions for extraction

					std::ifstream PRG_covered_regions_stream;
					PRG_covered_regions_stream.open(PRG_sequences_file.c_str());
					assert(PRG_covered_regions_stream.is_open());
					assert(PRG_covered_regions_stream.good());
					std::string line;
					std::getline(PRG_covered_regions_stream, line);
					std::string headerLine = line;
					Utilities::eraseNL(headerLine);
					std::vector<std::string> headerFields = Utilities::split(headerLine, "\t");

					assert(headerFields.at(0) == "SequenceID");
					assert(headerFields.at(1) == "Name");
					assert(headerFields.at(2) == "FASTAID");
					assert(headerFields.at(3) == "Chr");
					assert(headerFields.at(4) == "Start_1based");
					assert(headerFields.at(5) == "Stop_1based");

					while(PRG_covered_regions_stream.good())
					{
						std::getline(PRG_covered_regions_stream, line);
						Utilities::eraseNL(line);

						if(line.length() == 0)
						{
							continue;
						}

						std::vector<std::string> fields = Utilities::split(line, "\t");
						assert(fields.size() == headerFields.size());

						std::string BAMid;
						int startIndex_0based;
						int stopIndex_0based;

						std::string FastaID = fields.at(2);
						std::string Chr = fields.at(3);
						std::string start_str = fields.at(4);
						std::string stop_str = fields.at(5);
						if(Chr != "")
						{
							BAMid = Chr;
							assert(start_str.length());
							assert(stop_str.length());

							startIndex_0based = Utilities::StrtoI(start_str);
							stopIndex_0based = Utilities::StrtoI(stop_str);
							assert(startIndex_0based >= 0);
							assert(stopIndex_0based >= 0);
							assert(startIndex_0based <= stopIndex_0based);
						}
						else
						{
							assert(FastaID.length());
							BAMid = FastaID;
							assert(start_str.length() == 0);
							assert(stop_str.length() == 0);
							startIndex_0based = -1;
							stopIndex_0based = -1;
						}


						assert(BAMid.length());
						assert(((startIndex_0based == -1) && (stopIndex_0based == -1)) || ((startIndex_0based != -1) && (stopIndex_0based != -1) && (startIndex_0based <= stopIndex_0based)));

						regions_for_extraction[BAMid] = std::make_pair(startIndex_0based, stopIndex_0based);
					}

				}

				linearALTs::linearALTs::extractReadsFromBAM(dir_extractedReads_for_remapping, BAM, regions_for_extraction);
				std::string fastq_1 = dir_extractedReads_for_remapping + "/R_1.fq";
				std::string fastq_2 = dir_extractedReads_for_remapping + "/R_2.fq";
				assert(Utilities::fileExists(fastq_1));
				assert(Utilities::fileExists(fastq_2));

				//bwaMapper.map(extendedReferenceGenomePath, fastq_1, fastq_2, BAM_remapped, true);
				// todo this is new: map only against PRG
				bwaMapper.map(PRGonlyReferenceGenomePath, fastq_1, fastq_2, BAM_remapped, remap_with_a);
				remapped_against_extended_reference_genome = false;
				}

				std::cout << Utilities::timestamp() << "Remapping done.\n" << std::flush;
				assert(Utilities::fileExists(BAM_remapped));
				assert(Utilities::fileExists(BAM_remapped+".bai"));
			}
		}
		else
		{
			assert(arguments.count("FASTQ1"));
			assert(arguments.count("FASTQ2"));
			assert(arguments.count("mapAgainstCompleteGenome"));
			assert(arguments.count("longReads"));
						
			longReads = arguments.at("longReads");
			assert((longReads == "0") || (longReads == "ont2d") || (longReads == "pacbio"));

			if(longReads == "0")
			{
				longReads = "";
			}
			
			if(longReads.length())
			{
				assert(arguments.count("FASTQU"));				
			}
			
			std::string FASTQU = (arguments.count("FASTQU") ? arguments.at("FASTQU") : "");

			bool mapAgainstCompleteGenome = Utilities::StrtoB(arguments.at("mapAgainstCompleteGenome"));

			bwaMapper.createRemappedBAM_forGraph(
					PRG_graph_dir,
					arguments.at("FASTQ1"),
					arguments.at("FASTQ2"),
					FASTQU,
					mapAgainstCompleteGenome,
					longReads,
					BAM_remapped,
					remap_with_a
			);

			if(longReads.length() == 0)
			{
				IS_estimate = BAMprocessor.estimateInsertSize(BAM_remapped, mapAgainstCompleteGenome);
			}
			
			remapped_against_extended_reference_genome = mapAgainstCompleteGenome;
		}

		Graph* g = BAMprocessor.getGraph(); 
		hla::HLATyper HLAtyper(g, PRG_graph_dir, "");  
		
		// BAMprocessor.alignReads(BAM_remapped, 0, IS_estimate.first, IS_estimate.second, outputDirectory, false, &HLAtyper);
		BAMprocessor.alignReads_and_inferHLA(BAM_remapped, 0, IS_estimate.first, IS_estimate.second, outputDirectory, remapped_against_extended_reference_genome, &HLAtyper, 1, longReads);

		std::string expected_HLA_type_inference_output = outputDirectory + "/hla/R1_bestguess.txt";
		assert(Utilities::fileExists(expected_HLA_type_inference_output));
		hla::HLATyper::read_inferred_types(sampleID, inferredHLA, expected_HLA_type_inference_output);

		if(file_true_HLA_types.length())
		{
			std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>> trueHLA;
			hla::HLATyper::read_true_types(trueHLA, file_true_HLA_types);
			hla::HLATyper::evaluate_HLA_types(trueHLA, inferredHLA);
		}
	}
	else if(arguments.at("action") == "KIR")
	{
		std::string KIR_dir = "/gpfs1/well/gsk_hla/HLA-PRG-LA/linearReferenceALTs/KIRGRCh38";
		assert(Utilities::directoryExists(KIR_dir));

		std::string BAM = "/gpfs1/well/gsk_hla/temp_mapping_2/NA12878_PRG/merged.bam";
		assert(Utilities::fileExists(BAM));

		std::string sampleID = "NA12878";
		std::string baseDir = (Utilities::directoryExists("working/KIR") ? "working/KIR" : "../working/KIR");

		std::string workingDir = baseDir + "/" +  sampleID;
		Utilities::make_or_clearDirectory(workingDir);

		linearALTs::linearALTs KIRALTs(KIR_dir, pF);

	    // reads covering the KIR region
	    std::string dir_reads_region = workingDir + "/reads_extendedReferenceGenome";
	    std::string dir_reads_equalLengthHaplotypes = workingDir + "/reads_equalLengthHaplotypes";
	    std::string dir_reads_genes = workingDir + "/reads_explicitGenes";

		Utilities::make_or_clearDirectory(dir_reads_region);
		Utilities::make_or_clearDirectory(dir_reads_equalLengthHaplotypes);
		Utilities::make_or_clearDirectory(dir_reads_genes);

		std::cout << Utilities::timestamp() << "Action 'KIR': Init\n" << std::flush;
		
		std::string genePRGDir = KIRALTs.getGeneGraphPRGDirectory();
		mapper::processBAM BAMprocessor(genePRGDir);

		std::cout << Utilities::timestamp() << "Action 'KIR': Estimate Insert Size\n" << std::flush;
		
		std::pair<double, double> IS_estimate = BAMprocessor.estimateInsertSize_noGraph(BAM);
		std::cout << "Estimate IS from BAM as " << IS_estimate.first << " (mean) / " << IS_estimate.second << " (sd).\n" << std::flush;
		
		std::cout << Utilities::timestamp() << "Action 'KIR': Extract reads\n" << std::flush;
		
	    KIRALTs.extractReads_extendedReferenceGenome(dir_reads_region, BAM);
	    KIRALTs.extractReads_equalLengthHaplotypes(dir_reads_equalLengthHaplotypes, BAM);
		KIRALTs.extractReads_geneGraph(dir_reads_genes, BAM);
		
	    // Genes
	    {
			std::cout << Utilities::timestamp() << "Action 'KIR': Map to gene PRG\n" << std::flush;
		
		    std::string dir_inference_genes = workingDir + "/reads_inference_genes";
			Utilities::make_or_clearDirectory(dir_inference_genes);

		    std::string extendedReferenceGenomePath_genePRG = Utilities::getFirstLine(genePRGDir + "/extendedReferenceGenomePath.txt");
			mapper::bwa::BWAmapper bwaMapper(pF);
			bwaMapper.make_sure_ref_is_indexed(extendedReferenceGenomePath_genePRG);


			std::string FASTQ1 = dir_reads_genes + "/R_1.fq";
			std::string FASTQ2 = dir_reads_genes + "/R_2.fq";
			std::string mapped2GenePRG_BAM = dir_inference_genes + "/mapped_2_genePRG.bam";

			bwaMapper.map(
				extendedReferenceGenomePath_genePRG,
				FASTQ1,
				FASTQ2,
				mapped2GenePRG_BAM
			);

			assert(Utilities::fileExists(mapped2GenePRG_BAM));

			std::pair<double, double> IS_estimate_fromGeneGraph = BAMprocessor.estimateInsertSize(mapped2GenePRG_BAM, true);
			std::cout << "Estimated IS from gene graph (not used!): " << IS_estimate_fromGeneGraph.first << " (mean) / " << IS_estimate_fromGeneGraph.second << " (sd).\n";

			BAMprocessor.alignReads(mapped2GenePRG_BAM, 0, IS_estimate.first, IS_estimate.second, dir_inference_genes, true);
	    }


	    // Haplotypes
	    {
			std::cout << Utilities::timestamp() << "Action 'KIR': Haplotype estimation\n" << std::flush;
		
			mapper::processBAM BAMprocessor(genePRGDir);

			assert(IS_estimate.first > 0);
			
			std::cout << "KIR haplotype estimation, insert size: " << IS_estimate.first << " (mean) / " << IS_estimate.second << " (sd).\n" << std::flush;

			 std::pair<std::string, std::string> inferred_haplotypeIDs = KIRALTs.haplotypeLikelihoods(dir_reads_equalLengthHaplotypes, IS_estimate.first, IS_estimate.second, 101, std::make_pair("", ""));

			 std::cout << "KIR haplotype estimation:\n";
			 std::cout << "\t" << inferred_haplotypeIDs.first << "\n";
			 std::cout << "\t" << inferred_haplotypeIDs.second << "\n";
			 std::cout << std::flush;
	    }

		std::cout << Utilities::timestamp() << "Action 'KIR': Done\n" << std::flush;
				
	}

	else if(arguments.at("action") == "KIRhaplotypesSimulation")
	{
		std::string simulationBaseDir = (Utilities::directoryExists("tmp/KIR") ? "tmp/KIR" : "../tmp/KIR");
		std::string linearALTDir = (Utilities::directoryExists("linearReferenceALTs/KIRGRCh38") ? "linearReferenceALTs/KIRGRCh38" : "../linearReferenceALTs/KIRGRCh38");
		std::string qualityMatrixFile = (Utilities::fileExists("src/simulator/predefinedQualityMatrices/I101_NA12878.txt") ? "src/simulator/predefinedQualityMatrices/I101_NA12878.txt" : "simulator/predefinedQualityMatrices/I101_NA12878.txt");

		linearALTs::linearALTs KIRALTs(linearALTDir, pF);
		simulator::simulator S(qualityMatrixFile);

		std::vector<std::pair<std::map<std::string, std::set<std::string>>, std::map<std::string, std::set<std::string>>>> perSimulation_true_read2gene;
		std::vector< std::pair<std::string, std::string>> perSimulation_true_haplotypeIDs;
		std::vector< std::pair<std::string, std::string>> perSimulation_inferred_haplotypeIDs;

		size_t n_simulations = 1;
		for(unsigned int sI = 0; sI < n_simulations; sI++)
		{
			std::cout << "Simulation " << (sI+1) << " reads\n";
			std::cout << "===================================\n";

			std::string simulationOutputDir = simulationBaseDir + "/S"+ Utilities::ItoStr(sI);
			Utilities::make_or_clearDirectory(simulationOutputDir);

			std::map<std::string, std::set<std::string>> true_read2gene_R1;
			std::map<std::string, std::set<std::string>> true_read2gene_R2;

			 std::pair<std::string, std::string> true_haplotypeIDs = KIRALTs.simulateEqualLengthDiploidHaplotypes(simulationOutputDir, S, true_read2gene_R1, true_read2gene_R2, 5, false);

			 perSimulation_true_read2gene.push_back(make_pair(true_read2gene_R1, true_read2gene_R2));
			 perSimulation_true_haplotypeIDs.push_back(true_haplotypeIDs);

			 std::pair<std::string, std::string> inferred_haplotypeIDs = KIRALTs.haplotypeLikelihoods(simulationOutputDir, 100, 10, 101, true_haplotypeIDs);

			 perSimulation_inferred_haplotypeIDs.push_back(inferred_haplotypeIDs);
		}

		for(unsigned int tI = 0; tI <= 1; tI++)
		{
			double T = (double)tI/2.0;

			for(unsigned int useMLHaplotypes = 0; useMLHaplotypes <= 1; useMLHaplotypes++)
			{
				std::map<std::string, std::map<std::string, int>> reads_trueGene_2_inferredGene;
				std::map<std::string, std::map<std::string, int>> reads_inferredGene_2_trueGene;

				for(unsigned int sI = 0; sI < n_simulations; sI++)
				{
					std::cout << "Analyse simulation " << (sI+1) << " with useMLHaplotypes = " << useMLHaplotypes << " and T = " << T << "\n";
					std::cout << "===================================\n";

					std::string simulationOutputDir = simulationBaseDir + "/S"+ Utilities::ItoStr(sI);

					std::map<std::string, std::set<std::string>> true_read2gene_R1 = perSimulation_true_read2gene.at(sI).first;
					std::map<std::string, std::set<std::string>> true_read2gene_R2 = perSimulation_true_read2gene.at(sI).second;

					 std::pair<std::string, std::string> true_haplotypeIDs = perSimulation_true_haplotypeIDs.at(sI);

					std::pair<std::string, std::string> inferred_haplotypeIDs = perSimulation_inferred_haplotypeIDs.at(sI);

					 std::set<std::string> inferred_haplotypeIDs_set;

					 if(useMLHaplotypes)
					 {
						 inferred_haplotypeIDs_set.insert(inferred_haplotypeIDs.first);
						 inferred_haplotypeIDs_set.insert(inferred_haplotypeIDs.second);
					 }

					 KIRALTs.reads2Genes(simulationOutputDir, 100, 10, 101, true_read2gene_R1, true_read2gene_R2, reads_trueGene_2_inferredGene, reads_inferredGene_2_trueGene, inferred_haplotypeIDs_set, T);

				}
				 // KIRALTs.haplotypeLikelihoods(simulationOutputDir, 100, 10, 101, make_pair("GI568335980", "GI568335930"));

				std::set<std::string> genes_set;
				for(auto geneData : reads_trueGene_2_inferredGene)
				{
					genes_set.insert(geneData.first);
				}
				for(auto geneData : reads_inferredGene_2_trueGene)
				{
					genes_set.insert(geneData.first);
				}

				std::string output_reads_trueGene_2_inferredGene = "reads_trueGene_2_inferredGene_useML" + Utilities::ItoStr(useMLHaplotypes) + "_T" +  Utilities::DtoStr(T) + ".txt";
				std::string output_reads_inferredGene_2_trueGene = "reads_inferredGene_2_trueGene_useML" + Utilities::ItoStr(useMLHaplotypes) + "_T" +  Utilities::DtoStr(T) + ".txt";

				std::ofstream outputStream_trueGene_2_inferredGene;
				std::ofstream outputStream_inferredGene_2_trueGene;
				outputStream_trueGene_2_inferredGene.open(output_reads_trueGene_2_inferredGene.c_str());
				outputStream_inferredGene_2_trueGene.open(output_reads_inferredGene_2_trueGene.c_str());
				assert(outputStream_trueGene_2_inferredGene.is_open());
				assert(outputStream_inferredGene_2_trueGene.is_open());

				std::vector<std::string> genes(genes_set.begin(), genes_set.end());

				outputStream_trueGene_2_inferredGene << "SourceGene" << "\t" << "N" << "\t" << Utilities::join(genes, "\t") << "\n";
				outputStream_inferredGene_2_trueGene << "InferredGene" << "\t" << "N" << "\t" << Utilities::join(genes, "\t") << "\n";
				for(auto gene1 : genes)
				{
					std::vector<std::string> outputFields_trueGene_2_inferredGene;
					std::vector<std::string> outputFields_inferredGene_2_trueGene;

					outputFields_trueGene_2_inferredGene.push_back(gene1);
					outputFields_inferredGene_2_trueGene.push_back(gene1);

					int totalReads_trueGene = 0;
					int totalReads_inferredGene = 0;
					for(auto gene2: genes)
					{
						if(reads_trueGene_2_inferredGene[gene1].count(gene2))
						{
							totalReads_trueGene += reads_trueGene_2_inferredGene[gene1].at(gene2);
						}
						if(reads_inferredGene_2_trueGene[gene1].count(gene2))
						{
							totalReads_inferredGene += reads_inferredGene_2_trueGene[gene1].at(gene2);
						}
					}

					outputFields_trueGene_2_inferredGene.push_back(Utilities::ItoStr(totalReads_trueGene));
					outputFields_inferredGene_2_trueGene.push_back(Utilities::ItoStr(totalReads_inferredGene));

					for(auto gene2: genes)
					{
						double prop_trueGene = 0;
						if(reads_trueGene_2_inferredGene[gene1].count(gene2))
						{
							prop_trueGene = (double)reads_trueGene_2_inferredGene[gene1].at(gene2) / (double)totalReads_trueGene;
						}
						double prop_inferredGene = 0;
						if(reads_inferredGene_2_trueGene[gene1].count(gene2))
						{
							prop_inferredGene = (double)reads_inferredGene_2_trueGene[gene1].at(gene2) / (double)totalReads_inferredGene;
						}

						outputFields_trueGene_2_inferredGene.push_back(Utilities::DtoStr(prop_trueGene));
						outputFields_inferredGene_2_trueGene.push_back(Utilities::DtoStr(prop_inferredGene));
					}

					outputStream_trueGene_2_inferredGene << Utilities::join(outputFields_trueGene_2_inferredGene, "\t") << "\n";
					outputStream_inferredGene_2_trueGene << Utilities::join(outputFields_inferredGene_2_trueGene, "\t") << "\n";
				}
			}
		}
	}
	else if(arguments.at("action") == "testRealBAM")
	{
		
		/*
		std::string sampleID = "NA12878";
		std::string BAM = "/gpfs1/well/gsk_hla/temp_mapping_2/NA12878_PRG/merged.bam";
		std::string outputDirectory = "/gpfs1/well/gsk_hla/HLA-PRG-LA/working/NA12878_PLATINUM";
		*/
		
		/*
		std::string sampleID = "NA19240";
		std::string BAM = "/gpfs1/well/gsk_hla/temp_mapping_2/NA19240_PRG/merged.bam";
		std::string outputDirectory = "/gpfs1/well/gsk_hla/HLA-PRG-LA/working/NA19240_1000G";
		*/

		/*
		
		std::string sampleID = "NA19240_REDUCED";
		std::string BAM = "/gpfs1/well/gsk_hla/temp_mapping_2/NA19240_PRG_REDUCED/merged.bam";
		std::string outputDirectory = "/gpfs1/well/gsk_hla/HLA-PRG-LA/working/NA19240_1000G_RED";
		
		*/

		std::string sampleID = "NA12892_REDUCED";
		std::string BAM = "/gpfs1/well/gsk_hla/PRG_Remapping/BAMs/PLATINUM_NA12892/merged.bam";
		std::string outputDirectory = "/gpfs1/well/gsk_hla/HLA-PRG-LA/working/NA12892_PLATINUM_RED";
		
		
		std::string PRG_graph_dir = "/gpfs1/well/gsk_hla/HLA-PRG-LA/graphs/PRG_MHC_GRCh38_withIMGT";		
		std::string file_true_HLA_types = "/Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended";
		
		
		std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>> trueHLA;
		std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>> inferredHLA;
		hla::HLATyper::read_true_types(trueHLA, file_true_HLA_types);
		

		assert(Utilities::fileExists(BAM));

		mapper::processBAM BAMprocessor(PRG_graph_dir);
		std::pair<double, double> IS_estimate = BAMprocessor.estimateInsertSize(BAM, true);

		Utilities::make_or_clearDirectory(outputDirectory);

		Graph* g = BAMprocessor.getGraph();
		hla::HLATyper HLAtyper(g, PRG_graph_dir, "");

		BAMprocessor.alignReads(BAM, 0, IS_estimate.first, IS_estimate.second, outputDirectory, true, &HLAtyper);
		
		std::string expected_HLA_type_inference_output = outputDirectory + "/hla/R1_bestguess.txt";
		assert(Utilities::fileExists(expected_HLA_type_inference_output));		
		hla::HLATyper::read_inferred_types(sampleID, inferredHLA, expected_HLA_type_inference_output);

		hla::HLATyper::evaluate_HLA_types(trueHLA, inferredHLA);
		
	}	
	else if(arguments.at("action") == "testCheckPresence")
	{
		std::string allGraphsOutputDir = (Utilities::directoryExists("tmp/simulatedGraphs") ? "tmp/simulatedGraphs" : "../tmp/simulatedGraphs");
		std::string allGraphsWorkingDir = (Utilities::directoryExists("tmp/working") ? "tmp/working" : "../tmp/working");

		for(unsigned int outerIteration = 1; outerIteration <= 20; outerIteration++)
		{
			simpleGraphSimulator graphSimulator(false, pF);

			std::string graphName = "G" + Utilities::ItoStr(outerIteration);
			std::string thisGraphOutputDir = allGraphsOutputDir + "/" + graphName;
			std::cout << "Directory: " << thisGraphOutputDir << "\n" << std::flush;

			Utilities::make_or_clearDirectory(thisGraphOutputDir);
			assert(Utilities::directoryExists(thisGraphOutputDir));

			std::string readAlignmentOutputDir = allGraphsWorkingDir + "/" + graphName;
			Utilities::make_or_clearDirectory(readAlignmentOutputDir);

			graphSimulator.storeLikeRealPRG(thisGraphOutputDir);

			std::string FASTQ_outputDirectory = thisGraphOutputDir + "/FASTQ";
			Utilities::make_or_clearDirectory(FASTQ_outputDirectory);


			std::string BAM_outputDirectory = thisGraphOutputDir + "/BAM";
			Utilities::make_or_clearDirectory(BAM_outputDirectory);

			Graph g;
			g.readFromFile(thisGraphOutputDir + "/PRG/graph.txt");

			std::cout << "Sequence presence check:\n" << std::flush;
			int nHaplotypes = 10;
			std::vector<std::string> H = g.simulateHaplotypes(nHaplotypes, true);
			for(int haplotypeI = 0; haplotypeI < nHaplotypes; haplotypeI++)
			{
				std::string h = H.at(haplotypeI);
				bool isIn = g.checkSequencePresence(h);
				std::cout << "\t" << haplotypeI << "\t" << h.length() << "\t" << ((isIn) ? "YES" : "NO") << "\n" << std::flush;
				assert(isIn);
			}
		}
	}
	else if(arguments.at("action") == "checkKIRgraph")
	{
		std::string graphDir = "/home/jk/work/prg_input/lrc_kir";
		std::string alignment = "/gpfs1/well/gsk_hla/kir_all_minus_GI568335996.fa";

		std::map<std::string, std::string> inputAlignment = Utilities::readFASTA(alignment);

		Graph g;
		g.readFromFile(graphDir + "/graph.txt");

		std::cerr << "Loaded\n\tGraph " << graphDir << "\n\tAlignment " << alignment << "\n - now checking presence!\n" << std::flush;
		unsigned int L;
		for(std::map<std::string, std::string>::iterator seqIt = inputAlignment.begin(); seqIt != inputAlignment.end(); seqIt++ )
		{
			std::string seqID = seqIt->first;
			std::string S = seqIt->second;

			std::transform(S.begin(), S.end(), S.begin(), ::toupper);
						
			bool isIn = g.checkSequencePresence(S);
			std::cerr << "\t" << seqID << "\t" << S.length() << "\t" << isIn << "\n" << std::flush;

			if(! isIn)
			{
				g.checkSequencePresence(S, true);
			}
			
			if(seqIt == inputAlignment.begin())
			{
				L = S.length();
			}
			else
			{
				assert(L == S.length());
			}
		} 
	}
	else if(arguments.at("action") == "KIRgeneSimulation")
	{
		std::string qualityMatrixFile = (Utilities::fileExists("src/simulator/predefinedQualityMatrices/I101_NA12878.txt") ? "src/simulator/predefinedQualityMatrices/I101_NA12878.txt" : "simulator/predefinedQualityMatrices/I101_NA12878.txt");
		std::string linearALTDir = (Utilities::directoryExists("linearReferenceALTs/KIRGRCh38") ? "linearReferenceALTs/KIRGRCh38" : "../linearReferenceALTs/KIRGRCh38");

		linearALTs::linearALTs KIRALTs(linearALTDir, pF);

		bool simulateFromGraph = false;

		std::string tmpDir = (Utilities::directoryExists("tmp") ? "tmp" : "../tmp/");
		assert(Utilities::directoryExists(tmpDir));

		std::string genePRGDir = KIRALTs.getGeneGraphPRGDirectory();
	    std::string extendedReferenceGenomePath_genePRG = Utilities::getFirstLine(genePRGDir + "/extendedReferenceGenomePath.txt");
		// std::cerr << "extendedReferenceGenomePath_genePRG: " << extendedReferenceGenomePath_genePRG << "\n";
	    assert(Utilities::fileExists(extendedReferenceGenomePath_genePRG));
		mapper::bwa::BWAmapper bwaMapper(pF);
		bwaMapper.make_sure_ref_is_indexed(extendedReferenceGenomePath_genePRG);

		std::string tmpDirKIRSimulations = tmpDir + "/KIRSimulations";
		Utilities::make_or_clearDirectory(tmpDirKIRSimulations);

		mapper::processBAM BAMprocessor(genePRGDir);
		
		double IS_mean = 100;
		double IS_sd = 10;
		simulator::simulator S(qualityMatrixFile, 101, IS_mean, IS_sd);

		std::cout << "TestKIRTyping\n========================================\n";
		std::cout << "simulateFromGraph = " << simulateFromGraph << "\n";

		std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>> trueHLA;
		std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>> inferredHLA;
		for(unsigned int KIRsimulationI = 1; KIRsimulationI <= 5; KIRsimulationI++)
		{
			std::string tmpDirThisSimulation = tmpDirKIRSimulations + "/S" + Utilities::ItoStr(KIRsimulationI);
			Utilities::make_or_clearDirectory(tmpDirThisSimulation);
			
			std::string simulatedReadsDir = tmpDirThisSimulation + "/simulatedReads";
			Utilities::make_or_clearDirectory(simulatedReadsDir);

			KIRALTs.simulateGeneSequences(simulatedReadsDir, S, simulateFromGraph, 10, true);
			std::string FASTQ1 = simulatedReadsDir + "/R_1.fq";
			std::string FASTQ2 = simulatedReadsDir + "/R_2.fq";
			assert(Utilities::fileExists(FASTQ1));
			assert(Utilities::fileExists(FASTQ2));

			simulator::trueReadLevels trueLevels(
				simulatedReadsDir + "/R_1.levels",
				simulatedReadsDir + "/R_2.levels"
			);

			mapper::bwa::BWAmapper bwaMapper(pF);
			bwaMapper.map(
				extendedReferenceGenomePath_genePRG,
				FASTQ1,
				FASTQ2,
				tmpDirThisSimulation + "/output.bam"
			);

			std::string BAM = tmpDirThisSimulation + "/output.bam";
			assert(Utilities::fileExists(BAM));

			std::pair<double, double> IS_estimate = BAMprocessor.estimateInsertSize(BAM, true);

			std::cout << "Estimate IS as " << IS_estimate.first << " (mean) / " << IS_estimate.second << " (sd).\n";

			BAMprocessor.alignReads(BAM, &trueLevels, IS_estimate.first, IS_estimate.second, tmpDirThisSimulation, true);
			// BAMprocessor.alignReads(0, IS_estimate.fBirst, IS_estimate.second);

			assert(trueLevels.get_total_and_correct().first != 0);
			std::cout << "Quality:\n";
			std::cout << "\tGraph: " << trueLevels.get_total_and_correct().first << " " << ((double)trueLevels.get_total_and_correct().second / trueLevels.get_total_and_correct().first) << "\n";
			std::cout << std::flush;
		}
	}	
	else if(arguments.at("action") == "TestHLATyping")
	{
		std::string qualityMatrixFile = (Utilities::fileExists("src/simulator/predefinedQualityMatrices/I101_NA12878.txt") ? "src/simulator/predefinedQualityMatrices/I101_NA12878.txt" : "simulator/predefinedQualityMatrices/I101_NA12878.txt");

		std::string tmpDir = (Utilities::directoryExists("tmp") ? "tmp" : "../tmp/");
		assert(Utilities::directoryExists(tmpDir));

		// make sure reference genome is indexed
		std::string PRG_graph_dir = "/gpfs1/well/gsk_hla/HLA-PRG-LA/graphs/PRG_MHC_GRCh38_withIMGT";
	    std::string extendedReferenceGenomePath = Utilities::getFirstLine(PRG_graph_dir + "/extendedReferenceGenomePath.txt");
	    assert(Utilities::fileExists(extendedReferenceGenomePath));
		mapper::bwa::BWAmapper bwaMapper(pF);
		bwaMapper.make_sure_ref_is_indexed(extendedReferenceGenomePath);

		std::string tmpDirHLASimulations = tmpDir + "/HLATypeSimulations";
		Utilities::make_or_clearDirectory(tmpDirHLASimulations);

		mapper::processBAM BAMprocessor(PRG_graph_dir);
		Graph* g = BAMprocessor.getGraph();
		hla::HLATyper* HLAtyper = new hla::HLATyper(g, PRG_graph_dir, qualityMatrixFile);

		double IS_mean = 10;
		double IS_sd = 10;

		std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>> trueHLA;
		std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>> inferredHLA;
		for(unsigned int HLAsimulationI = 1; HLAsimulationI <= 1; HLAsimulationI++)
		{
			std::string tmpDirThisSimulation = tmpDirHLASimulations + "/S" + Utilities::ItoStr(HLAsimulationI);
			Utilities::make_or_clearDirectory(tmpDirThisSimulation);

			// get list of loci - intron/exon files
			// get list of locus alleles
			// for each locus:
			//	- diploid choice
			// 	- 2 strings
			//  - find graph level of first position of string (e.g. from search in complete list of loci)
			//  - simulate reads
			//  - amend simulated true levels by determined offset

			// decide on HLA type
			// simulate FASTQ
			// generate BAM from FASTQ

			std::string FASTQ1 = tmpDirThisSimulation + "/R_1.fq";   
			std::string FASTQ2 = tmpDirThisSimulation + "/R_2.fq";
			std::string trueTypes = tmpDirThisSimulation + "/HLAtypes.txt";

			std::cout << "Simulating one indivdual HLA reads - " << FASTQ1 << " - " << FASTQ2 << "\n" << std::flush;
			std::cout << "\t" << trueTypes << "\n" << std::flush;

			HLAtyper->simulateOneIndividual(tmpDirThisSimulation, IS_mean, IS_sd, true);

			std::string file_true_HLA_types = tmpDirThisSimulation + "/HLAtypes.txt";
			assert(Utilities::fileExists(file_true_HLA_types));
			hla::HLATyper::read_true_types(trueHLA, file_true_HLA_types);
			
			std::string BAM = tmpDirThisSimulation + "/reads.bam";

			bwaMapper.map(
						extendedReferenceGenomePath,
						FASTQ1,
						FASTQ2,
						BAM,
						true
			);

			std::cout << "Generated BAM: " << BAM << "\n" << std::flush;

			BAMprocessor.alignReads(BAM, 0, IS_mean, IS_sd, tmpDirThisSimulation, HLAtyper);
			
			std::string expected_HLA_type_inference_output = tmpDirThisSimulation + "/hla/R1_bestguess.txt";
			assert(Utilities::fileExists(expected_HLA_type_inference_output));
			
			hla::HLATyper::read_inferred_types(tmpDirThisSimulation, inferredHLA, expected_HLA_type_inference_output);
		}
		
		hla::HLATyper::evaluate_HLA_types(trueHLA, inferredHLA);
	}
	else if(arguments.at("action") == "prepareGraph")
	{
		std::cout << "prepareGraph\n" << std::flush;
		
		std::string PRG_graph_dir = arguments.at("PRG_graph_dir");
		Graph* g = new Graph();

		std::cout << Utilities::timestamp() << "Read graph from " << PRG_graph_dir << "\n" << std::flush;
		g->readFromFile(PRG_graph_dir + "/PRG/graph.txt");
		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

		g->checkStructure();

		// pre-gap-path-indexing serialization
		std::string graph_serialized_fn_preGapPathIndex = PRG_graph_dir + "/serializedGRAPH_preGapPathIndex";
		std::cout << Utilities::timestamp() << "Now serialize (pre-gap-pathindexed) graph to " << graph_serialized_fn_preGapPathIndex << "\n" << std::flush;

		std::ofstream serialization_preGapPath_ostream(graph_serialized_fn_preGapPathIndex);
		if(! serialization_preGapPath_ostream.is_open())
		{
				throw std::runtime_error(" deterministicAnalysis::kickOff: Cannot open file " + graph_serialized_fn_preGapPathIndex + " to read serialized genomeMix\n");
		}
		boost::archive::text_oarchive preGapPathIndex_archive(serialization_preGapPath_ostream);
		preGapPathIndex_archive & g;
		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;
		
		// compute gap-paths
		std::cout << Utilities::timestamp() << "Compute Gap-Edge paths for graph in " << PRG_graph_dir << "\n" << std::flush;
		g->computeGapEdgePaths();
		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;

		// final serialization
		std::string graph_serialized_fn = PRG_graph_dir + "/serializedGRAPH";
		std::cout << Utilities::timestamp() << "Now serialize graph to " << graph_serialized_fn << "\n" << std::flush;

		std::ofstream serialization_ostream(graph_serialized_fn);
		if(! serialization_ostream.is_open())
		{
				throw std::runtime_error(" deterministicAnalysis::kickOff: Cannot open file " + graph_serialized_fn + " to read serialized genomeMix\n");
		}
		boost::archive::text_oarchive archive(serialization_ostream);
		archive & g;

		std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;
	}
	else if(arguments.at("action") == "testPRGMappingUnpaired")
	{
		std::string qualityMatrixFile = (Utilities::fileExists("src/simulator/predefinedQualityMatrices/I101_NA12878.txt") ? "src/simulator/predefinedQualityMatrices/I101_NA12878.txt" : "simulator/predefinedQualityMatrices/I101_NA12878.txt");
		std::string allGraphsOutputDir = (Utilities::directoryExists("tmp/simulatedGraphs") ? "tmp/simulatedGraphs" : "../tmp/simulatedGraphs");
		std::string allGraphsWorkingDir = (Utilities::directoryExists("tmp/working") ? "tmp/working" : "../tmp/working");

		assert(Utilities::directoryExists(allGraphsOutputDir));
		assert(Utilities::directoryExists(allGraphsWorkingDir));

		bool switchContigs = false;
		for(unsigned int outerIteration = 1; outerIteration <= 10; outerIteration++)
		{
			simpleGraphSimulator graphSimulator((bool)switchContigs, pF);

			std::string graphName = "G" + Utilities::ItoStr(outerIteration);
			std::string thisGraphOutputDir = allGraphsOutputDir + "/" + graphName;
			std::cout << "Directory: " << thisGraphOutputDir << "\n" << std::flush;

			Utilities::make_or_clearDirectory(thisGraphOutputDir);
			assert(Utilities::directoryExists(thisGraphOutputDir));

			std::string readAlignmentOutputDir = allGraphsWorkingDir + "/" + graphName;
			Utilities::make_or_clearDirectory(readAlignmentOutputDir);

			graphSimulator.storeLikeRealPRG(thisGraphOutputDir);

			std::string FASTQ_outputDirectory = thisGraphOutputDir + "/FASTQ";
			Utilities::make_or_clearDirectory(FASTQ_outputDirectory);


			std::string BAM_outputDirectory = thisGraphOutputDir + "/BAM";
			Utilities::make_or_clearDirectory(BAM_outputDirectory);


			mapper::bowtie2::Bowtie2mapper bowtie2Mapper(pF);
			std::string bowtie_idxname = thisGraphOutputDir + "/bowtie2IDX";
			bowtie2Mapper.createIndex(thisGraphOutputDir + "/referenceGenome/ref.fa", bowtie_idxname);

			Graph g;
			g.readFromFile(thisGraphOutputDir + "/PRG/graph.txt");
			g.computeGapEdgePaths();
			{
				std::cerr << "Start serialization\n" << std::flush;
				std::string graph_serialized_fn = thisGraphOutputDir + "/graphSerialized";
				std::ofstream serialization_ostream(graph_serialized_fn);
				if(! serialization_ostream.is_open())
				{
						throw std::runtime_error(" deterministicAnalysis::kickOff: Cannot open file " + graph_serialized_fn + " to read serialized genomeMix\n");
				}
				boost::archive::text_oarchive archive(serialization_ostream);
				archive & g;

				/*
				GraphAndEdgeIndex gI(&g, 25);

				std::cerr << "Start serialization index\n" << std::flush;
				std::string index_serialized_fn = thisGraphOutputDir + "/indexSerialized";
				std::ofstream index_serialization_ostream(index_serialized_fn);
				if(! index_serialization_ostream.is_open())
				{
						throw std::runtime_error(" deterministicAnalysis::kickOff: Cannot open file " + index_serialized_fn + " to read serialized genomeMix\n");
				}
				boost::archive::text_oarchive archive2(index_serialization_ostream);
				archive2 & gI;
				*/
			}
			std::cerr << "\n\nDone\n\n" << std::flush;

			std::string qualityMatrixFile = (Utilities::fileExists("src/simulator/predefinedQualityMatrices/I101_NA12878.txt") ? "src/simulator/predefinedQualityMatrices/I101_NA12878.txt" : "simulator/predefinedQualityMatrices/I101_NA12878.txt");
			simulator::simulator S(qualityMatrixFile);
			S.simulateFromGraph(&g, 1, FASTQ_outputDirectory, 5, true);

			// assert( 1 == 0 );

			simulator::trueReadLevels trueLevels(
					thisGraphOutputDir + "/FASTQ/R_1.levels",
					thisGraphOutputDir + "/FASTQ/R_2.levels"
			);

//			assert( 1 == 0 );

			mapper::bwa::BWAmapper bwaMapper(pF);
			bwaMapper.map(
						thisGraphOutputDir + "/referenceGenome/ref.fa",
						thisGraphOutputDir + "/FASTQ/R_1.fq",
						thisGraphOutputDir + "/FASTQ/R_2.fq",
						BAM_outputDirectory + "/output.bam",
						true
			);
			bwaMapper.map(
						thisGraphOutputDir + "/referenceGenome/ref.fa",
						thisGraphOutputDir + "/FASTQ/R_1.fq",
						thisGraphOutputDir + "/FASTQ/R_2.fq",
						BAM_outputDirectory + "/output.bam",
						true
			);

			bowtie2Mapper.mapUnpaired(bowtie_idxname,
					thisGraphOutputDir + "/FASTQ/R_1.fq",
					BAM_outputDirectory + "/output_1.bam"
			);

			bowtie2Mapper.mapUnpaired(bowtie_idxname,
					thisGraphOutputDir + "/FASTQ/R_2.fq",
					BAM_outputDirectory + "/output_2.bam"
			);

			/*
			bwaMapper.mapUnpaired(
						thisGraphOutputDir + "/referenceGenome/ref.fa",
						thisGraphOutputDir + "/FASTQ/R_1.fq",
						BAM_outputDirectory + "/output_1.bam",
						true
			);

			bwaMapper.mapUnpaired(
						thisGraphOutputDir + "/referenceGenome/ref.fa",
						thisGraphOutputDir + "/FASTQ/R_2.fq",
						BAM_outputDirectory + "/output_2.bam",
						true
			);

			*/

			std::string BAM = thisGraphOutputDir + "/BAM/output.bam";
			std::string BAM1 = thisGraphOutputDir + "/BAM/output_1.bam";
			std::string BAM2 = thisGraphOutputDir + "/BAM/output_2.bam";
			assert(Utilities::fileExists(BAM));
			assert(Utilities::fileExists(BAM1));
			assert(Utilities::fileExists(BAM2));

			mapper::processBAM BAMprocessor(thisGraphOutputDir);
			std::pair<double, double> IS_estimate = BAMprocessor.estimateInsertSize(BAM, true);

			//BAMprocessor.alignReads(BAM, &trueLevels, IS_estimate.first, IS_estimate.second, readAlignmentOutputDir);
			BAMprocessor.alignReads2(BAM1, BAM2, &trueLevels, IS_estimate.first, IS_estimate.second, readAlignmentOutputDir, true);

			// BAMprocessor.alignReads(0, IS_estimate.fBirst, IS_estimate.second);

			assert(trueLevels.get_total_and_correct().first != 0);
			std::cout << "Quality:\n";
			std::cout << "\tGraph: " << trueLevels.get_total_and_correct().first << " " << ((double)trueLevels.get_total_and_correct().second / trueLevels.get_total_and_correct().first) << "\n";
			std::cout << std::flush;

//			 assert(1  == 0 );
		}
	}
	else if(arguments.at("action") == "testPRGMapping")
	{
		std::string qualityMatrixFile = (Utilities::fileExists("src/simulator/predefinedQualityMatrices/I101_NA12878.txt") ? "src/simulator/predefinedQualityMatrices/I101_NA12878.txt" : "simulator/predefinedQualityMatrices/I101_NA12878.txt");
		std::string allGraphsOutputDir = (Utilities::directoryExists("tmp/simulatedGraphs") ? "tmp/simulatedGraphs" : "../tmp/simulatedGraphs");
		std::string allGraphsWorkingDir = (Utilities::directoryExists("tmp/working") ? "tmp/working" : "../tmp/working");

		assert(Utilities::directoryExists(allGraphsOutputDir));
		assert(Utilities::directoryExists(allGraphsWorkingDir));

		bool switchContigs = false;
		for(unsigned int outerIteration = 1; outerIteration <= 1; outerIteration++)
		{
			simpleGraphSimulator graphSimulator((bool)switchContigs, pF);

			std::string graphName = "G" + Utilities::ItoStr(outerIteration);
			std::string thisGraphOutputDir = allGraphsOutputDir + "/" + graphName;
			std::cout << "Directory: " << thisGraphOutputDir << "\n" << std::flush;

			Utilities::make_or_clearDirectory(thisGraphOutputDir);
			assert(Utilities::directoryExists(thisGraphOutputDir));

			std::string readAlignmentOutputDir = allGraphsWorkingDir + "/" + graphName;
			Utilities::make_or_clearDirectory(readAlignmentOutputDir);

			graphSimulator.storeLikeRealPRG(thisGraphOutputDir);

			std::string FASTQ_outputDirectory = thisGraphOutputDir + "/FASTQ";
			Utilities::make_or_clearDirectory(FASTQ_outputDirectory);


			std::string BAM_outputDirectory = thisGraphOutputDir + "/BAM";
			Utilities::make_or_clearDirectory(BAM_outputDirectory);

			Graph g;
			g.readFromFile(thisGraphOutputDir + "/PRG/graph.txt");
			g.computeGapEdgePaths();
			{
				std::cerr << "Start serialization\n" << std::flush;
				std::string graph_serialized_fn = thisGraphOutputDir + "/graphSerialized";
				std::ofstream serialization_ostream(graph_serialized_fn);
				if(! serialization_ostream.is_open())
				{
						throw std::runtime_error(" deterministicAnalysis::kickOff: Cannot open file " + graph_serialized_fn + " to read serialized genomeMix\n");
				}
				boost::archive::text_oarchive archive(serialization_ostream);
				archive & g;
			}
			std::cerr << "\n\nDone\n\n" << std::flush;

//			assert(1 == 0);
			std::string qualityMatrixFile = (Utilities::fileExists("src/simulator/predefinedQualityMatrices/I101_NA12878.txt") ? "src/simulator/predefinedQualityMatrices/I101_NA12878.txt" : "simulator/predefinedQualityMatrices/I101_NA12878.txt");
			simulator::simulator S(qualityMatrixFile);
			S.simulateFromGraph(&g, 1, FASTQ_outputDirectory, 2, true);

			// assert( 1 == 0 );

			simulator::trueReadLevels trueLevels(
					thisGraphOutputDir + "/FASTQ/R_1.levels",
					thisGraphOutputDir + "/FASTQ/R_2.levels"
			);

//			assert( 1 == 0 );

			mapper::bwa::BWAmapper bwaMapper(pF);
			bwaMapper.map(
						thisGraphOutputDir + "/referenceGenome/ref.fa",
						thisGraphOutputDir + "/FASTQ/R_1.fq",
						thisGraphOutputDir + "/FASTQ/R_2.fq",
						BAM_outputDirectory + "/output.bam",
						true
			);

			std::string BAM = thisGraphOutputDir + "/BAM/output.bam";
			assert(Utilities::fileExists(BAM));

			mapper::processBAM BAMprocessor(thisGraphOutputDir);
			std::pair<double, double> IS_estimate = BAMprocessor.estimateInsertSize(BAM, true);

			BAMprocessor.alignReads(BAM, &trueLevels, IS_estimate.first, IS_estimate.second, readAlignmentOutputDir, true);
			// BAMprocessor.alignReads(0, IS_estimate.fBirst, IS_estimate.second);

			assert(trueLevels.get_total_and_correct().first != 0);
			std::cout << "Quality:\n";
			std::cout << "\tGraph: " << trueLevels.get_total_and_correct().first << " " << ((double)trueLevels.get_total_and_correct().second / trueLevels.get_total_and_correct().first) << "\n";
			std::cout << std::flush;

//			 assert(1  == 0 );
		}
	}
	else if(arguments.at("action") == "testAlignments2Chains")
	{
		std::string qualityMatrixFile = (Utilities::fileExists("src/simulator/predefinedQualityMatrices/I101_NA12878.txt") ? "src/simulator/predefinedQualityMatrices/I101_NA12878.txt" : "simulator/predefinedQualityMatrices/I101_NA12878.txt");

		for(unsigned int outerIteration = 1; outerIteration <= 2; outerIteration++)
		{
			for(int switchContigs = 0; switchContigs <= 1; switchContigs++)
			{
				simpleGraphSimulator graphSimulator((bool)switchContigs, pF);
				Graph* g = graphSimulator.getGraph();

				for(int withError = 0; withError <= ((switchContigs == 0) ? 1 : 0); withError++)
				{

					std::string title = "switchContigs" + Utilities::ItoStr(switchContigs) + "_" + "withError" + Utilities::ItoStr(withError);


					std::cout << "Test run: " << title << "\n" << std::flush;

					std::map<std::string, std::vector<mapper::reads::PRGContigBAMAlignment>> simulatedBAMAlignments = graphSimulator.simulateBAMAlignments(
								50,
								qualityMatrixFile,
								101,
								10,
								1,
								(int)withError
					);

					chrono::milliseconds ms_before = chrono::duration_cast< chrono::milliseconds >(
						chrono::system_clock::now().time_since_epoch() );

			//		std::cout << "Start miliseconds: " << ms_before.count() << "\n" << std::flush;

					size_t transformedAlignments = 0;
					size_t totalL_graph = 0;
					size_t totalM_graph = 0;
					size_t totalL_sequence = 0;
					size_t totalM_sequence = 0;
					// std::cout << "Now examine alignment transformation...\n" << std::flush;
					for(std::map<std::string, std::vector<mapper::reads::PRGContigBAMAlignment>>::iterator categoryIt = simulatedBAMAlignments.begin(); categoryIt != simulatedBAMAlignments.end(); categoryIt++)
					{
						std::string category = categoryIt->first;
						// std::cout << "\t" << category << "\n" << std::flush;

						for(unsigned int bamAlignmentI = 0; bamAlignmentI < simulatedBAMAlignments.at(category).size(); bamAlignmentI++)
						{
		//					std::cout << "\t\tbamAlignmentI: " << bamAlignmentI << "\n" << std::flush;
							mapper::reads::PRGContigBAMAlignment al = simulatedBAMAlignments.at(category).at(bamAlignmentI);

							std::cout << "Alignment" << "\n";
							std::cout << "\t" << "al.graph_aligned_levels" << ": " << Utilities::join(Utilities::ItoStr(al.graph_aligned_levels), ", ") << "\n";
							std::cout << "\t" << "al.graph_aligned       " << ": " << al.graph_aligned << "\n";
							std::cout << "\t" << "al.sequence_aligned    " << ": " << al.sequence_aligned << "\n";
							std::cout << "\n" << std::flush;

			//				std::cout << Utilities::join(Utilities::ItoStr(al.graph_aligned_levels), ", ") << "\n";
			//				std::cout << al.graph_aligned << "\n";
			//				std::cout << al.sequence_aligned << "\n";
			//				std::cout << "\n" << std::flush;
			//				assert(1 == 0);

							mapper::reads::verboseSeedChain graphSeed;
							mapper::reads::verboseSeedChain sequenceSeed;

							std::vector<bool> inGraphGapStretch;
							inGraphGapStretch.resize(g->NodesPerLevel.size() - 1, false);
							mapper::processBAM::PRGContigAlignment2Seed(g, al, true, graphSeed, sequenceSeed, inGraphGapStretch);

							std::vector<int> originalLevels;
							std::string originalSequence;
							assert(al.sequence_aligned.length() == al.graph_aligned_levels.size());
							for(unsigned int i = 0; i < al.sequence_aligned.length(); i++)
							{
								if(al.sequence_aligned.at(i) != '_')
								{
									originalLevels.push_back(al.graph_aligned_levels.at(i));
									originalSequence.push_back(al.sequence_aligned.at(i));
								}
							}

							totalL_graph += graphSeed.qualityLevels(originalLevels, originalSequence).first;
							totalM_graph += graphSeed.qualityLevels(originalLevels, originalSequence).second;
							totalL_sequence += sequenceSeed.quality().first;
							totalM_sequence += sequenceSeed.quality().second;
							transformedAlignments++;
						}
					}

					// std::cout << "Transformed alignments: " << transformedAlignments << "\n";

					chrono::milliseconds ms_after= chrono::duration_cast< chrono::milliseconds >(
							chrono::system_clock::now().time_since_epoch()
					);
			//		std::cout << "Stop miliseconds: " << ms_after.count() << "\n" << std::flush;

					chrono::milliseconds ms_duration = ms_after - ms_before;

					// std::cout << "Duration: " << ms_duration.count() << " ms\n" << std::flush;

					double time_per_alignment_ms = (double)transformedAlignments / (double)ms_duration.count();
					double alignments_per_s = time_per_alignment_ms * 1000;

					std::cout << "Alignments per second " << title << ": " << alignments_per_s << "\n" << std::flush;
					std::cout << "Quality:\n";
					std::cout << "\tGraph: " << totalL_graph << " " << ((double)totalM_graph / totalL_graph) << "\n";
					std::cout << "\tSequence: " << totalL_sequence << " " << ((double)totalM_sequence / totalL_sequence) << "\n";
					std::cout << "\n" << std::flush;
				}
			}
		}
	}
	else if(arguments.at("action") == "testChainExtension")
	{
		std::string qualityMatrixFile = (Utilities::fileExists("src/simulator/predefinedQualityMatrices/I101_NA12878.txt") ? "src/simulator/predefinedQualityMatrices/I101_NA12878.txt" : "simulator/predefinedQualityMatrices/I101_NA12878.txt");

		for(unsigned int outerIteration = 1; outerIteration < 20; outerIteration++)
		{
			for(int switchContigs = 0; switchContigs <= 1; switchContigs++)
			{
				simpleGraphSimulator graphSimulator((bool)switchContigs, pF);
				Graph* g = graphSimulator.getGraph();

				mapper::aligner::extensionAligner g_extensionAligner(g);

				for(int withError = 0; withError <= ((switchContigs == 0) ? 1 : 0); withError++)
				{

					std::string title = "switchContigs" + Utilities::ItoStr(switchContigs) + "_" + "withError" + Utilities::ItoStr(withError);


					std::cout << "Test run: " << title << "\n" << std::flush;

					std::map<std::string, std::vector<mapper::reads::PRGContigBAMAlignment>> simulatedBAMAlignments = graphSimulator.simulateBAMAlignments(
								5,
								qualityMatrixFile,
								101,
								10,
								1,
								(int)withError
					);

					chrono::milliseconds ms_before = chrono::duration_cast< chrono::milliseconds >(
						chrono::system_clock::now().time_since_epoch() );

					size_t transformedAlignments = 0;
					size_t totalL_graph = 0;
					size_t totalM_graph = 0;
					size_t totalL_sequence = 0;
					size_t totalM_sequence = 0;
					// std::cout << "Now examine alignment transformation...\n" << std::flush;
					for(std::map<std::string, std::vector<mapper::reads::PRGContigBAMAlignment>>::iterator categoryIt = simulatedBAMAlignments.begin(); categoryIt != simulatedBAMAlignments.end(); categoryIt++)
					{
						std::string category = categoryIt->first;
						// std::cout << "\t" << category << "\n" << std::flush;

						for(unsigned int bamAlignmentI = 0; bamAlignmentI < simulatedBAMAlignments.at(category).size(); bamAlignmentI++)
						{
							mapper::reads::PRGContigBAMAlignment al = simulatedBAMAlignments.at(category).at(bamAlignmentI);

							// unsigned int al_originalLength = al.graph_aligned_levels.size();

							if((bamAlignmentI % 500) == 0)
							{
								std::cout << "Alignment: " << bamAlignmentI << "\n" << std::flush;
							}

							std::string sequence = al.sequence_aligned;
							assert(sequence.length() > 30);

							std::vector<int> originalLevels;
							std::string originalSequence;
							assert(al.sequence_aligned.length() == al.graph_aligned_levels.size());
							for(unsigned int i = 0; i < al.sequence_aligned.length(); i++)
							{
								if(al.sequence_aligned.at(i) != '_')
								{
									originalLevels.push_back(al.graph_aligned_levels.at(i));
									originalSequence.push_back(al.sequence_aligned.at(i));
								}
							}

							al.checkAlignmentConcordanceWithSequence(originalSequence);

							al.removeSequenceCharacters(true, 10);
							al.removeSequenceCharacters(false, 10);

							al.checkAlignmentConcordanceWithSequence(originalSequence);

							// std::cout << "----\n\n" << std::flush;

							mapper::reads::verboseSeedChain graphSeed;
							mapper::reads::verboseSeedChain sequenceSeed;

							std::vector<bool> inGraphGapStretch;
							inGraphGapStretch.resize(g->NodesPerLevel.size() - 1, false);
							mapper::processBAM::PRGContigAlignment2Seed(g, al, true, graphSeed, sequenceSeed, inGraphGapStretch);
							sequenceSeed.checkChainConcordanceWithSequence(originalSequence);

							mapper::reads::verboseSeedChain sequenceSeed_extended = g_extensionAligner.extendSeedChain(originalSequence, sequenceSeed);

							// unsigned int extended_alignment_levels = sequenceSeed_extended.graph_aligned_levels.size();

							std::string extendedSeed_sequence_noGaps = filter(sequenceSeed_extended.sequence_aligned,[](unsigned char c){return (c != '_');});

							if(!(extendedSeed_sequence_noGaps == originalSequence))
							{
								std::cerr << "! (extendedSeed_sequence_noGaps == originalSequence)\n";
								std::cerr << "\t" << "extendedSeed_sequence_noGaps" << ": " << extendedSeed_sequence_noGaps << "\n";
								std::cerr << "\t" << "originalSequence" << ": " << originalSequence << "\n";
								std::cerr << "\t" << "sequenceSeed_extended.graph_aligned" << ": " << sequenceSeed_extended.graph_aligned << "\n";
								std::cerr << "\t" << "sequenceSeed_extended.sequence_aligned" << ": " << sequenceSeed_extended.sequence_aligned << "\n";
								std::cerr << "\n" << std::flush;
							}
							assert(extendedSeed_sequence_noGaps == originalSequence);
							totalL_graph += graphSeed.quality().first;
							totalM_graph += graphSeed.quality().second;
							totalL_sequence += sequenceSeed_extended.quality().first;
							totalM_sequence += sequenceSeed_extended.quality().second;
							transformedAlignments++;
						}
					}

					chrono::milliseconds ms_after= chrono::duration_cast< chrono::milliseconds >(
							chrono::system_clock::now().time_since_epoch()
					);

					chrono::milliseconds ms_duration = ms_after - ms_before;

					double time_per_alignment_ms = (double)transformedAlignments / (double)ms_duration.count();
					double alignments_per_s = time_per_alignment_ms * 1000;

					std::cout << "Alignments per second " << title << ": " << alignments_per_s << "\n" << std::flush;
					std::cout << "Quality:\n";
					std::cout << "\tGraph: " << totalL_graph << " " << ((double)totalM_graph / totalL_graph) << "\n";
					std::cout << "\tSequence: " << totalL_sequence << " " << ((double)totalM_sequence / totalL_sequence) << "\n";
					std::cout << "\n" << std::flush;
				}
			}
		}
	}
	else if(arguments.at("action") == "oneSimulationFromPRG")
	{
		/*
		arguments["action"] = "oneSimulationFromPRG";
		arguments["simulation_outputDirectory"] = "tmp/simulated";
		arguments["graph"] = "C:\\Users\\AlexanderDilthey\\OneDrive\\PRG-BWA\\testData\\TEST_40_from_HLA-P_to_HLA-G\\";
		arguments["simulations_haploidCoverage"] = "10";
		arguments["simulations_withError"] = "1";
		*/

		assert(arguments.count("graph"));
		assert(arguments.count("simulation_outputDirectory"));
		assert(arguments.count("simulations_haploidCoverage"));
		assert(arguments.count("simulations_withError"));

		Utilities::clearDirectory(arguments.at("simulation_outputDirectory"));

		std::string mockGraph = arguments.at("graph") + std::string("/PRG/graph.txt");
		Graph* g = new Graph();
		g->readFromFile(mockGraph);

		std::cout << "Simulator: have graph with " << g->NodesPerLevel.size() << " levels\n";

		double simulations_haploidCoverage = Utilities::StrtoD(arguments.at("simulations_haploidCoverage"));
		bool withError = Utilities::StrtoI(arguments.at("simulations_withError"));

		std::string qualityMatrixFile = (Utilities::fileExists("src/simulator/predefinedQualityMatrices/I101_NA12878.txt") ? "src/simulator/predefinedQualityMatrices/I101_NA12878.txt" : "simulator/predefinedQualityMatrices/I101_NA12878.txt");
		simulator::simulator S(qualityMatrixFile);
		S.simulateFromGraph(g, 1, arguments.at("simulation_outputDirectory"), simulations_haploidCoverage, withError);

	}
	else if(arguments.at("action") == "simulateFromNormalGenome")
	{
		/*
		arguments["action"] = "simulateFromNormalGenome";
		arguments["simulation_outputDirectory"] = "tmp/simulatedNormalGenome";
		arguments["graph"] = "C:\\Users\\AlexanderDilthey\\OneDrive\\PRG-BWA\\output";
		arguments["referenceGenomePath"] = "C:\\Users\\AlexanderDilthey\\OneDrive\\PRG-BWA\\testData\\Reference6Only\\_GRCh38_Chr6.fa";
		arguments["simulations_haploidCoverage"] = "10";
		arguments["simulations_withError"] = "1";
		*/
		
		assert(arguments.count("graph"));
		assert(arguments.count("simulation_outputDirectory"));
		assert(arguments.count("referenceGenomePath"));
		assert(arguments.count("simulations_haploidCoverage"));
		assert(arguments.count("simulations_withError"));

		Utilities::clearDirectory(arguments.at("simulation_outputDirectory"));


		double simulations_haploidCoverage = Utilities::StrtoD(arguments.at("simulations_haploidCoverage"));
		bool withError = Utilities::StrtoI(arguments.at("simulations_withError"));
		int threads = 1;
		if(arguments.count("threads"))
		{
			threads = Utilities::StrtoI(arguments.at("threads"));
		}
		
		std::string qualityMatrixFile = (Utilities::fileExists("src/simulator/predefinedQualityMatrices/I101_NA12878.txt") ? "src/simulator/predefinedQualityMatrices/I101_NA12878.txt" : "simulator/predefinedQualityMatrices/I101_NA12878.txt");
		simulator::simulator S(qualityMatrixFile);
		S.simulateNormalGenome(arguments.at("graph"), arguments.at("referenceGenomePath"), arguments.at("simulation_outputDirectory"), simulations_haploidCoverage, withError, threads);
	}
	else
	{
		throw std::runtime_error("Invalid --action: " + arguments.at("action"));
	}

	return 0;
}

