/*
 * simulator.cpp
 *
 *  Created on: 16.09.2015
 *      Author: AlexanderDilthey
 */

#include "simulator.h"

#include "../Utilities.h"

#include <assert.h>
#include <fstream>
#include <map>
#include <vector>
#include <stdexcept>
#include <exception>
#include <omp.h>

namespace simulator {


simulator::simulator(std::string qualityMatrixFile_, int read_length_, double insertSize_mean_, double insertSize_sd_) {
	read_length = read_length_;
	insertSize_mean = insertSize_mean_;
	insertSize_sd= insertSize_sd_;
	qualityMatrixFile = qualityMatrixFile_;
	rS = new readSimulator(qualityMatrixFile, read_length);
}

void simulator::simulateFromHaploidAlignedSequence(const std::string& S, const std::string& outputDirectory, double haploidCoverage, bool withError) const
{
	for(unsigned int i = 0; i  < S.length(); i++)
	{
		std::string c = S.substr(i, 1);
		bool valid = (
			(c == "A") ||
			(c == "C") ||
			(c == "G") ||
			(c == "T") ||
			(c == "_")
		);
		
		if(! valid)
		{
			std::cerr << "Unexpected character: -" << c << "-\n" << std::flush;
		}
		assert(valid);
	}

	std::string parametersFile = outputDirectory + "/parameters.txt";
	std::ofstream parametersStream;
	parametersStream.open(parametersFile.c_str());
	parametersStream << "qualityMatrixFile: " << qualityMatrixFile << "\n";
	parametersStream << "read_length: " << read_length << "\n";
	parametersStream << "insertSize_mean: " << insertSize_mean << "\n";
	parametersStream << "insertSize_sd: " << insertSize_sd << "\n";
	parametersStream << "haploidCoverage: " << haploidCoverage << "\n";
	parametersStream << "withError: " << withError << "\n";
	parametersStream << "rS average error rates: " << rS->averageErrorRate_R1_R2().first << "\t" <<  rS->averageErrorRate_R1_R2().second <<  "\n";
	parametersStream.close();


	std::cout << "simulator::simulateFromHaploidAlignedSequence(..): Start read simulation.\n" << std::flush;
	std::cout << "\t" << "Sequence length: " << S.length() << "\n" << std::flush;

	std::vector<oneReadPair> reads = rS->simulate_paired_reads_from_string(S, haploidCoverage, insertSize_mean, insertSize_sd, ! withError, "H1");

	std::string fn_output_FASTQ_r1 = outputDirectory + "/R_1.fq";
	std::string fn_output_FASTQ_r2 = outputDirectory + "/R_2.fq";

	std::string fn_output_levels_r1 = outputDirectory + "/R_1.levels";
	std::string fn_output_levels_r2 = outputDirectory + "/R_2.levels";

	std::ofstream stream_FASTQ_r1;
	stream_FASTQ_r1.open(fn_output_FASTQ_r1.c_str());
	assert(stream_FASTQ_r1.is_open());

	std::ofstream stream_FASTQ_r2;
	stream_FASTQ_r2.open(fn_output_FASTQ_r2.c_str());
	assert(stream_FASTQ_r2.is_open());

	std::ofstream stream_levels_r1;
	stream_levels_r1.open(fn_output_levels_r1.c_str());
	assert(stream_levels_r1.is_open());

	std::ofstream stream_levels_r2;
	stream_levels_r2.open(fn_output_levels_r2.c_str());
	assert(stream_levels_r2.is_open());

	/*
	auto fullAlignment_coordinates_edgePath_onlyNonGaps = [](std::vector<int> fullAlignment_coordinates_edgePath, std::string underlyingEdgeLabels) -> std::vector<int> {
		assert(fullAlignment_coordinates_edgePath.size() == underlyingEdgeLabels.size());
		std::vector<int> forReturn;
		for(unsigned int i = 0; i < fullAlignment_coordinates_edgePath.size(); i++)
		{
			std::string c = underlyingEdgeLabels.substr(i, 1);
			if(c != "_")
			{
				forReturn.push_back();
			}
		}
		return forReturn;
	};
	*/
	
	auto print_one_readPair = [&] (const oneReadPair& rP) -> void
	{
		stream_levels_r1 << "@" << rP.reads.first.name << "\n";
		// std::vector<int> r1_fullAlignment_coordinates_edgePath_onlyNonGaps = fullAlignment_coordinates_edgePath_onlyNonGaps(rP.reads.first.fullAlignment_coordinates_edgePath, rP.reads.first.underlyingEdgeLabels);
		stream_levels_r1 << Utilities::join(Utilities::ItoStr(rP.reads.first.coordinates_edgePath), " ") << "\n";
		stream_levels_r1 << rP.reads.first.underlyingEdgeLabels << "\n";
		stream_levels_r1 << Utilities::join(Utilities::ItoStr(rP.reads.first.fullAlignment_coordinates_edgePath), " ") << "\n";
		stream_levels_r1 << rP.reads.first.fullAlignment_underlyingEdgeLabels << "\n";
		stream_levels_r1 << rP.reads.first.fullAlignment_sequence << "\n";

		
		
		stream_levels_r2 << "@" << rP.reads.second.name << "\n";
		//std::vector<int> r2_fullAlignment_coordinates_edgePath_onlyNonGaps = fullAlignment_coordinates_edgePath_onlyNonGaps(rP.reads.second.fullAlignment_coordinates_edgePath, rP.reads.second.underlyingEdgeLabels);		
		stream_levels_r2 << Utilities::join(Utilities::ItoStr(rP.reads.second.coordinates_edgePath), " ") << "\n";
		stream_levels_r2 << rP.reads.second.underlyingEdgeLabels << "\n";
		stream_levels_r2 << Utilities::join(Utilities::ItoStr(rP.reads.second.fullAlignment_coordinates_edgePath), " ") << "\n";
		stream_levels_r2 << rP.reads.second.fullAlignment_underlyingEdgeLabels << "\n";
		stream_levels_r2 << rP.reads.second.fullAlignment_sequence << "\n";
		
		
		stream_FASTQ_r1 << "@" << rP.reads.first.name << "\n";
		stream_FASTQ_r1 << rP.reads.first.sequence << "\n";
		stream_FASTQ_r1 << "+" << "\n";
		stream_FASTQ_r1 << rP.reads.first.quality << "\n";

		stream_FASTQ_r2 << "@" << rP.reads.second.name << "\n";
		stream_FASTQ_r2 << rP.reads.second.sequence << "\n";
		stream_FASTQ_r2 << "+" << "\n";
		stream_FASTQ_r2 << rP.reads.second.quality << "\n";

		/*
		stream_levels_r1 << "@" << rP.reads.first.name << "\n";
		stream_levels_r1 << Utilities::join(Utilities::ItoStr(rP.reads.first.coordinates_string), " ") << "\n";

		stream_levels_r2 << "@" << rP.reads.second.name << "\n";
		stream_levels_r2 << Utilities::join(Utilities::ItoStr(rP.reads.second.coordinates_string), " ") << "\n";
		*/
//		if(rP.reads.first.fullAlignment_underlyingEdgeLabels != rP.reads.first.fullAlignment_sequence)
//		{
//			if(rP.reads.first.fullAlignment_sequence.find("_") != std::string::npos)
//			{
//				std::cerr << "GOOD 1!\n";
//				std::cerr << "\t" << Utilities::join(Utilities::ItoStr(rP.reads.first.coordinates_edgePath), " ") << "\n";
//				std::cerr << "\t" << rP.reads.first.fullAlignment_underlyingEdgeLabels << "\n";
//				std::cerr << "\t" << rP.reads.first.fullAlignment_sequence << "\n" << std::flush;
//			}
//		}

//		if(rP.reads.second.fullAlignment_underlyingEdgeLabels != rP.reads.second.fullAlignment_sequence)
//		{
//			std::cerr << "GOOD 2!\n";
//			std::cerr << "\t" << rP.reads.second.fullAlignment_sequence << "\n"  << std::flush;
//		}
	};

	for(unsigned int rI = 0; rI < reads.size(); rI++)
	{
		print_one_readPair(reads.at(rI));
	}
}

void simulator::simulateFromDiploidHaplotypes(const std::string& h1, const std::string& h2, const std::string& outputDirectory, double haploidCoverage, bool withError) const
{

	std::string parametersFile = outputDirectory + "/parameters.txt";
	std::ofstream parametersStream;
	parametersStream.open(parametersFile.c_str());
	parametersStream << "qualityMatrixFile: " << qualityMatrixFile << "\n";
	parametersStream << "read_length: " << read_length << "\n";
	parametersStream << "insertSize_mean: " << insertSize_mean << "\n";
	parametersStream << "insertSize_sd: " << insertSize_sd << "\n";
	parametersStream << "haploidCoverage: " << haploidCoverage << "\n";
	parametersStream << "withError: " << withError << "\n";
	parametersStream << "rS average error rates: " << rS->averageErrorRate_R1_R2().first << "\t" <<  rS->averageErrorRate_R1_R2().second <<  "\n";
	parametersStream.close();


	std::cout << "simulator::simulateFromDiploidHaplotypes(..): Start read simulation.\n" << std::flush;
	std::cout << "\t" << "H1 length: " << h1.length() << "\n" << std::flush;
	std::cout << "\t" << "H2 length: " << h2.length() << "\n" << std::flush;

	std::vector<oneReadPair> reads_h1 = rS->simulate_paired_reads_from_string(h1, haploidCoverage, insertSize_mean, insertSize_sd, ! withError, "H1");
	std::vector<oneReadPair> reads_h2 = rS->simulate_paired_reads_from_string(h2, haploidCoverage, insertSize_mean, insertSize_sd, ! withError, "H2");

	std::string fn_output_FASTQ_r1 = outputDirectory + "/R_1.fq";
	std::string fn_output_FASTQ_r2 = outputDirectory + "/R_2.fq";

	std::string fn_output_levels_r1 = outputDirectory + "/R_1.levels";
	std::string fn_output_levels_r2 = outputDirectory + "/R_2.levels";

	std::ofstream stream_FASTQ_r1;
	stream_FASTQ_r1.open(fn_output_FASTQ_r1.c_str());
	assert(stream_FASTQ_r1.is_open());

	std::ofstream stream_FASTQ_r2;
	stream_FASTQ_r2.open(fn_output_FASTQ_r2.c_str());
	assert(stream_FASTQ_r2.is_open());

	std::ofstream stream_levels_r1;
	stream_levels_r1.open(fn_output_levels_r1.c_str());
	assert(stream_levels_r1.is_open());

	std::ofstream stream_levels_r2;
	stream_levels_r2.open(fn_output_levels_r2.c_str());
	assert(stream_levels_r2.is_open());

	auto print_one_readPair = [&] (const oneReadPair& rP) -> void
	{
		stream_FASTQ_r1 << "@" << rP.reads.first.name << "\n";
		stream_FASTQ_r1 << rP.reads.first.sequence << "\n";
		stream_FASTQ_r1 << "+" << "\n";
		stream_FASTQ_r1 << rP.reads.first.quality << "\n";

		stream_FASTQ_r2 << "@" << rP.reads.second.name << "\n";
		stream_FASTQ_r2 << rP.reads.second.sequence << "\n";
		stream_FASTQ_r2 << "+" << "\n";
		stream_FASTQ_r2 << rP.reads.second.quality << "\n";

		stream_levels_r1 << "@" << rP.reads.first.name << "\n";
		stream_levels_r1 << Utilities::join(Utilities::ItoStr(rP.reads.first.coordinates_string), " ") << "\n";

		stream_levels_r2 << "@" << rP.reads.second.name << "\n";
		stream_levels_r2 << Utilities::join(Utilities::ItoStr(rP.reads.second.coordinates_string), " ") << "\n";

//		if(rP.reads.first.fullAlignment_underlyingEdgeLabels != rP.reads.first.fullAlignment_sequence)
//		{
//			if(rP.reads.first.fullAlignment_sequence.find("_") != std::string::npos)
//			{
//				std::cerr << "GOOD 1!\n";
//				std::cerr << "\t" << Utilities::join(Utilities::ItoStr(rP.reads.first.coordinates_edgePath), " ") << "\n";
//				std::cerr << "\t" << rP.reads.first.fullAlignment_underlyingEdgeLabels << "\n";
//				std::cerr << "\t" << rP.reads.first.fullAlignment_sequence << "\n" << std::flush;
//			}
//		}

//		if(rP.reads.second.fullAlignment_underlyingEdgeLabels != rP.reads.second.fullAlignment_sequence)
//		{
//			std::cerr << "GOOD 2!\n";
//			std::cerr << "\t" << rP.reads.second.fullAlignment_sequence << "\n"  << std::flush;
//		}
	};

	for(unsigned int rI = 0; rI < reads_h1.size(); rI++)
	{
		print_one_readPair(reads_h1.at(rI));
	}

	for(unsigned int rI = 0; rI < reads_h2.size(); rI++)
	{
		print_one_readPair(reads_h2.at(rI));
	}
}

void simulator::simulateFromGraph(Graph* g, int simulatedGraphGenomes, std::string outputDirectory, double haploidCoverage, bool withError) const {


	std::string parametersFile = outputDirectory + "/parameters.txt";
	std::ofstream parametersStream;
	parametersStream.open(parametersFile.c_str());
	parametersStream << "Graph: " << g->filename_last_read << "\n";
	parametersStream << "simulatedGraphGenomes: " << simulatedGraphGenomes << "\n";
	parametersStream << "qualityMatrixFile: " << qualityMatrixFile << "\n";
	parametersStream << "read_length: " << read_length << "\n";
	parametersStream << "insertSize_mean: " << insertSize_mean << "\n";
	parametersStream << "insertSize_sd: " << insertSize_sd << "\n";
	parametersStream << "haploidCoverage: " << haploidCoverage << "\n";
	parametersStream << "withError: " << withError << "\n";
	parametersStream << "rS average error rates: " << rS->averageErrorRate_R1_R2().first << "\t" <<  rS->averageErrorRate_R1_R2().second <<  "\n";
	parametersStream.close();

	std::string fn_output_FASTQ_r1 = outputDirectory + "/R_1.fq";
	std::string fn_output_FASTQ_r2 = outputDirectory + "/R_2.fq";

	std::string fn_output_levels_r1 = outputDirectory + "/R_1.levels";
	std::string fn_output_levels_r2 = outputDirectory + "/R_2.levels";

	std::ofstream stream_FASTQ_r1;
	stream_FASTQ_r1.open(fn_output_FASTQ_r1.c_str());
	assert(stream_FASTQ_r1.is_open());

	std::ofstream stream_FASTQ_r2;
	stream_FASTQ_r2.open(fn_output_FASTQ_r2.c_str());
	assert(stream_FASTQ_r2.is_open());

	std::ofstream stream_levels_r1;
	stream_levels_r1.open(fn_output_levels_r1.c_str());
	assert(stream_levels_r1.is_open());

	std::ofstream stream_levels_r2;
	stream_levels_r2.open(fn_output_levels_r2.c_str());
	assert(stream_levels_r2.is_open());

	auto print_one_readPair = [&] (const oneReadPair& rP) -> void
	{
		stream_FASTQ_r1 << "@" << rP.reads.first.name << "\n";
		stream_FASTQ_r1 << rP.reads.first.sequence << "\n";
		stream_FASTQ_r1 << "+" << "\n";
		stream_FASTQ_r1 << rP.reads.first.quality << "\n";

		stream_FASTQ_r2 << "@" << rP.reads.second.name << "\n";
		stream_FASTQ_r2 << rP.reads.second.sequence << "\n";
		stream_FASTQ_r2 << "+" << "\n";
		stream_FASTQ_r2 << rP.reads.second.quality << "\n";

		stream_levels_r1 << "@" << rP.reads.first.name << "\n";
		stream_levels_r1 << Utilities::join(Utilities::ItoStr(rP.reads.first.coordinates_edgePath), " ") << "\n";
		stream_levels_r1 << rP.reads.first.underlyingEdgeLabels << "\n";
		stream_levels_r1 << Utilities::join(Utilities::ItoStr(rP.reads.first.fullAlignment_coordinates_edgePath), " ") << "\n";
		stream_levels_r1 << rP.reads.first.fullAlignment_underlyingEdgeLabels << "\n";
		stream_levels_r1 << rP.reads.first.fullAlignment_sequence << "\n";

		stream_levels_r2 << "@" << rP.reads.second.name << "\n";
		stream_levels_r2 << Utilities::join(Utilities::ItoStr(rP.reads.second.coordinates_edgePath), " ") << "\n";
		stream_levels_r2 << rP.reads.second.underlyingEdgeLabels << "\n";
		stream_levels_r2 << Utilities::join(Utilities::ItoStr(rP.reads.second.fullAlignment_coordinates_edgePath), " ") << "\n";
		stream_levels_r2 << rP.reads.second.fullAlignment_underlyingEdgeLabels << "\n";
		stream_levels_r2 << rP.reads.second.fullAlignment_sequence << "\n";

//		if(rP.reads.first.fullAlignment_underlyingEdgeLabels != rP.reads.first.fullAlignment_sequence)
//		{
//			if(rP.reads.first.fullAlignment_sequence.find("_") != std::string::npos)
//			{
//				std::cerr << "GOOD 1!\n";
//				std::cerr << "\t" << Utilities::join(Utilities::ItoStr(rP.reads.first.coordinates_edgePath), " ") << "\n";
//				std::cerr << "\t" << rP.reads.first.fullAlignment_underlyingEdgeLabels << "\n";
//				std::cerr << "\t" << rP.reads.first.fullAlignment_sequence << "\n" << std::flush;
//			}
//		}

//		if(rP.reads.second.fullAlignment_underlyingEdgeLabels != rP.reads.second.fullAlignment_sequence)
//		{
//			std::cerr << "GOOD 2!\n";
//			std::cerr << "\t" << rP.reads.second.fullAlignment_sequence << "\n"  << std::flush;
//		}
	};

	size_t combined_read_pairs = 0;


	for(int genomeI = 0; genomeI < simulatedGraphGenomes; genomeI++)
	{
		diploidEdgePointerPath simulatedGraphGenome = g->simulateRandomDiploidPath();
		std::vector<oneReadPair> reads_h1 = rS->simulate_paired_reads_from_edgePath(simulatedGraphGenome.h1, haploidCoverage, insertSize_mean, insertSize_sd, ! withError, "PRG_"+Utilities::ItoStr(genomeI)+"h1_");
		std::vector<oneReadPair> reads_h2 = rS->simulate_paired_reads_from_edgePath(simulatedGraphGenome.h2, haploidCoverage, insertSize_mean, insertSize_sd, ! withError, "PRG_"+Utilities::ItoStr(genomeI)+"h2_");


		for(size_t i = 0; i < reads_h1.size(); i++)
		{
			print_one_readPair(reads_h1.at(i));
		}

		for(size_t i = 0; i < reads_h2.size(); i++)
		{
			print_one_readPair(reads_h2.at(i));
		}

		combined_read_pairs += (reads_h1.size() + reads_h2.size());
	}

	stream_FASTQ_r1.close();
	stream_FASTQ_r2.close();
	stream_levels_r1.close();
	stream_levels_r2.close();

	std::cout << "simulator::simulateFromGraph(..): Printed " << combined_read_pairs << " reads pairs.\n" << std::flush;
}

void simulator::simulateNormalGenome(std::string graphDir, std::string referenceGenomePath, std::string outputDirectory, double haploidCoverage, bool withError, unsigned int threads)
{
	std::string sequencesFile = graphDir + "/sequences.txt";
	assert(Utilities::fileExists(sequencesFile));

	std::string parametersFile = outputDirectory + "/parameters.txt";
	std::ofstream parametersStream;
	parametersStream.open(parametersFile.c_str());
	parametersStream << "referenceGenomePath: " << referenceGenomePath << "\n";
	parametersStream << "qualityMatrixFile: " << qualityMatrixFile << "\n";
	parametersStream << "read_length: " << read_length << "\n";
	parametersStream << "insertSize_mean: " << insertSize_mean << "\n";
	parametersStream << "insertSize_sd: " << insertSize_sd << "\n";
	parametersStream << "haploidCoverage: " << haploidCoverage << "\n";
	parametersStream << "withError: " << withError << "\n";
	parametersStream << "rS average error rates: " << rS->averageErrorRate_R1_R2().first << "\t" <<  rS->averageErrorRate_R1_R2().second <<  "\n";
	parametersStream << "threads: " << threads << "\n";
	parametersStream.close();

	std::cout << "simulator::simulateNormalGenome(..): Reading reference genome..\n" << std::flush;
	std::map<std::string, std::string> referenceGenome = Utilities::readFASTA(referenceGenomePath);
	std::cout << "\tDone. Have " << referenceGenome.size() << " chromosomes.\n\n" << std::flush;

	std::ifstream sequencesStream;
	sequencesStream.open(sequencesFile.c_str());
	assert(sequencesStream.is_open());

	std::cout << "simulator::simulateNormalGenome(..): Block PRG areas.\n" << std::flush;
	std::string headerLine;
	assert(sequencesStream.good());
	std::getline(sequencesStream, headerLine);
	Utilities::eraseNL(headerLine);
	std::vector<std::string> headerFields = Utilities::split(headerLine, "\t");

	std::string line;
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

		if(L.at("Chr") != "")
		{
			// in reference genome
			if(! referenceGenome.count(L.at("Chr")))
			{
				throw std::runtime_error("FASTA ID " +  L.at("Chr") + " not in reference genome "+referenceGenomePath);
			}
			assert(referenceGenome.count(L.at("Chr")));

			int start_1based = Utilities::StrtoI(L.at("Start_1based"));
			int stop_1based = Utilities::StrtoI(L.at("Stop_1based"));
			assert(stop_1based > start_1based);

			assert(start_1based >= 1);
			assert(stop_1based <= (int)referenceGenome.at(L.at("Chr")).length());

			for(int i = start_1based; i <= stop_1based; i++)
			{
				referenceGenome.at(L.at("Chr")).at(i - 1) = 'N';
			}

			std::cout << "\tBlocking " << L.at("Chr") << " from " << start_1based << " to " << stop_1based << " (1-based).\n" << std::flush;
		}
	}

	std::vector<std::string> chromosomeIDs;
	for(std::map<std::string, std::string>::iterator refGenomeIt = referenceGenome.begin(); refGenomeIt != referenceGenome.end(); refGenomeIt++)
	{
		chromosomeIDs.push_back(refGenomeIt->first);
	}

	std::cout << "simulator::simulateNormalGenome(..): Start read simulation.\n" << std::flush;


	for(unsigned int chromosomeI = 0; chromosomeI < chromosomeIDs.size(); chromosomeI++)
	{
		std::string chromosomeID = chromosomeIDs.at(chromosomeI);
		std::string fn_output_prefix = outputDirectory + "/R_Chr" + Utilities::ItoStr(chromosomeI);

		rS->simulate_paired_reads_from_string_mt_immediateOutput(referenceGenome.at(chromosomeID), fn_output_prefix, 2 * haploidCoverage, insertSize_mean, insertSize_sd, ! withError, threads, "REF");
	}

	sequencesStream.close();
}

simulator::~simulator() {
	delete(rS);
}


} /* namespace simulator */
