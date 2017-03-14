/*
 * readSimulator.cpp
 *
 *  Created on: 21.05.2013
 *      Author: AlexanderDilthey
 */

#include <omp.h>
#include <algorithm>
#include <utility>
#include <fstream>
#include <cmath>

#include "readSimulator.h"
#include "../Utilities.h"


#include <boost/random.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/random_device.hpp>

namespace simulator {

std::string readName_field_separator = "|||";

double readSimulator::averageErrorRate(std::vector<std::map<char, double> > q_freq, std::vector<std::map<char, double> > error_freq_conditional_q)
{
	unsigned int readLength = q_freq.size();
	assert(error_freq_conditional_q.size() == readLength);

	double forReturn = 0;

	for(unsigned int positionInRead = 0; positionInRead < readLength; positionInRead++)
	{
		double thisPos_errorRate = 0;
		assert(q_freq.at(positionInRead).size() == error_freq_conditional_q.at(positionInRead).size());
		Utilities::check_map_is_normalized(q_freq.at(positionInRead));
		for(std::map<char, double>::iterator qIt = q_freq.at(positionInRead).begin(); qIt != q_freq.at(positionInRead).end(); qIt++)
		{
			char quality = qIt->first;
			double p_quality = qIt->second;
			double p_error = 1 - error_freq_conditional_q.at(positionInRead).at(quality);
			assert(p_error >= 0);
			assert(p_error <= 1);
			thisPos_errorRate += p_quality * p_error;
		}

		forReturn += thisPos_errorRate;
	}

	forReturn = forReturn / (double) readLength;

	return forReturn;
}

readSimulator::readSimulator(std::string qualityMatrixFile, unsigned int readLength, bool interpolateLength_, char removeUpperBaseQualityIndices, char additional_2ndRead_removeUpperBaseQualityIndices) {

	quiet = true;

	// attention - these are chars - make sure that there are no overflows
	assert(removeUpperBaseQualityIndices >= 0);
	assert(removeUpperBaseQualityIndices <= 20);
	assert(additional_2ndRead_removeUpperBaseQualityIndices >= 0);
	assert(additional_2ndRead_removeUpperBaseQualityIndices <= 20);
	
	std::ifstream matrixStream;
	matrixStream.open (qualityMatrixFile.c_str(), std::ios::in);

	interpolateLength = interpolateLength_;


	read_quality_frequencies.resize(readLength);
	read_quality_correctness.resize(readLength);
	read_INDEL_freq.resize(readLength, 0);

	std::map< int, std::vector<std::map<char, double> > > read_quality_frequencies_perLength;
	std::map< int, std::vector<std::map<char, double> > > read_quality_correctness_perLength;
	std::map< int, std::vector<double> > read_INDEL_freq_perLength;
	std::map< int, std::vector<double> > read_nonINDEL_freq_perLength;

	if(matrixStream.is_open())
	{
		std::string line;
		std::string currentIdentifier;
		std::string currentSequence;

		size_t lineCounter = 0;

		while(matrixStream.good())
		{
			std::getline(matrixStream, line);
			lineCounter++;

			Utilities::eraseNL(line);
			if(line.length() == 0)
				continue;

			std::vector<std::string> line_fields = Utilities::split(line, "\t");

			if(lineCounter == 1)
			{
				if(!(line_fields.size() == 6))
				{
					throw std::runtime_error("readSimulator::readSimulator(): file "+qualityMatrixFile+", expect 6 fields, but have "+Utilities::ItoStr(line_fields.size())+"; line "+Utilities::ItoStr(lineCounter)+".\n"+line+"\n");
				}
				assert(line_fields.at(0) == "readLength");
				assert(line_fields.at(1) == "qualityScore");
				assert(line_fields.at(2) == "positionInRead");
				assert(line_fields.at(3) == "N");
				assert(line_fields.at(4) == "ExpectedCorrect");
				assert(line_fields.at(5) == "EmpiricalCorrect");
			}
			else
			{
				int thisLine_readLength = Utilities::StrtoI(line_fields.at(0));

				if(read_quality_frequencies_perLength.count(thisLine_readLength) == 0)
				{
					read_quality_frequencies_perLength[thisLine_readLength].resize(thisLine_readLength);
					read_quality_correctness_perLength[thisLine_readLength].resize(thisLine_readLength);
					read_INDEL_freq_perLength[thisLine_readLength].resize(thisLine_readLength, 0);
					read_nonINDEL_freq_perLength[thisLine_readLength].resize(thisLine_readLength, 0);
				}

				char quality = line_fields.at(1).at(0);
				unsigned int positionInRead = Utilities::StrtoI(line_fields.at(2));
				int datapoints = Utilities::StrtoI(line_fields.at(3));

				assert(positionInRead < readLength);


				if(quality == 0)
				{
					read_INDEL_freq_perLength.at(thisLine_readLength).at(positionInRead) += datapoints;
				}
				else
				{
					read_nonINDEL_freq_perLength.at(thisLine_readLength).at(positionInRead) += datapoints;

					if(read_quality_frequencies_perLength.at(thisLine_readLength).at(positionInRead).count(quality) == 0)
					{
						read_quality_frequencies_perLength.at(thisLine_readLength).at(positionInRead)[quality] = 0;
					}
					read_quality_frequencies_perLength.at(thisLine_readLength).at(positionInRead)[quality] += datapoints;

					if(!(read_quality_correctness_perLength.at(thisLine_readLength).at(positionInRead).count(quality) == 0))
					{
						throw std::runtime_error("readSimulator::readSimulator(): file "+qualityMatrixFile+", double-defined quality. Line "+Utilities::ItoStr(lineCounter)+".\n"+line+"\n");
					}
					assert(read_quality_correctness_perLength.at(thisLine_readLength).at(positionInRead).count(quality) == 0);
					read_quality_correctness_perLength.at(thisLine_readLength).at(positionInRead)[quality] += Utilities::StrtoD(line_fields.at(5));
				}
			}
		}
		matrixStream.close();

		std::vector<double> read_nonINDEL_freq;


		if(read_quality_correctness_perLength.count(readLength))
		{
			read_quality_frequencies = read_quality_frequencies_perLength.at(readLength);
			read_quality_correctness = read_quality_correctness_perLength.at(readLength);
			read_INDEL_freq = read_INDEL_freq_perLength.at(readLength);
			read_nonINDEL_freq = read_nonINDEL_freq_perLength.at(readLength);
		}
		else
		{
			assert(interpolateLength);

			assert(read_quality_frequencies_perLength.size() > 0);
			int closestNeighbour_which = -1;
			int closestNeighbour_distance = -1;

			for(std::map< int, std::vector<std::map<char, double> > >::iterator correctnessIt = read_quality_correctness_perLength.begin(); correctnessIt != read_quality_correctness_perLength.end(); correctnessIt++ )
			{
				int thisDist = abs(correctnessIt->first - (int)readLength);
				if((correctnessIt == read_quality_correctness_perLength.begin()) || (thisDist < closestNeighbour_distance))
				{
					closestNeighbour_which = correctnessIt->first;
					closestNeighbour_distance = thisDist;
				}
			}
			assert(closestNeighbour_distance >= 0);

			read_quality_frequencies.resize(readLength);
			read_quality_correctness.resize(readLength);
			read_INDEL_freq.resize(readLength);
			read_nonINDEL_freq.resize(readLength);

			for(unsigned int posInRead = 0; posInRead < readLength; posInRead++)
			{
				double fraction = (double)posInRead / (double) readLength;
				double idx_target = fraction * (closestNeighbour_which - 1);
				int idx_target_int = round(idx_target);
				assert(idx_target_int >= 0);
				assert(idx_target_int < (int)read_quality_correctness_perLength.at(closestNeighbour_which).size());

				read_quality_frequencies.at(posInRead) = read_quality_frequencies_perLength.at(closestNeighbour_which).at(idx_target_int);
				read_quality_correctness.at(posInRead) = read_quality_correctness_perLength.at(closestNeighbour_which).at(idx_target_int);
				read_INDEL_freq.at(posInRead) = read_INDEL_freq_perLength.at(closestNeighbour_which).at(idx_target_int);
				read_nonINDEL_freq.at(posInRead) = read_nonINDEL_freq_perLength.at(closestNeighbour_which).at(idx_target_int);
			}
		}

		for(unsigned int i = 0; i < readLength; i++)
		{
			assert(read_quality_correctness.at(i).size() > 0);
			assert(read_quality_frequencies.at(i).size() > 0);
			read_quality_frequencies.at(i) = Utilities::normalize_map(read_quality_frequencies.at(i));

			double total_data_position = read_INDEL_freq.at(i) + read_nonINDEL_freq.at(i);
			assert(total_data_position != 0);
			// std::cout << i << "\t" << read_INDEL_freq.at(i) << "\t" << total_data_position << "\n" << std::flush;
			read_INDEL_freq.at(i) = read_INDEL_freq.at(i) / total_data_position;
			if(read_INDEL_freq.at(i) == 0)
			{
				read_INDEL_freq.at(i) = 1e-4;
			}
		}
		
		read_quality_frequencies_2nd = read_quality_frequencies;
		read_quality_correctness_2nd = read_quality_correctness;
		read_INDEL_freq_2nd = read_INDEL_freq;		
		

		if(removeUpperBaseQualityIndices || additional_2ndRead_removeUpperBaseQualityIndices)
		{
			if(! quiet)
			{
				std::cout << "Before adjusting base qualities by " << (int)removeUpperBaseQualityIndices << " and " << (int)additional_2ndRead_removeUpperBaseQualityIndices << ", have average error:\n";
				std::cout << "\tR1: " << averageErrorRate(read_quality_frequencies, read_quality_correctness) << "\n";
				std::cout << "\tR2: " << averageErrorRate(read_quality_frequencies_2nd, read_quality_correctness_2nd) << "\n";
				std::cout << std::flush;
			}
		}

		if(removeUpperBaseQualityIndices)
		{
			for(unsigned int positionInRead = 0; positionInRead < readLength; positionInRead++)
			{
				std::map<char, double> new_qualityFrequencies;
				std::map<char, double> new_correctness;
				std::vector<char> existingQualities;
				for(std::map<char, double>::iterator existingFrequencyIt = read_quality_frequencies.at(positionInRead).begin(); existingFrequencyIt != read_quality_frequencies.at(positionInRead).end(); existingFrequencyIt++ )
				{
					existingQualities.push_back(existingFrequencyIt->first);
				}
				std::sort(existingQualities.begin(), existingQualities.end());
				if(existingQualities.size() > 1)
				{
					assert(existingQualities.at(0) <= existingQualities.at(1));
				}
				assert(removeUpperBaseQualityIndices >= 0);
				assert(removeUpperBaseQualityIndices < (int)(existingQualities.size()-1));
				for(int baseQualityI = 0; baseQualityI < (int)(existingQualities.size() - removeUpperBaseQualityIndices); baseQualityI++)
				{
					char qualityCharacter = existingQualities.at(baseQualityI);
					double f = read_quality_frequencies.at(positionInRead).at(qualityCharacter);
					new_qualityFrequencies[qualityCharacter] = f;
					new_correctness[qualityCharacter] =  read_quality_correctness.at(positionInRead).at(qualityCharacter);
				}
				new_qualityFrequencies = Utilities::normalize_map(new_qualityFrequencies);
				read_quality_frequencies.at(positionInRead) = new_qualityFrequencies;
				read_quality_correctness.at(positionInRead) = new_correctness;
			}
		}


		if(removeUpperBaseQualityIndices || additional_2ndRead_removeUpperBaseQualityIndices)
		{
			char shift_read_2 = removeUpperBaseQualityIndices + additional_2ndRead_removeUpperBaseQualityIndices;
			
			for(unsigned int positionInRead = 0; positionInRead < readLength; positionInRead++)
			{
				std::map<char, double> new_qualityFrequencies;
				std::map<char, double> new_correctness;

				std::vector<char> existingQualities;
				for(std::map<char, double>::iterator existingFrequencyIt = read_quality_frequencies_2nd.at(positionInRead).begin(); existingFrequencyIt != read_quality_frequencies_2nd.at(positionInRead).end(); existingFrequencyIt++ )
				{
					existingQualities.push_back(existingFrequencyIt->first);
				}
				std::sort(existingQualities.begin(), existingQualities.end());
				if(existingQualities.size() > 1)
				{
					assert(existingQualities.at(0) <= existingQualities.at(1));
				}
				assert(shift_read_2 >= 0);
				assert(shift_read_2 < (int)(existingQualities.size()-1));
				for(int baseQualityI = 0; baseQualityI < (int)(existingQualities.size() - shift_read_2); baseQualityI++)
				{
					char qualityCharacter = existingQualities.at(baseQualityI);
					double f = read_quality_frequencies_2nd.at(positionInRead).at(qualityCharacter);
					new_qualityFrequencies[qualityCharacter] = f;
					new_correctness[qualityCharacter] =  read_quality_correctness_2nd.at(positionInRead).at(qualityCharacter);
				}

				new_qualityFrequencies = Utilities::normalize_map(new_qualityFrequencies);
				read_quality_frequencies_2nd.at(positionInRead) = new_qualityFrequencies;
				read_quality_correctness_2nd.at(positionInRead) = new_correctness;

			}
		}

	}
	else
	{
		throw std::runtime_error("readSimulator::readSimulator(): Cannot open file "+qualityMatrixFile);
	}

	for(unsigned i = 0; i < readLength; i++)
	{
		Utilities::check_map_is_normalized(read_quality_frequencies.at(i));
		Utilities::check_map_is_normalized(read_quality_frequencies_2nd.at(i));
	}
	
	if(1 == 1)
	{
		if(! quiet)
		{
			std::cout << "Summary error probabilities:\n";
			std::cout << "After adjusting base qualities by " << (int)removeUpperBaseQualityIndices << " and " << (int)additional_2ndRead_removeUpperBaseQualityIndices << ", have average error:\n";
			std::cout << "\tR1: " << averageErrorRate(read_quality_frequencies, read_quality_correctness) << "\n";
			std::cout << "\tR2: " << averageErrorRate(read_quality_frequencies_2nd, read_quality_correctness_2nd) << "\n";
			std::cout << std::flush;
		}
	}

	read_length = readLength;
	paranoid = true;
}

std::pair<double, double> readSimulator::averageErrorRate_R1_R2()
{
	return make_pair(averageErrorRate(read_quality_frequencies, read_quality_correctness), averageErrorRate(read_quality_frequencies_2nd, read_quality_correctness_2nd));
}


std::vector<oneReadPair> readSimulator::simulate_paired_reads_from_string(std::string S, double expected_haploid_coverage, double starting_coordinates_diff_mean, double starting_coordinates_diff_sd, bool perfectly, std::string readIDPrefix)
{
	std::vector<Edge*> edgePath;
	for(unsigned int eI = 0; eI < S.size(); eI++)
	{
		std::string edgeEmission = S.substr(eI, 1);
		Edge* e = new Edge();
		e->emission = edgeEmission.at(0);
		edgePath.push_back(e);
	}	
	std::vector<oneReadPair> forReturn = simulate_paired_reads_from_edgePath(edgePath, expected_haploid_coverage, starting_coordinates_diff_mean, starting_coordinates_diff_sd, perfectly, readIDPrefix, false);
	
	for(auto edge : edgePath)
	{
		delete edge;
	}
	
	return forReturn;
	forReturn.clear();

	std::string edgePath_string;
	edgePath_string.reserve(S.size());
	std::vector<unsigned int> edgePath_string_originLevel;

	assert(edgePath_string.size() == edgePath_string_originLevel.size());


	double poissonStartRate = expected_haploid_coverage / ( 2.0 * (double)(read_length)); // this is of reads and their pairs, thus / 2

	long long firstPosition = 0;
	long long lastPosition = edgePath_string.length() - read_length;

	if(!(lastPosition >= firstPosition))
	{
		return forReturn;
		throw std::runtime_error("readSimulator::simulate_paired_reads_from_string(): Problem -- lastPosition < firstPosition -- the supplied first string is not long enough!\n");
	}

	double global_indel_events = 0;
	double global_generated_bases = 0;
	double global_generated_errors = 0;
	std::map<char, double> global_error_NUC;

	{
		boost::mt19937 rnd_gen;

//		todo activate
//		auto seed = boost::random::random_device()();
//		rnd_gen.seed(seed);

		boost::random::poisson_distribution<> rnd_starting_reads ( poissonStartRate );
		std::vector< boost::random::poisson_distribution<> > rnd_INs;
		std::vector< boost::random::poisson_distribution<> > rnd_DELs;
		boost::random::normal_distribution<> rnd_jumpSize (starting_coordinates_diff_mean, starting_coordinates_diff_sd);

		for(unsigned int i = 0; i < read_length; i++)
		{
			rnd_INs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
			rnd_DELs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
		}

		double thread_indel_events = 0;
		double thread_generated_bases = 0;
		double thread_generated_errors = 0;
		size_t thread_read_pairs = 0;

		auto sampleOneBase = [&](unsigned int position_in_read, char underlyingBase, char& returnedBase, char& returnedQuality) -> void {

			assert(position_in_read < this->read_quality_frequencies.size());
			assert(position_in_read < this->read_quality_correctness.size());

			returnedQuality = Utilities::choose_from_normalized_map(this->read_quality_frequencies.at(position_in_read), rnd_gen);

			assert(returnedQuality > 0);

			bool generateError = Utilities::oneBernoulliTrial( 1 - this->read_quality_correctness.at(position_in_read).at(returnedQuality), rnd_gen);

			if(generateError && (! perfectly))
			{
				returnedBase =  Utilities::randomNucleotide(rnd_gen);
				if(global_error_NUC.count(returnedBase) == 0)
				{
					global_error_NUC[returnedBase] = 0;
				}
				global_error_NUC[returnedBase]++;
				thread_generated_errors++;
			}
			else
			{
				returnedBase = underlyingBase;
			}

			thread_generated_bases++;

			returnedQuality += 32;
		};

		auto sampleRead = [&](long long index_into_baseString, std::string& read, std::string& read_qualities, std::vector<int>& coordinates_string, bool& success) -> void {

			read.resize(this->read_length, 0);
			read_qualities.resize(this->read_length, 0);

			coordinates_string.clear();
			success = true;

			int INDEL_events = 0;

			for(unsigned int base = 0; base < this->read_length; base++)
			{
				int insertions = rnd_INs.at(base)(rnd_gen);
				int deletions = rnd_DELs.at(base)(rnd_gen);

				if(perfectly)
				{
					insertions = 0;
					deletions = 0;
				}

				// std::cout << "\tbase " << base << " " << insertions << " " << deletions << "\n" << std::flush;

				INDEL_events += (insertions + deletions);

				if(insertions > 0)
				{
					for(int insertionI = 0; insertionI < insertions; insertionI++)
					{
						char baseChar = Utilities::randomNucleotide(rnd_gen);

						char base_for_read; char quality_for_read;
						sampleOneBase(base, baseChar, base_for_read, quality_for_read);
						read.at(base) = base_for_read;
						read_qualities.at(base) = quality_for_read;
						base++;

						coordinates_string.push_back(-1);

						if(base >= this->read_length)
							break;
					}

					if(base >= this->read_length)
						break;
				}

				if(deletions > 0)
				{
					index_into_baseString += deletions;
				}

				if(index_into_baseString >= edgePath_string.length())
				{
					success = false;
					break;
				}

				char base_for_read; char quality_for_read;
				assert(index_into_baseString < edgePath_string.size());

				sampleOneBase(base, edgePath_string.at(index_into_baseString), base_for_read, quality_for_read);

				read.at(base) = base_for_read;
				read_qualities.at(base) = quality_for_read;
				coordinates_string.push_back(index_into_baseString);

				// std::cout << "index_into_baseString: " << index_into_baseString << "; base: " << base << "; coordinates_string.size(): " << coordinates_string.size() << "\n" << std::flush;
				index_into_baseString++;
			}


			thread_indel_events += INDEL_events;
			// std::cout << "INDEL events: " << INDEL_events << "\n";

			if(!((! success) || (coordinates_string.size() == read.size())))
			{
				std::cerr << "success: " << success << "\n";
				std::cerr << "coordinates_string.size(): " << coordinates_string.size() << "\n";
				std::cerr << "read.size(): " << read.size() << "\n" << std::flush;
			}

			assert((! success) || (coordinates_string.size() == read.size()));

			if(paranoid && success)
			{
				assert(std::find(read.begin(), read.end(), 0) == read.end());
				assert(std::find(read_qualities.begin(), read_qualities.end(), 0) == read_qualities.end());
			}
		};

		long long rPI = 0;
		for(long long i = 0; i < lastPosition; i++)
		{
			int starting_reads = rnd_starting_reads(rnd_gen);

			for(int readI = 0; readI < starting_reads; readI++)
			{
				int jumpSize = floor(rnd_jumpSize(rnd_gen));

				rPI++;

				std::string read1; std::string read1_qualities; std::vector<int> read1_coordinates_string; bool read1_success;
				std::string read2; std::string read2_qualities; std::vector<int> read2_coordinates_string; bool read2_success;

				sampleRead(i,            read1, read1_qualities, read1_coordinates_string, read1_success);
				sampleRead(i + this->read_length + jumpSize, read2, read2_qualities, read2_coordinates_string, read2_success);

				if(read1_success && read2_success)
				{
					thread_read_pairs++;

					read2 = Utilities::seq_reverse_complement(read2);
					std::reverse(read2_qualities.begin(), read2_qualities.end());

					// std::string read1_name = "p1" + readName_field_separator + Utilities::ItoStr(i) + readName_field_separator + Utilities::ItoStr(rPI);
					// std::string read2_name = "p2" + readName_field_separator + Utilities::ItoStr(i + jumpSize) + readName_field_separator + Utilities::ItoStr(rPI);

					std::string read1_name = readIDPrefix + "p" + readName_field_separator + Utilities::ItoStr(i) + readName_field_separator + Utilities::ItoStr(jumpSize) + readName_field_separator + Utilities::ItoStr(rPI);
					std::string read2_name = read1_name;
					
					oneRead r1(read1_name, read1, read1_qualities);
					oneRead r2(read2_name, read2, read2_qualities);

					std::vector<int> read1_coordinates_edgePath;
					std::vector<int> read2_coordinates_edgePath;
//					for(unsigned int cI = 0; cI < read1_coordinates_string.size(); cI++)
//					{
//						int c = read1_coordinates_string.at(cI);
//					}
//					for(unsigned int cI = 0; cI < read2_coordinates_string.size(); cI++)
//					{
//						int c = read2_coordinates_string.at(cI);
//					}

					for(unsigned int cI = 0; cI < read1_coordinates_string.size(); cI++)
					{
						int c = read1_coordinates_string.at(cI);
						if(c != -1)
						{
							read1_coordinates_string.at(cI) = edgePath_string_originLevel.at(c);
						}
					}

					for(unsigned int cI = 0; cI < read2_coordinates_string.size(); cI++)
					{
						int c = read2_coordinates_string.at(cI);
						if(c != -1)
						{
							read2_coordinates_string.at(cI) = edgePath_string_originLevel.at(c);
						}
					}

					r1.coordinates_string = read1_coordinates_string;
					r1.coordinates_edgePath = read1_coordinates_edgePath;

					std::reverse(read2_coordinates_string.begin(), read2_coordinates_string.end());
					std::reverse(read2_coordinates_edgePath.begin(), read2_coordinates_edgePath.end());

					r2.coordinates_string = read2_coordinates_string;
					r2.coordinates_edgePath = read2_coordinates_edgePath;

					assert(r1.coordinates_string.size() == r1.sequence.size());
					assert(r2.coordinates_string.size() == r2.sequence.size());

					oneReadPair rP(r1, r2, jumpSize);

					if(Utilities::oneBernoulliTrial(0.5, rnd_gen))
					{
						rP.invert();
					}

					forReturn.push_back(rP);
				}
			}
		}

		{
			global_generated_bases += thread_generated_bases;
			global_generated_errors += thread_generated_errors;
			global_indel_events += thread_indel_events;
		}
	}

	std::cout << "readSimulator::simulate_paired_reads_from_string(..): Simulated " << forReturn.size() << " read pairs.\n";
	std::cout << "\t" << "global_generated_bases" << ": " << global_generated_bases << "\n";
	std::cout << "\t" << "global_generated_errors" << ": " << global_generated_errors << "\n";
	std::cout << "\t" << "global_indel_events" << ": " << global_indel_events << "\n\n";
	std::cout << "\t" << "error base counts: \n";
	for(std::map<char, double>::iterator bIt = global_error_NUC.begin(); bIt != global_error_NUC.end(); bIt++)
	{
		std::cout << "\t\t" << bIt->first << ": " << bIt->second << "\n";
	}
	std::cout << "\n" << std::flush;

	return forReturn;
}


void readSimulator::simulate_paired_reads_from_string_mt_immediateOutput(std::string S, std::string fn_output_prefix, double expected_haploid_coverage, double starting_coordinates_diff_mean, double starting_coordinates_diff_sd, bool perfectly, unsigned int threads, std::string readIDPrefix)
{
	const std::string edgePath_string = S;

	std::vector<std::ofstream*> stream_FASTQ_r1;
	std::vector<std::ofstream*> stream_FASTQ_r2;

	assert(threads > 0);
	stream_FASTQ_r1.resize(threads, 0);
	stream_FASTQ_r2.resize(threads, 0);
	for(unsigned int tI = 0; tI < threads; tI++)
	{
		stream_FASTQ_r1.at(tI) = new std::ofstream();
		stream_FASTQ_r2.at(tI) = new std::ofstream();

		std::string fn_R1 = fn_output_prefix + "T" + Utilities::ItoStr(tI) + "_1.fq";
		std::string fn_R2 = fn_output_prefix + "T" + Utilities::ItoStr(tI) + "_2.fq";

		stream_FASTQ_r1.at(tI)->open(fn_R1.c_str());
		assert(stream_FASTQ_r1.at(tI)->is_open());

		stream_FASTQ_r2.at(tI)->open(fn_R2.c_str());
		assert(stream_FASTQ_r2.at(tI)->is_open());
	}

	double poissonStartRate = expected_haploid_coverage / ( 2.0 * (double)(read_length)); // this is of reads and their pairs, thus / 2

	long long firstPosition = 0;
	long long lastPosition = edgePath_string.length() - read_length;

	if(!(lastPosition >= firstPosition))
	{
		throw std::runtime_error("readSimulator::simulate_paired_reads_from_string_immediateOutput(): Problem -- lastPosition < firstPosition -- the supplied first string is not long enough!\n");
	}

	double global_indel_events = 0;
	double global_generated_bases = 0;
	double global_generated_errors = 0;
	std::map<char, size_t> global_error_NUC;

	size_t global_reads_simulated = 0;

	long long chunkSize = lastPosition / threads;
	omp_set_num_threads(threads);

	#pragma omp parallel
	{
		unsigned int tI = omp_get_thread_num();

		auto print_one_readPair = [&] (const oneReadPair& rP) -> void
		{
			(*stream_FASTQ_r1.at(tI)) << "@" << rP.reads.first.name << "\n";
			(*stream_FASTQ_r1.at(tI)) << rP.reads.first.sequence << "\n";
			(*stream_FASTQ_r1.at(tI)) << "+" << "\n";
			(*stream_FASTQ_r1.at(tI)) << rP.reads.first.quality << "\n";

			(*stream_FASTQ_r2.at(tI)) << "@" << rP.reads.second.name << "\n";
			(*stream_FASTQ_r2.at(tI)) << rP.reads.second.sequence << "\n";
			(*stream_FASTQ_r2.at(tI)) << "+" << "\n";
			(*stream_FASTQ_r2.at(tI)) << rP.reads.second.quality << "\n";
		};

		long long thisChunk_start = tI * chunkSize;
		long long thisChunk_stop = (tI+1) * chunkSize - 1;
		if(tI == (threads - 1))
		{
			thisChunk_stop = lastPosition;
		}
		assert(thisChunk_start >= 0);
		assert(thisChunk_stop > thisChunk_start);
		assert(thisChunk_stop <= lastPosition);

		boost::mt19937 rnd_gen;

		//		todo activate
		//		auto seed = boost::random::random_device()();
		//		rnd_gen.seed(seed);

		boost::random::poisson_distribution<> rnd_starting_reads ( poissonStartRate );
		std::vector< boost::random::poisson_distribution<> > rnd_INs;
		std::vector< boost::random::poisson_distribution<> > rnd_DELs;
		boost::random::normal_distribution<> rnd_jumpSize (starting_coordinates_diff_mean, starting_coordinates_diff_sd);

		for(unsigned int i = 0; i < read_length; i++)
		{
			rnd_INs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
			rnd_DELs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
		}

		double thread_indel_events = 0;
		double thread_generated_bases = 0;
		double thread_generated_errors = 0;
		size_t thread_read_pairs = 0;
		std::map<char, size_t> thread_error_NUC;

		auto sampleOneBase = [&](unsigned int position_in_read, char underlyingBase, char& returnedBase, char& returnedQuality) -> void {

			assert(position_in_read < this->read_quality_frequencies.size());
			assert(position_in_read < this->read_quality_correctness.size());

			returnedQuality = Utilities::choose_from_normalized_map(this->read_quality_frequencies.at(position_in_read), rnd_gen);

			assert(returnedQuality > 0);

			bool generateError = Utilities::oneBernoulliTrial( 1 - this->read_quality_correctness.at(position_in_read).at(returnedQuality), rnd_gen);

			if(generateError && (! perfectly))
			{
				returnedBase =  Utilities::randomNucleotide(rnd_gen);
				if(thread_error_NUC.count(returnedBase) == 0)
				{
					thread_error_NUC[returnedBase] = 0;
				}
				thread_error_NUC[returnedBase]++;
				thread_generated_errors++;
			}
			else
			{
				returnedBase = underlyingBase;
			}

			thread_generated_bases++;

			returnedQuality += 32;
		};

		auto sampleRead = [&](long long index_into_baseString, std::string& read, std::string& read_qualities, std::vector<int>& coordinates_string, bool& success) -> void {

			read.resize(this->read_length, 0);
			read_qualities.resize(this->read_length, 0);

			coordinates_string.clear();
			success = true;

			int INDEL_events = 0;

			for(unsigned int base = 0; base < this->read_length; base++)
			{
				int insertions = rnd_INs.at(base)(rnd_gen);
				int deletions = rnd_DELs.at(base)(rnd_gen);

				if(perfectly)
				{
					insertions = 0;
					deletions = 0;
				}

				// std::cout << "\tbase " << base << " " << insertions << " " << deletions << "\n" << std::flush;

				INDEL_events += (insertions + deletions);

				if(insertions > 0)
				{
					for(int insertionI = 0; insertionI < insertions; insertionI++)
					{
						char baseChar = Utilities::randomNucleotide(rnd_gen);

						char base_for_read; char quality_for_read;
						sampleOneBase(base, baseChar, base_for_read, quality_for_read);
						read.at(base) = base_for_read;
						read_qualities.at(base) = quality_for_read;
						base++;

						coordinates_string.push_back(-1);

						if(base >= this->read_length)
							break;
					}

					if(base >= this->read_length)
						break;
				}

				if(deletions > 0)
				{
					index_into_baseString += deletions;
				}

				if(index_into_baseString >= edgePath_string.length())
				{
					success = false;
					break;
				}

				char base_for_read; char quality_for_read;
				assert(index_into_baseString < edgePath_string.size());

				sampleOneBase(base, edgePath_string.at(index_into_baseString), base_for_read, quality_for_read);

				read.at(base) = base_for_read;
				read_qualities.at(base) = quality_for_read;
				coordinates_string.push_back(index_into_baseString);

				// std::cout << "index_into_baseString: " << index_into_baseString << "; base: " << base << "; coordinates_string.size(): " << coordinates_string.size() << "\n" << std::flush;
				index_into_baseString++;
			}


			thread_indel_events += INDEL_events;
			// std::cout << "INDEL events: " << INDEL_events << "\n";

			if(!((! success) || (coordinates_string.size() == read.size())))
			{
				std::cerr << "success: " << success << "\n";
				std::cerr << "coordinates_string.size(): " << coordinates_string.size() << "\n";
				std::cerr << "read.size(): " << read.size() << "\n" << std::flush;
			}

			assert((! success) || (coordinates_string.size() == read.size()));

			if(paranoid && success)
			{
				assert(std::find(read.begin(), read.end(), 0) == read.end());
				assert(std::find(read_qualities.begin(), read_qualities.end(), 0) == read_qualities.end());
			}
		};

		long long rPI = 0;
		for(long long i = thisChunk_start; i <= thisChunk_stop; i++)
		{
			int starting_reads = rnd_starting_reads(rnd_gen);

			for(int readI = 0; readI < starting_reads; readI++)
			{
				int jumpSize = floor(rnd_jumpSize(rnd_gen));

				rPI++;

				std::string read1; std::string read1_qualities; std::vector<int> read1_coordinates_string; bool read1_success;
				std::string read2; std::string read2_qualities; std::vector<int> read2_coordinates_string; bool read2_success;

				sampleRead(i,            read1, read1_qualities, read1_coordinates_string, read1_success);
				sampleRead(i + this->read_length + jumpSize, read2, read2_qualities, read2_coordinates_string, read2_success);

				if(read1_success && read2_success)
				{
					thread_read_pairs++;

					read2 = Utilities::seq_reverse_complement(read2);
					std::reverse(read2_qualities.begin(), read2_qualities.end());

					// std::string read1_name = "p1" + readName_field_separator + Utilities::ItoStr(i) + readName_field_separator + Utilities::ItoStr(rPI);
					// std::string read2_name = "p2" + readName_field_separator + Utilities::ItoStr(i + jumpSize) + readName_field_separator + Utilities::ItoStr(rPI);

					std::string read1_name = readIDPrefix + "p" + readName_field_separator + Utilities::ItoStr(i) + readName_field_separator + Utilities::ItoStr(jumpSize) + readName_field_separator + "t" + Utilities::ItoStr(tI) + readName_field_separator + Utilities::ItoStr(rPI);
					std::string read2_name = read1_name;

					oneRead r1(read1_name, read1, read1_qualities);
					oneRead r2(read2_name, read2, read2_qualities);

					std::vector<int> read1_coordinates_edgePath;
					std::vector<int> read2_coordinates_edgePath;
//					for(unsigned int cI = 0; cI < read1_coordinates_string.size(); cI++)
//					{
//						int c = read1_coordinates_string.at(cI);
//					}
//					for(unsigned int cI = 0; cI < read2_coordinates_string.size(); cI++)
//					{
//						int c = read2_coordinates_string.at(cI);
//					}

					r1.coordinates_string = read1_coordinates_string;
					r1.coordinates_edgePath = read1_coordinates_edgePath;

					std::reverse(read2_coordinates_string.begin(), read2_coordinates_string.end());
					std::reverse(read2_coordinates_edgePath.begin(), read2_coordinates_edgePath.end());

					r2.coordinates_string = read2_coordinates_string;
					r2.coordinates_edgePath = read2_coordinates_edgePath;

					assert(r1.coordinates_string.size() == r1.sequence.size());
					assert(r2.coordinates_string.size() == r2.sequence.size());

					oneReadPair rP(r1, r2, jumpSize);

					if(Utilities::oneBernoulliTrial(0.5, rnd_gen))
					{
						rP.invert();
					}

					print_one_readPair(rP);
				}
			}
		}

		#pragma omp critical
		{
			global_reads_simulated += thread_read_pairs;
			global_generated_bases += thread_generated_bases;
			global_generated_errors += thread_generated_errors;
			global_indel_events += thread_indel_events;

			for(std::map<char, size_t>::iterator nucIt = thread_error_NUC.begin(); nucIt != thread_error_NUC.end(); nucIt++)
			{
				if(global_error_NUC.count(nucIt->first) == 0)
				{
					global_error_NUC[nucIt->first] = 0;
				}
				global_error_NUC.at(nucIt->first) +=  nucIt->second;
			}

		}
	}

	std::cout << "readSimulator::simulate_paired_reads_from_string_immediateOutput(..): Simulated " << global_reads_simulated << " read pairs.\n";
	std::cout << "\t" << "global_generated_bases" << ": " << global_generated_bases << "\n";
	std::cout << "\t" << "global_generated_errors" << ": " << global_generated_errors << "\n";
	std::cout << "\t" << "global_indel_events" << ": " << global_indel_events << "\n\n";
	std::cout << "\t" << "error base counts: \n";
	for(std::map<char, size_t>::iterator bIt = global_error_NUC.begin(); bIt != global_error_NUC.end(); bIt++)
	{
		std::cout << "\t\t" << bIt->first << ": " << bIt->second << "\n";
	}
	std::cout << "\n" << std::flush;

	for(unsigned int tI = 0; tI < threads; tI++)
	{
		stream_FASTQ_r1.at(tI)->close();
		stream_FASTQ_r2.at(tI)->close();
		delete(stream_FASTQ_r1.at(tI));
		delete(stream_FASTQ_r2.at(tI));
	}
}



std::vector<oneRead> readSimulator::simulate_unpaired_reads_from_string(std::string S, double expected_haploid_coverage, bool perfectly)
{
	std::vector<oneRead> forReturn;

	std::string edgePath_string = S;

	double poissonStartRate = expected_haploid_coverage / (double)read_length; // this is of reads and their pairs, thus / 2

	long long firstPosition = 0;
	long long lastPosition = edgePath_string.length() - read_length;

	if(!(lastPosition >= firstPosition))
	{
		throw std::runtime_error("readSimulator::simulate_unpaired_reads_from_string(): Problem -- lastPosition < firstPosition -- the supplied first string is not long enough!\n");
	}

	double global_indel_events = 0;
	double global_generated_bases = 0;
	double global_generated_errors = 0;
	std::map<char, double> global_error_NUC;

	{
		boost::mt19937 rnd_gen;

		//		todo activate
		//		auto seed = boost::random::random_device()();
		//		rnd_gen.seed(seed);

		boost::random::poisson_distribution<> rnd_starting_reads ( poissonStartRate );
		std::vector< boost::random::poisson_distribution<> > rnd_INs;
		std::vector< boost::random::poisson_distribution<> > rnd_DELs;

		for(unsigned int i = 0; i < read_length; i++)
		{
			rnd_INs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
			rnd_DELs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
		}

		double thread_indel_events = 0;
		double thread_generated_bases = 0;
		double thread_generated_errors = 0;
		size_t thread_read_pairs = 0;

		auto sampleOneBase = [&](unsigned int position_in_read, char underlyingBase, char& returnedBase, char& returnedQuality) -> void {

			assert(position_in_read < this->read_quality_frequencies.size());
			assert(position_in_read < this->read_quality_correctness.size());

			returnedQuality = Utilities::choose_from_normalized_map(this->read_quality_frequencies.at(position_in_read), rnd_gen);

			assert(returnedQuality > 0);

			bool generateError = Utilities::oneBernoulliTrial( 1 - this->read_quality_correctness.at(position_in_read).at(returnedQuality), rnd_gen);

			if(generateError && (! perfectly))
			{
				returnedBase =  Utilities::randomNucleotide(rnd_gen);
				if(global_error_NUC.count(returnedBase) == 0)
				{
					global_error_NUC[returnedBase] = 0;
				}
				global_error_NUC[returnedBase]++;
				thread_generated_errors++;
			}
			else
			{
				returnedBase = underlyingBase;
			}

			thread_generated_bases++;

			returnedQuality += 32;
		};

		auto sampleRead = [&](long long index_into_baseString, std::string& read, std::string& read_qualities, std::vector<int>& coordinates_string, bool& success) -> void {

			read.resize(this->read_length, 0);
			read_qualities.resize(this->read_length, 0);

			coordinates_string.clear();
			success = true;

			int INDEL_events = 0;

			for(unsigned int base = 0; base < this->read_length; base++)
			{
				int insertions = rnd_INs.at(base)(rnd_gen);
				int deletions = rnd_DELs.at(base)(rnd_gen);

				if(perfectly)
				{
					insertions = 0;
					deletions = 0;
				}

				// std::cout << "\tbase " << base << " " << insertions << " " << deletions << "\n" << std::flush;

				INDEL_events += (insertions + deletions);

				if(insertions > 0)
				{
					for(int insertionI = 0; insertionI < insertions; insertionI++)
					{
						char baseChar = Utilities::randomNucleotide(rnd_gen);

						char base_for_read; char quality_for_read;
						sampleOneBase(base, baseChar, base_for_read, quality_for_read);
						read.at(base) = base_for_read;
						read_qualities.at(base) = quality_for_read;
						base++;

						coordinates_string.push_back(-1);

						if(base >= this->read_length)
							break;
					}

					if(base >= this->read_length)
						break;
				}

				if(deletions > 0)
				{
					index_into_baseString += deletions;
				}

				if(index_into_baseString >= edgePath_string.length())
				{
					success = false;
					break;
				}

				char base_for_read; char quality_for_read;
				assert(index_into_baseString < edgePath_string.size());

				sampleOneBase(base, edgePath_string.at(index_into_baseString), base_for_read, quality_for_read);

				read.at(base) = base_for_read;
				read_qualities.at(base) = quality_for_read;
				coordinates_string.push_back(index_into_baseString);

				// std::cout << "index_into_baseString: " << index_into_baseString << "; base: " << base << "; coordinates_string.size(): " << coordinates_string.size() << "\n" << std::flush;
				index_into_baseString++;
			}


			thread_indel_events += INDEL_events;
			// std::cout << "INDEL events: " << INDEL_events << "\n";

			if(!((! success) || (coordinates_string.size() == read.size())))
			{
				std::cerr << "success: " << success << "\n";
				std::cerr << "coordinates_string.size(): " << coordinates_string.size() << "\n";
				std::cerr << "read.size(): " << read.size() << "\n" << std::flush;
			}

			assert((! success) || (coordinates_string.size() == read.size()));

			if(paranoid && success)
			{
				assert(std::find(read.begin(), read.end(), 0) == read.end());
				assert(std::find(read_qualities.begin(), read_qualities.end(), 0) == read_qualities.end());
			}
		};

		long long rPI = 0;
		for(long long i = 0; i < lastPosition; i++)
		{
			int starting_reads = rnd_starting_reads(rnd_gen);

			for(int readI = 0; readI < starting_reads; readI++)
			{
				rPI++;

				std::string read1; std::string read1_qualities; std::vector<int> read1_coordinates_string; bool read1_success;

				sampleRead(i,            read1, read1_qualities, read1_coordinates_string, read1_success);

				if(read1_success)
				{
					thread_read_pairs++;

					// std::string read1_name = "p1" + readName_field_separator + Utilities::ItoStr(i) + readName_field_separator + Utilities::ItoStr(rPI);
					// std::string read2_name = "p2" + readName_field_separator + Utilities::ItoStr(i + jumpSize) + readName_field_separator + Utilities::ItoStr(rPI);

					std::string read1_name = "p" + readName_field_separator + Utilities::ItoStr(i) + readName_field_separator + Utilities::ItoStr(rPI);

					oneRead r1(read1_name, read1, read1_qualities);

					std::vector<int> read1_coordinates_edgePath;
//					for(unsigned int cI = 0; cI < read1_coordinates_string.size(); cI++)
//					{
//						int c = read1_coordinates_string.at(cI);
//					}
//					for(unsigned int cI = 0; cI < read2_coordinates_string.size(); cI++)
//					{
//						int c = read2_coordinates_string.at(cI);
//					}

					r1.coordinates_string = read1_coordinates_string;
					r1.coordinates_edgePath = read1_coordinates_edgePath;

					assert(r1.coordinates_string.size() == r1.sequence.size());

					forReturn.push_back(r1);
				}
			}
		}

		{
			global_generated_bases += thread_generated_bases;
			global_generated_errors += thread_generated_errors;
			global_indel_events += thread_indel_events;
		}
	}

	std::cout << "readSimulator::simulate_paired_reads_from_string(..): Simulated " << forReturn.size() << " read pairs.\n";
	std::cout << "\t" << "global_generated_bases" << ": " << global_generated_bases << "\n";
	std::cout << "\t" << "global_generated_errors" << ": " << global_generated_errors << "\n";
	std::cout << "\t" << "global_indel_events" << ": " << global_indel_events << "\n\n";
	std::cout << "\t" << "error base counts: \n";
	for(std::map<char, double>::iterator bIt = global_error_NUC.begin(); bIt != global_error_NUC.end(); bIt++)
	{
		std::cout << "\t\t" << bIt->first << ": " << bIt->second << "\n";
	}
	std::cout << "\n" << std::flush;

	return forReturn;
}

std::vector<oneReadPair> readSimulator::simulate_paired_reads_from_edgePath(std::vector<Edge*> edgePath, double expected_haploid_coverage, double starting_coordinates_diff_mean, double starting_coordinates_diff_sd, bool perfectly, std::string readIDPrefix, bool includeDeletions)
{
	std::vector<oneReadPair> forReturn;

	std::string edgePath_string;
	std::vector<unsigned int> edgePath_string_originLevel;
	for(unsigned int eI = 0; eI < edgePath.size(); eI++)
	{
		std::string edgeEmission = edgePath.at(eI)->getEmission();
		/*
		if(edgeEmission == "*")
		{
			std::cerr << "Edge emission " << edgeEmission << " at level " << eI << " / " << edgePath.size() << "\n" << std::flush;
		}
		assert(edgeEmission != "*");
		*/ 
		assert(edgeEmission.length() == 1);
		if(edgeEmission != "_")
		{
			edgePath_string.append(edgeEmission);
			edgePath_string_originLevel.push_back(eI);
		}
	}
	assert(edgePath_string.size() == edgePath_string_originLevel.size());

	double poissonStartRate = expected_haploid_coverage / ( 2.0 * (double)(read_length)); // this is of reads and their pairs, thus / 2

	long long firstPosition = 0;
	long long lastPosition = edgePath_string.length() - read_length;

	if(!(lastPosition >= firstPosition))
	{
		// return forReturn;
		throw std::runtime_error("readSimulator::simulate_paired_reads_from_edgePath(): Problem -- lastPosition < firstPosition -- the supplied first string is not long enough!\n");
	}

	double global_indel_events = 0;
	double global_generated_bases = 0;
	double global_generated_errors = 0;
	std::map<char, double> global_error_NUC;
	
	{
		boost::mt19937 rnd_gen;

		//		todo activate
		//		auto seed = boost::random::random_device()();
		//		rnd_gen.seed(seed);

		
		boost::random::poisson_distribution<> rnd_starting_reads ( poissonStartRate );
		std::vector< boost::random::poisson_distribution<> > rnd_INs;
		std::vector< boost::random::poisson_distribution<> > rnd_DELs;
		boost::random::normal_distribution<> rnd_jumpSize (starting_coordinates_diff_mean, starting_coordinates_diff_sd);

		for(unsigned int i = 0; i < read_length; i++)
		{
			rnd_INs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
			rnd_DELs.push_back(boost::random::poisson_distribution<>( read_INDEL_freq.at(i) ));
		}

		double thread_indel_events = 0;
		double thread_generated_bases = 0;
		double thread_generated_errors = 0;
		size_t thread_read_pairs = 0;

		auto sampleOneBase = [&](unsigned int position_in_read, char underlyingBase, char& returnedBase, char& returnedQuality, bool is2nd) -> void {

			assert(position_in_read < this->read_quality_frequencies.size());
			assert(position_in_read < this->read_quality_correctness.size());

			if(is2nd)
			{
				returnedQuality = Utilities::choose_from_normalized_map(this->read_quality_frequencies_2nd.at(position_in_read), rnd_gen);
			}
			else
			{
				returnedQuality = Utilities::choose_from_normalized_map(this->read_quality_frequencies.at(position_in_read), rnd_gen);
			}


			assert(returnedQuality > 0);

			bool generateError;

			if(is2nd)
			{
				generateError = Utilities::oneBernoulliTrial( 1 - this->read_quality_correctness_2nd.at(position_in_read).at(returnedQuality), rnd_gen);
			}
			else
			{
				generateError = Utilities::oneBernoulliTrial( 1 - this->read_quality_correctness.at(position_in_read).at(returnedQuality), rnd_gen);
			}
						
			if(generateError && (! perfectly))
			{
				returnedBase =  Utilities::randomNucleotide(rnd_gen);
				if(global_error_NUC.count(returnedBase) == 0)
				{
					global_error_NUC[returnedBase] = 0;
				}
				global_error_NUC[returnedBase]++;
				thread_generated_errors++;
			}
			else
			{
				returnedBase = underlyingBase;
			}

			thread_generated_bases++;

			returnedQuality += 32;
		};
		
		auto sampleRead = [&](long long index_into_baseString, std::string& read, std::string& read_qualities, std::vector<int>& coordinates_string, std::string& read_underlyingEdgeLabels, std::vector<int>& fA_coordinates_edgePath, std::string& fA_underlyingEdgeLabels, std::string& fA_sequence, bool& success, bool is2nd) -> void {

			read.resize(this->read_length, 0);
			read_qualities.resize(this->read_length, 0);
			read_underlyingEdgeLabels.resize(this->read_length, 0);

			coordinates_string.clear();
			fA_coordinates_edgePath.clear();
			fA_underlyingEdgeLabels.clear();
			fA_sequence.clear();

			coordinates_string.reserve(this->read_length);
			fA_coordinates_edgePath.reserve(this->read_length);
			fA_underlyingEdgeLabels.reserve(this->read_length);
			fA_sequence.reserve(this->read_length);

			success = true;

			int INDEL_events = 0;

			int outstandingDeletions = 0;
			for(unsigned int base = 0; base < this->read_length; base++)
			{
				assert(outstandingDeletions >= 0);

				int insertions = rnd_INs.at(base)(rnd_gen);
				int deletions = rnd_DELs.at(base)(rnd_gen);

				if(perfectly)
				{
					insertions = 0;
					deletions = 0;
				}
				
				// std::cout << "\tbase " << base << " " << insertions << " " << deletions << "\n" << std::flush;

				INDEL_events += (insertions + deletions);

				if(insertions > 0)
				{
					for(int insertionI = 0; insertionI < insertions; insertionI++)
					{
						char baseChar = Utilities::randomNucleotide(rnd_gen);

						char base_for_read; char quality_for_read;
						sampleOneBase(base, baseChar, base_for_read, quality_for_read, is2nd);
						read.at(base) = base_for_read;
						read_qualities.at(base) = quality_for_read;
						read_underlyingEdgeLabels.at(base) = '_';
						base++;

						coordinates_string.push_back(-1);
						
						fA_coordinates_edgePath.push_back(-1);
						fA_underlyingEdgeLabels.push_back('_');
						fA_sequence.push_back(base_for_read);

						if(base >= this->read_length)
							break;
					}

					if(base >= this->read_length)
						break;
				}

				if(includeDeletions)
				{
					outstandingDeletions += deletions;
				}
				else
				{
					if(deletions > 0)
					{
						for(int dI = 0; dI < deletions; dI++)
						{
							fA_coordinates_edgePath.push_back(index_into_baseString);
							fA_underlyingEdgeLabels.push_back(edgePath_string.at(index_into_baseString));
							fA_sequence.push_back('_');
							index_into_baseString++;
						}
						//index_into_baseString += deletions;
					}
				}

				if(index_into_baseString >= edgePath_string.length())
				{
					success = false;
					break;
				}

				char base_for_read; char quality_for_read;
				assert(index_into_baseString < edgePath_string.size());

				if(includeDeletions && outstandingDeletions)
				{
					read.at(base) = '_';
					read_qualities.at(base) = 1;
					read_underlyingEdgeLabels.at(base) = edgePath_string.at(index_into_baseString);
					coordinates_string.push_back(index_into_baseString);

					fA_coordinates_edgePath.push_back(index_into_baseString);
					fA_underlyingEdgeLabels.push_back(edgePath_string.at(index_into_baseString));
					fA_sequence.push_back('_');

					// std::cout << "index_into_baseString: " << index_into_baseString << "; base: " << base << "; coordinates_string.size(): " << coordinates_string.size() << "\n" << std::flush;
					index_into_baseString++;
					outstandingDeletions--;
				}
				else
				{
					assert(outstandingDeletions == 0);
					sampleOneBase(base, edgePath_string.at(index_into_baseString), base_for_read, quality_for_read, is2nd);

					read.at(base) = base_for_read;
					read_qualities.at(base) = quality_for_read;
					read_underlyingEdgeLabels.at(base) = edgePath_string.at(index_into_baseString);
					coordinates_string.push_back(index_into_baseString);

					fA_coordinates_edgePath.push_back(index_into_baseString);
					fA_underlyingEdgeLabels.push_back(edgePath_string.at(index_into_baseString));
					fA_sequence.push_back(base_for_read);

					// std::cout << "index_into_baseString: " << index_into_baseString << "; base: " << base << "; coordinates_string.size(): " << coordinates_string.size() << "\n" << std::flush;
					index_into_baseString++;
				}
			}
		

			thread_indel_events += INDEL_events;
			// std::cout << "INDEL events: " << INDEL_events << "\n";

			if(!((! success) || (coordinates_string.size() == read.size())))
			{
				std::cerr << "success: " << success << "\n";
				std::cerr << "coordinates_string.size(): " << coordinates_string.size() << "\n";
				std::cerr << "read.size(): " << read.size() << "\n" << std::flush;
			}
			
			assert((! success) || (coordinates_string.size() == read.size()));
			
			if(paranoid && success)
			{
				assert(std::find(read.begin(), read.end(), 0) == read.end());
				if(! includeDeletions)
				assert(std::find(read_qualities.begin(), read_qualities.end(), 0) == read_qualities.end());
			}

			unsigned int remove_Gaps_beginningSequence = 0;
			for(unsigned int i = 0; i  < fA_sequence.size(); i++)
			{
				if(fA_sequence.at(i) == '_')
				{
					remove_Gaps_beginningSequence++;
				}
				else
				{
					break;
				}
			}

			unsigned int remove_Gaps_endOfSequence = 0;
			for(int i = (fA_sequence.size() - 1); i  >= 0; i++)
			{
				if(fA_sequence.at(i) == '_')
				{
					remove_Gaps_endOfSequence++;
				}
				else
				{
					break;
				}
			}

			if(remove_Gaps_beginningSequence || remove_Gaps_beginningSequence)
			{
//				std::cerr << "Before modification:\n";
//				std::cerr << "\t" << Utilities::join(Utilities::ItoStr(fA_coordinates_edgePath), " ") << "\n";
//				std::cerr << "\t" << fA_underlyingEdgeLabels << "\n";
//				std::cerr << "\t" << fA_sequence << "\n" << std::flush;
			}

			if(remove_Gaps_beginningSequence)
			{
				assert(remove_Gaps_beginningSequence < fA_coordinates_edgePath.size());
				fA_coordinates_edgePath = std::vector<int>(fA_coordinates_edgePath.begin()+remove_Gaps_beginningSequence, fA_coordinates_edgePath.end());
				fA_sequence = fA_sequence.substr(remove_Gaps_beginningSequence);
				fA_underlyingEdgeLabels = fA_underlyingEdgeLabels.substr(remove_Gaps_beginningSequence);
			}

			if(remove_Gaps_endOfSequence)
			{
				assert(remove_Gaps_endOfSequence < fA_coordinates_edgePath.size());
				fA_coordinates_edgePath = std::vector<int>(fA_coordinates_edgePath.begin(), fA_coordinates_edgePath.begin() + (fA_coordinates_edgePath.size() - remove_Gaps_endOfSequence));
				fA_sequence = fA_sequence.substr(0, fA_sequence.length() - remove_Gaps_endOfSequence);
				fA_underlyingEdgeLabels = fA_underlyingEdgeLabels.substr(0, fA_underlyingEdgeLabels.length() - remove_Gaps_endOfSequence);
			}

			if(remove_Gaps_beginningSequence || remove_Gaps_beginningSequence)
			{
//				std::cerr << "AFTER modification:\n";
//				std::cerr << "\t" << Utilities::join(Utilities::ItoStr(fA_coordinates_edgePath), " ") << "\n";
//				std::cerr << "\t" << fA_underlyingEdgeLabels << "\n";
//				std::cerr << "\t" << fA_sequence << "\n" << std::flush;
			}

			assert(fA_coordinates_edgePath.size() == fA_sequence.size());
			assert(fA_coordinates_edgePath.size() == fA_underlyingEdgeLabels.size());
		};

		for(long long i = 0; i < lastPosition; i++)
		{
			int starting_reads = rnd_starting_reads(rnd_gen);

			for(int readI = 0; readI < starting_reads; readI++)
			{
				int jumpSize = floor(rnd_jumpSize(rnd_gen));

				bool willInvert = false;
				bool firstReadSampled_second = false;
				bool secondReadSampled_second = true;

				if(Utilities::oneBernoulliTrial(0.5, rnd_gen))
				{
					willInvert = true;
					firstReadSampled_second = true;
					secondReadSampled_second = false;
				}
				assert(!(firstReadSampled_second && secondReadSampled_second));

				std::string read1; std::string read1_qualities; std::vector<int> read1_coordinates_string; std::string read1_underlyingEdgeLabels; bool read1_success;
				std::string read2; std::string read2_qualities; std::vector<int> read2_coordinates_string; std::string read2_underlyingEdgeLabels; bool read2_success;

				std::vector<int> read1_fA_coordinates_edgePath; std::string r1_fA_underlyingEdgeLabels; std::string r1_fA_sequence;
				std::vector<int> read2_fA_coordinates_edgePath; std::string r2_fA_underlyingEdgeLabels; std::string r2_fA_sequence;

				sampleRead(i,            read1, read1_qualities, read1_coordinates_string, read1_underlyingEdgeLabels, read1_fA_coordinates_edgePath, r1_fA_underlyingEdgeLabels, r1_fA_sequence, read1_success, firstReadSampled_second);
				sampleRead(i + this->read_length + jumpSize, read2, read2_qualities, read2_coordinates_string, read2_underlyingEdgeLabels, read2_fA_coordinates_edgePath, r2_fA_underlyingEdgeLabels, r2_fA_sequence, read2_success, secondReadSampled_second);

				if(read1_success && read2_success)
				{
					thread_read_pairs++;

					read2 = Utilities::seq_reverse_complement(read2);
					std::reverse(read2_qualities.begin(), read2_qualities.end());
					
					std::string read1_name = readIDPrefix + "r" + Utilities::ItoStr(thread_read_pairs)  + readName_field_separator + Utilities::ItoStr(i) + readName_field_separator + Utilities::ItoStr(jumpSize);
					// std::string read2_name = readIDPrefix + "p2" + readName_field_separator + Utilities::ItoStr(i + jumpSize);
					std::string read2_name = read1_name; 

					oneRead r1(read1_name, read1, read1_qualities);
					oneRead r2(read2_name, read2, read2_qualities);

					std::vector<int> read1_coordinates_edgePath;
					std::vector<int> read2_coordinates_edgePath;
					for(unsigned int cI = 0; cI < read1_coordinates_string.size(); cI++)
					{
						int c = read1_coordinates_string.at(cI);
						if(c == -1)
						{
							read1_coordinates_edgePath.push_back(c);
						}
						else
						{
							read1_coordinates_edgePath.push_back(edgePath_string_originLevel.at(c));
						}
					}
					for(unsigned int cI = 0; cI < read2_coordinates_string.size(); cI++)
					{
						int c = read2_coordinates_string.at(cI);
						if(c == -1)
						{
							read2_coordinates_edgePath.push_back(c);
						}
						else
						{
							read2_coordinates_edgePath.push_back(edgePath_string_originLevel.at(c));
						}
					}

					for(unsigned int cI = 0; cI < read1_fA_coordinates_edgePath.size(); cI++)
					{
						int c = read1_fA_coordinates_edgePath.at(cI);
						if(c == -1)
						{
							read1_fA_coordinates_edgePath.at(cI) = c;
						}
						else
						{
							read1_fA_coordinates_edgePath.at(cI) = edgePath_string_originLevel.at(c);
						}
					}
					for(unsigned int cI = 0; cI < read2_fA_coordinates_edgePath.size(); cI++)
					{
						int c = read2_fA_coordinates_edgePath.at(cI);
						if(c == -1)
						{
							read2_fA_coordinates_edgePath.at(cI) = c;
						}
						else
						{
							read2_fA_coordinates_edgePath.at(cI) = edgePath_string_originLevel.at(c);
						}
					}
					r1.coordinates_string = read1_coordinates_string;
					r1.coordinates_edgePath = read1_coordinates_edgePath;
					r1.underlyingEdgeLabels = read1_underlyingEdgeLabels;
					r1.fullAlignment_coordinates_edgePath = read1_fA_coordinates_edgePath;
					r1.fullAlignment_underlyingEdgeLabels = r1_fA_underlyingEdgeLabels;
					r1.fullAlignment_sequence = r1_fA_sequence;

					
					std::reverse(read2_coordinates_string.begin(), read2_coordinates_string.end());
					std::reverse(read2_coordinates_edgePath.begin(), read2_coordinates_edgePath.end());
					std::reverse(read2_underlyingEdgeLabels.begin(), read2_underlyingEdgeLabels.end());
					std::reverse(read2_fA_coordinates_edgePath.begin(), read2_fA_coordinates_edgePath.end());
					std::reverse(r2_fA_underlyingEdgeLabels.begin(), r2_fA_underlyingEdgeLabels.end());
					std::reverse(r2_fA_sequence.begin(), r2_fA_sequence.end());

					r2.coordinates_string = read2_coordinates_string;
					r2.coordinates_edgePath = read2_coordinates_edgePath;
					r2.underlyingEdgeLabels = read2_underlyingEdgeLabels;
					r2.fullAlignment_coordinates_edgePath = read2_fA_coordinates_edgePath;
					r2.fullAlignment_underlyingEdgeLabels = r2_fA_underlyingEdgeLabels;
					r2.fullAlignment_sequence = r2_fA_sequence;


					assert(r1.coordinates_string.size() == r1.sequence.size());
					assert(r2.coordinates_string.size() == r2.sequence.size());

					assert(r1.coordinates_string.size() == r1.coordinates_edgePath.size());
					assert(r2.coordinates_string.size() == r2.coordinates_edgePath.size());

					if(!(r1.coordinates_string.size() == r1.underlyingEdgeLabels.size()))
					{
						std::cerr << "r1.coordinates_string.size()" << ": " << r1.coordinates_string.size() << "\n";
						std::cerr << "r1.underlyingEdgeLabels.size()" << ": " << r1.underlyingEdgeLabels.size() << "\n";
						std::cerr << "read1_underlyingEdgeLabels.size()" << ": " << read1_underlyingEdgeLabels.size() << "\n";
						std::cerr << "\n" << std::flush;
					}
					assert(r1.coordinates_string.size() == r1.underlyingEdgeLabels.size());
					assert(r2.coordinates_string.size() == r2.underlyingEdgeLabels.size());

					assert(r1.fullAlignment_coordinates_edgePath.size() == r1.fullAlignment_underlyingEdgeLabels.size());
					assert(r1.fullAlignment_coordinates_edgePath.size() == r1.fullAlignment_sequence.size());
					assert(r1.fullAlignment_coordinates_edgePath.size() >= r1.sequence.size());
					assert(r2.fullAlignment_coordinates_edgePath.size() == r2.fullAlignment_underlyingEdgeLabels.size());
					assert(r2.fullAlignment_coordinates_edgePath.size() == r2.fullAlignment_sequence.size());
					assert(r2.fullAlignment_coordinates_edgePath.size() >= r2.sequence.size());



					oneReadPair rP(r1, r2, jumpSize);

					if(willInvert)
					{
						rP.invert();
					}

					rP.firstRead_minusStrand = willInvert;

					forReturn.push_back(rP);
				}
			}
		}
		  
		{
			global_generated_bases += thread_generated_bases;
			global_generated_errors += thread_generated_errors;
			global_indel_events += thread_indel_events;
		}
	}

	if(! quiet)
	{
		std::cout << "readSimulator::simulate_paired_reads_from_edgePath(..): Simulated " << forReturn.size() << " read pairs.\n";
		std::cout << "\t" << "global_generated_bases" << ": " << global_generated_bases << "\n";
		std::cout << "\t" << "global_generated_errors" << ": " << global_generated_errors << "\n";
		std::cout << "\t" << "global_indel_events" << ": " << global_indel_events << "\n\n";
		std::cout << "\t" << "error base counts: \n";
		for(std::map<char, double>::iterator bIt = global_error_NUC.begin(); bIt != global_error_NUC.end(); bIt++)
		{
			std::cout << "\t\t" << bIt->first << ": " << bIt->second << "\n";
		}
		std::cout << "\n" << std::flush;
	}

	return forReturn;
}

}
