/*
 * Utilities.cpp
 *
 *  Created on: 23.05.2011
 *      Author: Alexander Dilthey
 */

#include "Utilities.h"

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <exception>
#include <stdexcept>
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/poisson_distribution.hpp>

#include "dirent.h"

// todo at some point activate...
// unsigned int globalRandRSeed = time(NULL);
unsigned int globalRandRSeed = 0;

using namespace boost::filesystem;


Utilities::Utilities() {
}

Utilities::~Utilities() {
}

int Utilities::randomPoisson(double mean)
{
	boost::mt19937 rnd_gen;
	boost::random::poisson_distribution<> rnd_starting_reads ( mean );
	return rnd_starting_reads(rnd_gen);
}


bool Utilities::directoryExists(std::string dir)
{
	path p(dir);
	return (exists(p) && is_directory(p));
}

void Utilities::makeDir(std::string dir)
{
	if(! directoryExists(dir))
	{
		path p(dir);
		bool success = create_directory(p);
		if(! success)
		{
			throw std::runtime_error("Cannot create directory "+dir);
		}
	}
}

std::string Utilities::getCWD()
{
	boost::filesystem::path full_path( boost::filesystem::current_path() );
	return full_path.string();
}

void Utilities::setCWD(std::string path)
{
	boost::filesystem::current_path(path);
}

std::string Utilities::findFileFromAlternatives(std::vector<std::string> alternatives)
{	
	for(unsigned int i = 0; i < alternatives.size(); i++)
	{
		if(fileExists(alternatives.at(i)))
		{
			return alternatives.at(i);
		}
	}
	
	throw std::runtime_error("Utilities::findFileFromAlternatives(): No alternative present from " +Utilities::join(alternatives, ", "));
}

void Utilities::make_or_clearDirectory(std::string dir)
{
	if(directoryExists(dir))
	{
		clearDirectory(dir);
	}
	else
	{
		makeDir(dir);
	}
}

void Utilities::deleteFile(std::string file)
{
	assert(fileExists(file));
    boost::filesystem::wpath F(file);
    assert(boost::filesystem::remove(F));
}

size_t Utilities::fileLastWrite(std::string file)
{
	assert(fileExists(file));
    boost::filesystem::wpath F(file);
    return last_write_time(F);
}

void Utilities::clearDirectory(std::string dir)
{
	if(! directoryExists(dir))
	{
		makeDir(dir);
	}

	for(boost::filesystem::directory_iterator itr(dir); itr!=directory_iterator(); ++itr)
	{
	    if(is_regular_file(itr->status()))
	    {
	    	boost::filesystem::remove_all(itr->path());
	    }
	}

}

bool Utilities::fileExists(std::string filepath)
{
	std::ifstream ifile(filepath);
	if(ifile.good())
	{
		ifile.close();
		return true;
	}
	else
	{
		return false;
	}
}

bool Utilities::fileReadable(std::string file)
{
	std::ifstream fS;
	fS.open(file.c_str());
	bool forReturn = fS.is_open();
	fS.close();
	return forReturn;
}

std::vector<int> Utilities::StrtoI(std::vector<std::string> s)
{
	std::vector<int> forReturn;
	for(unsigned int i = 0; i < s.size(); i++)
	{
		forReturn.push_back(StrtoI(s.at(i)));
	}
	return forReturn;
}

bool Utilities::intervalsOverlap(int x1, int x2, int y1, int y2)
{
	assert(x1 <= x2);
	assert(y1 <= y2);
	return (x1 >= y1 && x1 <= y2) ||
		 (x2 >= y1 && x2 <= y2) ||
		 (y1 >= x1 && y1 <= x2) ||
		 (y2 >= x1 && y2 <= x2);
}

unsigned char Utilities::PCorrectToPhred(double PCorrect)
{
	assert(PCorrect >= 0);
	assert(PCorrect <= 1);

	double pWrong = 1 - PCorrect;
	if(pWrong == 0)
	{
		pWrong = 1e-100;
	}
	
	double phred1 = -10.0 * log10(pWrong);
	
	if((phred1 + 33) > 255)
	{
		phred1 = 255 - 33;
	}
	
	assert((phred1 + 33) >= 0);
	
	int forReturn = round(phred1 + 33);
	assert(forReturn >= 0);
	assert(forReturn <= 255);
	
	return (unsigned char)(forReturn);
}

std::string Utilities::JoinMapUInt2Str(std::map<std::string, unsigned int> M)
{
	std::vector<std::string> forReturn_parts;

	for(std::map<std::string, unsigned int>::iterator mIt = M.begin(); mIt != M.end(); mIt++)
	{
		std::string thisE = mIt->first + ":" + ItoStr(mIt->second);
		forReturn_parts.push_back(thisE);
	}

	return join(forReturn_parts, ",");
}


double Utilities::LogSumLogPs (const std::vector<double>& v)
{
	std::pair<double, unsigned int> maxPair = findVectorMax(v);
	double max = maxPair.first;
	double secondTerm = 0;
	for(size_t i = 0; i < v.size(); i++)
	{
		if(i == maxPair.second)
		{
			continue;
		}
		else
		{
			secondTerm += exp(v.at(i) - max);
		}
	}

	double forReturn = (max + log(1 + secondTerm));
	return forReturn;
}

std::pair<double, std::string> Utilities::findStringMapMax(const std::map<std::string, double>& m)
{
	assert(m.size() > 0);
	double max;
	std::string iMax;
	for(std::map<std::string, double>::const_iterator mIt = m.begin(); mIt != m.end(); mIt++)
	{
		if((mIt == m.begin()) || (mIt->second > max))
		{
			max = mIt->second;
			iMax = mIt->first;
		}
	}

	return std::pair<double, std::string>(max, iMax);
}

std::pair<double, int> Utilities::findIntMapMax(std::map<int, double>& m)
{
	assert(m.size() > 0);
	double max;
	unsigned iMax;
	for(std::map<int, double>::iterator mIt = m.begin(); mIt != m.end(); mIt++)
	{
		if((mIt == m.begin()) || (mIt->second > max))
		{
			max = mIt->second;
			iMax = mIt->first;
		}
	}

	return std::pair<double, int>(max, iMax);
}


std::pair<double, int> Utilities::findIntMapMaxP_nonCritical(std::map<int, double>& m, unsigned int* thisSeed)
{
	assert(m.size() > 0);
	double max;
	std::vector<int> iMax;
	for(std::map<int, double>::iterator mIt = m.begin(); mIt != m.end(); mIt++)
	{
		if((mIt == m.begin()) || (mIt->second > max))
		{
			max = mIt->second;
			iMax.clear();
			iMax.push_back(mIt->first);
		}
	}

	if(iMax.size() > 1)
	{
		int selectedI = randomNumber_nonCritical(iMax.size() - 1, thisSeed);
		assert((selectedI >= 0) && (selectedI < (int)iMax.size()));
		return std::pair<double, int>(max, iMax.at(selectedI));
	}
	else
	{
		return std::pair<double, int>(max, iMax.at(0));
	}
}

bool Utilities::extractBit(unsigned int number, unsigned int bit)
{
	unsigned int bitmask = pow((unsigned int)2, bit);
	bool bitValue =  ((number & bitmask) == bitmask);
	return bitValue;
}

std::pair<double, unsigned int> Utilities::findVectorMax(const std::vector<double>& v)
{
	assert(v.size() > 0);
	double max;
	unsigned iMax;
	for(unsigned int i = 0; i < v.size(); i++)
	{
		if((i == 0) || (v.at(i) > max))
		{
			iMax = i;
			max = v.at(i);
		}
	}

	return std::pair<double, unsigned int>(max, iMax);
}


std::pair<double, unsigned int> Utilities::findVectorMaxP(std::vector<double>& v)
{
	assert(v.size() > 0);
	double max;
	// unsigned iMax;
	std::vector<unsigned int> iMaxs;
	for(unsigned int i = 0; i < v.size(); i++)
	{
		if((i == 0) || (v.at(i) > max))
		{
			iMaxs.clear();
			iMaxs.push_back(i);
			max = v.at(i);
		}
	}

	unsigned int selectedI;
	if(iMaxs.size() == 1)
	{
		selectedI = iMaxs.at(0);
	}
	else
	{
		selectedI = iMaxs.at(randomNumber(iMaxs.size() - 1));
	}

	return std::pair<double, unsigned int>(max, selectedI);
}

// Converts char to probability of correct genotype, according to Illumina scheme
double Utilities::PhredToPCorrect(unsigned char nucleotideQuality)
{
	if(nucleotideQuality == 0)
	{
		return -1;
	}
	int illuminaPhred = (int)nucleotideQuality - 33;
	assert(illuminaPhred >= 0);

	//char illuminaPhred = nucleotideQuality;
	double log10_pWrong = (double)illuminaPhred / (double)-10;
	double pWrong = exp(log(10) * log10_pWrong);
	assert(pWrong >= 0);
	assert(pWrong <= 1);

	double pRight = 1 - pWrong;
	assert(pRight >= 0);
	assert(pRight <= 1);
	
	return pRight;
}

std::pair<double, unsigned int> Utilities::findVectorMaxP_nonCritical(std::vector<double>& v, unsigned int* thisSeed)
{
	assert(v.size() > 0);
	double max;
	// unsigned iMax;
	std::vector<unsigned int> iMaxs;
	for(unsigned int i = 0; i < v.size(); i++)
	{
		if((i == 0) || (v.at(i) > max))
		{
			iMaxs.clear();
			iMaxs.push_back(i);
			max = v.at(i);
		}
	}

	unsigned int selectedI;
	if(iMaxs.size() == 1)
	{
		selectedI = iMaxs.at(0);
	}
	else
	{
		selectedI = iMaxs.at(randomNumber_nonCritical(iMaxs.size() - 1, thisSeed));
	}

	return std::pair<double, unsigned int>(max, selectedI);
}

vector<string> Utilities::ItoStr(vector<int> i)
{
	vector<string> forReturn;
	for(int I = 0; I < (int)i.size(); I++)
	{
		forReturn.push_back(ItoStr(i.at(I)));
	}
	return forReturn;
}

std::string Utilities::BtoStr(bool b)
{
	std::stringstream sstm;
	sstm << b;
	return sstm.str();
}

std::vector<std::string> Utilities::BtoStr(std::vector<bool> b)
{
	std::vector<string> forReturn;
	for(int I = 0; I < (int)b.size(); I++)
	{
		forReturn.push_back(BtoStr(b.at(I)));
	}
	return forReturn;
}

double Utilities::randomDouble()
{
	assert(RAND_MAX != 0);
	double f = (double)rand() / RAND_MAX;
	assert(f >= 0);
	assert(f <= 1);
	return f;
}


int Utilities::chooseFromVector(vector<double>& v)
{
	//TODO activate at later point
	// srand ( time(NULL) );
	
	assert(RAND_MAX != 0);
	double f = (double)rand() / RAND_MAX;
	double resolution = 1 / RAND_MAX;
	assert(f >= 0);
	assert(f <= 1);
	
	double sum = 0.0;
	for(unsigned int i = 0; i < v.size(); i++)
	{
		sum += v.at(i);
	}
	assert(sum != 0);
	
	int problematic_items = 0;
	double problematic_proportion = 0.0;
	double running_sum = 0.0;
	int selected_i = -1;
	
	assert(v.size() > 0);
	for(unsigned int i = 0; i < v.size(); i++)
	{
		double proportion = v.at(i)/sum;
		running_sum += proportion;
		if(proportion < resolution)
		{
			problematic_items++;
			problematic_proportion += proportion;
		}
		
		if(f <= running_sum)
		{
			selected_i = i;
			break;
		}
		if(i == (v.size() - 1))
		{
			selected_i = i;
			cout << "WARNING!!\n";
			cout << "selected_i: " << selected_i << "\n";
			cout << "f: " << f << "\n";
			cout << "running_sum: " << running_sum << "\n";
			cout << "resolution: " << resolution << "\n";
			cout << "problematic_items: " << problematic_items << "\n";
			cout << "sum: " << sum << "\n";

			double deb_running_sum = 0.00;
			for(unsigned int k = 0; k < v.size(); k++)
			{
				double deb_proportion = v.at(k)/sum;
				cout << "Vector v element " << k << "\n";
				cout << "\tv.at(k): " << v.at(k) << "\n";
				cout << "\tdeb_proportion: "<< deb_proportion << "\n";
				cout << "\tdeb_running_sum before: "<< deb_running_sum << "\n";

				deb_running_sum += deb_proportion;

				cout << "\tdeb_running_sum after: "<< deb_running_sum << "\n";


			}

		}
	}
	
	assert(selected_i != -1);
	assert(problematic_proportion < 0.25);
	
	return selected_i;
}

vector<string> Utilities::split(const string &s, char delim, vector<string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    return elems;
}

vector<string> Utilities::split(string input, string delimiter)
{
	vector<string> output;
	if(input.length() == 0)
	{
		return output;
	}

	if(delimiter == "")
	{
		output.reserve(input.size());
		for(unsigned int i = 0; i < input.length(); i++)
		{
			output.push_back(input.substr(i, 1));
		}
	}
	else
	{
		if(input.find(delimiter) == string::npos)
		{
			output.push_back(input);
		}
		else
		{
			int s = 0;
			int p = input.find(delimiter);

			do {
				output.push_back(input.substr(s, p - s));
				s = p + delimiter.size();
				p = input.find(delimiter, s);
			} while (p != (int)string::npos);
			output.push_back(input.substr(s));
		}
	}

	return output;
}


vector<string> Utilities::split(const string &s, char delim) {
    vector<string> elems;
    return split(s, delim, elems);
}

string Utilities::ItoStr(int i)
{
	std::stringstream sstm;
	sstm << i;
	return sstm.str();
}


string Utilities::PtoStr(void* p)
{
	std::stringstream sstm;
	sstm << p;
	return sstm.str();
}


string Utilities::DtoStr(double d)
{
	std::stringstream sstm;
	sstm << d;
	return sstm.str();
}

string Utilities::join(vector<string> parts, string delim)
{
	if(parts.size() == 0)
		return "";

	string ret = parts.at(0);

	for(unsigned int i = 1; i < parts.size(); i++)
	{
		ret.append(delim);
		ret.append(parts.at(i));
	}

	return ret;
}

double Utilities::proportionN(std::string S)
{
	size_t c = 0;
	for(size_t i = 0; i < S.length(); i++)
	{
		if((S.at(i) == 'N') || (S.at(i) == 'n'))
		{
			c++;
		}
	}

	return ( (double)c / (double) S.length());
}

void Utilities::eraseNL(string& s)
{
	while (!s.empty() && ((s[s.length()-1] == '\r') || (s[s.length()-1] == '\n'))) {
	    s.erase(s.length()-1);
	}
}

long long Utilities::StrtoLongLong(string s)
{
	  stringstream ss(s);
	  long long i;
	  ss >> i;
	  return i;
}

int Utilities::StrtoI(string s)
{
	  stringstream ss(s);
	  int i;
	  ss >> i;
	  return i;
}

double Utilities::StrtoD(string s)
{
	  stringstream ss(s);
	  double d;
	  ss >> d;
	  return d;
}

bool Utilities::StrtoB(string s)
{
	  stringstream ss(s);
	  bool b;
	  ss >> b;
	  return b;
}



std::string Utilities::timestamp()
{
    time_t ltime; /* calendar time */
    ltime=time(NULL); /* get current cal time */
    char* timeCString = asctime( localtime(&ltime) );
    std::string forReturn(timeCString);
    return " [ "+ forReturn.substr(0, forReturn.length() - 1)+" ] ";
}

void Utilities::writeFASTA(std::string file, const std::map<std::string, std::string>& sequences)
{
	std::ofstream FASTAstream;

	FASTAstream.open(file.c_str());
	if(! FASTAstream.is_open())
	{
		throw std::runtime_error("writeFASTA(): Cannot open file "+file);
	}

	for(std::map<std::string, std::string>::const_iterator seqIt = sequences.begin(); seqIt != sequences.end(); seqIt++)
	{
		std::string sequenceID = seqIt->first;
		std::string sequence = seqIt->second;

		FASTAstream << ">" << sequenceID << "\n";

		if(sequence.length() == 0)
		{
			continue;
		}

		size_t line_length = 50;
		size_t chunks = sequence.length() / line_length;
		for(size_t chunkI = 0; chunkI < (chunks + 1); chunkI++)
		{
			size_t chunk_start = chunkI * line_length;
			size_t chunk_stop = (chunkI + 1) * line_length - 1;

			if(chunk_start > (sequence.length() - 1))
			{
				break;
			}

			if(chunk_stop > (sequence.length() - 1))
			{
				chunk_stop = sequence.length() - 1;
			}

			assert(chunk_start <= (sequence.length() - 1));
			int L = chunk_stop - chunk_start + 1;
			assert(L > 0);
			std::string line = sequence.substr(chunk_start, L);
			assert(line.length() > 0);
			assert(line.length() <= line_length);

			FASTAstream << line << "\n";
		}
	}


	FASTAstream.close();

	bool paranoid = true;
	if(paranoid)
	{
		std::map<std::string, std::string> reconstructed_sequences = readFASTA(file, true);
		assert(reconstructed_sequences.size() == sequences.size());
		for(std::map<std::string, std::string>::const_iterator seqIt = sequences.begin(); seqIt != sequences.end(); seqIt++)
		{
			assert(seqIt->second == reconstructed_sequences.at(seqIt->first));
		}
	}
}

bool Utilities::sequence_all_Ns(const std::string S)
{
	for(size_t i = 0; i < S.length(); i++)
	{
		if(!((S.at(i) == 'N') || (S.at(i) == 'n')))
		{
			return false;
		}
	}
	return true;
}


map<string, string> Utilities::readFASTA(std::string file, bool fullIdentifier)
{
	map<string, string> forReturn;

	std::ifstream FASTAstream;

	FASTAstream.open(file.c_str());
	if(! FASTAstream.is_open())
	{
		throw std::runtime_error("readFASTA(): Cannot open file "+file);
	}
	
	while(FASTAstream.good())
	{
		std::string line;
		size_t lineCounter = 0;

		std::string currentSequenceIdentifier;
		while(FASTAstream.good())
		{
			std::getline(FASTAstream, line);
			Utilities::eraseNL(line);

			lineCounter++;

			if(line.substr(0, 1) == ">")
			{
				std::string ident = line.substr(1);

				if(! fullIdentifier)
				{
					for(size_t i = 0; i < ident.size(); i++)
					{
						if(ident.at(i) == ' ')
						{
							ident = ident.substr(0, i);
							break;
						}
					}
				}
				currentSequenceIdentifier = ident;
				assert(forReturn.count(ident) == 0);
			}
			else
			{
				forReturn[currentSequenceIdentifier] += line;
			}
		}
	}

	return forReturn;
}

std::string Utilities::repeatString(std::string s, int repeatNumber)
{
	std::string forReturn;
	assert(repeatNumber >= 0);
	for(int i = 1; i <= repeatNumber; i++)
	{
		forReturn.append(s);
	}
	return forReturn;
}


std::pair<std::string,std::vector<int>> Utilities::modifySequence(std::string sequence, std::vector<int> positionOrigin, double mutationFrequence, double insertionFrequence, int insertionMaxLength, double deletionFrequence, int deletionMaxLength)
{
	assert(sequence.size() == positionOrigin.size());

	assert((mutationFrequence >= 0) && (mutationFrequence <= 1));
	assert((insertionFrequence >= 0) && (insertionFrequence <= 1));
	assert((deletionFrequence >= 0) && (deletionFrequence <= 1));
	assert(insertionMaxLength > 1);
	assert(deletionMaxLength > 1);

	std::pair<std::string,std::vector<int>> forReturn;
	for(int i = 0; i < (int)sequence.length(); i++)
	{
		if(randomDouble() <= mutationFrequence)
		{
			char newNuc = randomNucleotide();
			forReturn.first.push_back(newNuc);
			if(newNuc == sequence.at(i))
			{
				forReturn.second.push_back(positionOrigin.at(i));
			}
			else
			{
				forReturn.second.push_back(-1);
			}
		}
		else if(randomDouble() <= deletionFrequence)
		{
			int deletionLength = randomNumber(deletionMaxLength - 1);
			for(int dI = 0; dI < deletionLength; dI++)
			{
				forReturn.first.push_back('_');
				forReturn.second.push_back(-1);
			}
			i += deletionLength;
		}
		else if(randomDouble() <= insertionFrequence)
		{
			int insertionLength = randomNumber(insertionMaxLength - 1);
			forReturn.first.append(generateRandomSequence(insertionLength));
			for(int dI = 0; dI < insertionLength; dI++)
			{
				forReturn.second.push_back(-1);
			}
		}
		else
		{
			forReturn.first.push_back(sequence.at(i));
			forReturn.second.push_back(positionOrigin.at(i));

		}
	}

	assert(forReturn.first.size() == forReturn.second.size());

	return forReturn;
}

std::string Utilities::generateRandomSequenceWithGaps(int length, double gapFrequency)
{
	std::string forReturn;
	forReturn.resize(length);
	for(int i = 0; i < length; i++)
	{
		if(randomDouble() <= gapFrequency)
		{
			forReturn.at(i) = '_';
		}
		else
		{
			forReturn.at(i) = randomNucleotide();
		}
	}
	return forReturn;
}

std::string Utilities::generateRandomSequence(int length)
{
	std::string forReturn;
	forReturn.resize(length);
	for(int i = 0; i < length; i++)
	{
		forReturn.at(i) = randomNucleotide();
	}
	return forReturn;
}

int Utilities::randomNumber(int max)
{
	assert(max >= 0);
	int n;
	#pragma omp critical
	{
		// n = rand() % (max + 1);
		n = rand_r(&globalRandRSeed) % (max + 1);
	}
	assert((n >= 0) && (n <= max));
	return n;
}

int Utilities::randomNumber_nonCritical(int max, unsigned int* thisSeed)
{
	int n = rand_r(thisSeed) % (max + 1);
	assert((n >= 0) && (n <= max));
	return n;
}

char Utilities::randomNucleotide()
{
	char nucleotides[4] = {'A', 'C', 'G', 'T'};
	int n = rand() % 4;
	assert((n >= 0) && (n <= 3));
	return nucleotides[n];
}


void Utilities::writeStatus(std::string statusFile, int status)
{
	std::ofstream statusStream;
	statusStream.open(statusFile.c_str(), std::ios::out);
	if(! statusStream.is_open())
	{
		throw std::runtime_error("Cannot open file for writing: "+statusFile);
	}
	statusStream << status << "\n";
	statusStream.close();
}


int Utilities::readStatus(std::string statusFile)
{
	std::ifstream statusStream;
	statusStream.open(statusFile.c_str(), std::ios::in);
	if(! statusStream.is_open())
	{
		return 0;
	}
	assert(statusStream.good());
	std::string line;
	getline (statusStream, line);
	Utilities::eraseNL(line);
	if(line.length() == 0)
	{
		return 0;
	}
	int status = Utilities::StrtoI(line);
	assert(status >= 0);

	return status;
}

std::vector<double> Utilities::normalize_log_vector(const std::vector<double>& v)
{
	std::vector<double> forReturn = v;
	std::pair<double, unsigned int> max_v = findVectorMax(forReturn);
	for(size_t i = 0; i < forReturn.size(); i++)
	{
		forReturn.at(i) = exp(forReturn.at(i) - max_v.first);
		assert(forReturn.at(i) >= 0);
		assert(forReturn.at(i) <= 1);
	}
	normalize_vector(forReturn);
	return forReturn;
}

void Utilities::normalize_vector(std::vector<double>& v)
{
	double S = 0;
	for(unsigned int i = 0; i < v.size(); i++)
	{
		S += v.at(i);
	}

	for(unsigned int i = 0; i < v.size(); i++)
	{
		v.at(i) = v.at(i) / S;
	}
}

void Utilities::check_map_is_normalized(std::map<char, double> m)
{
	double sum = 0;
	for(std::map<char, double>::iterator cIt = m.begin(); cIt != m.end(); cIt++)
	{
		assert(cIt->second >= 0);
		sum += cIt->second;
	}
	assert(abs(sum - 1) < 1e-5);
}

std::map<char, double> Utilities::normalize_map(std::map<char, double> m)
{
	std::map<char, double> forReturn;
	double sum = 0;
	for(std::map<char, double>::iterator cIt = m.begin(); cIt != m.end(); cIt++)
	{
		assert(cIt->second >= 0);
		sum += cIt->second;
	}
	assert(sum > 0);
	for(std::map<char, double>::iterator cIt = m.begin(); cIt != m.end(); cIt++)
	{
		forReturn[cIt->first] = cIt->second/sum;
	}
	return forReturn;
}


std::string Utilities::choose_uniformly_from_vector(std::vector<std::string> v)
{
	assert(v.size() > 0);
	int R = rand();
	assert(R >= 0);
	R = (R % v.size());
	return v.at(R);
}

int Utilities::choose_from_nonnormalized_map(std::map<int, int> m)
{
	//TODO activate at later point
	// srand ( time(NULL) );

	assert(m.size());
	assert(RAND_MAX != 0);
	double f = (double)rand() / RAND_MAX;
	// double resolution = 1 / RAND_MAX;
	assert(f >= 0);
	assert(f <= 1);

	int forReturn;

	// std::cout << "Sum!\n";
	long long m_sum = 0;
	for(auto e : m)
	{
		assert(e.second >= 0);
		m_sum += e.second;
		// std::cout << e.first << "\t" << e.second << "\t" << m_sum << "\n";
	}
	assert(m_sum > 0);
	
	// std::cout << "Normalize.\n" << std::flush;
	std::map<int, double> normalized_m;
	for(auto e : m)
	{
		// std::cout << e.first << "\t" << e.second << "\t" << m_sum << "\n";
		
		double normalized_value = (double)e.second/(double)m_sum;
		if(!((normalized_value >= 0) && (normalized_value <= 1)))
		{
			std::cerr << "normalized_value" << ": " << normalized_value << "\n";
			std::cerr << "m_sum" << ": " << m_sum << "\n";
			std::cerr << "e.second" << ": " << e.second << "\n";
			std::cerr << std::flush;
		}
		assert((normalized_value >= 0) && (normalized_value <= 1));		
		normalized_m[e.first] = normalized_value;
		assert((normalized_m[e.first] >= 0) && (normalized_m[e.first] <= 1));
	}

	double current_sum = 0;
	for(std::map<int, double>::iterator mIt = normalized_m.begin(); mIt != normalized_m.end(); mIt++)
	{
		current_sum += mIt->second;
		if(f <= current_sum)
		{
			forReturn = mIt->first;
			current_sum = -1;
			break;
		}
	}
	if(current_sum != -1)
	{
		std::cerr << "current_sum" << ": " << current_sum << "\n";
		std::cerr << std::flush;
	}
	assert(current_sum == -1);

	return forReturn;
}

char Utilities::choose_from_normalized_map(std::map<char, double> m)
{
	//TODO activate at later point
	// srand ( time(NULL) );

	assert(RAND_MAX != 0);
	double f = (double)rand() / RAND_MAX;
	// double resolution = 1 / RAND_MAX;
	assert(f >= 0);
	assert(f <= 1);

	char forReturn;

	double current_sum = 0;
	for(std::map<char, double>::iterator mIt = m.begin(); mIt != m.end(); mIt++)
	{
		current_sum += mIt->second;
		if(f <= current_sum)
		{
			forReturn = mIt->first;
			current_sum = -1;
			break;
		}
	}
	assert(current_sum == -1);

	return forReturn;
}


char Utilities::choose_from_normalized_map(std::map<char, double> m, boost::mt19937& rng)
{
	boost::random::uniform_01<boost::mt19937&> f_gen(rng);
	double f = f_gen();
	assert(f >= 0);
	assert(f <= 1);

	char forReturn;

	double current_sum = 0;
	for(std::map<char, double>::iterator mIt = m.begin(); mIt != m.end(); mIt++)
	{
		current_sum += mIt->second;
		if(f <= current_sum)
		{
			forReturn = mIt->first;
			current_sum = -1;
			break;
		}
	}
	assert(current_sum = -1);

	return forReturn;
}

std::string Utilities::seq_reverse_complement(std::string sequence)
{
	int length = sequence.size();
	std::string forReturn;
	forReturn.resize(length);
    for(int k=0; k < length; k++)
    {
        forReturn[k] = reverse_char_nucleotide(sequence.at(length-k-1));
    }
    return forReturn;
}


char Utilities::reverse_char_nucleotide(char c)
{
    switch (c)
    {
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		case 'N':
			return 'N';
		case 'a':
			return 't';
		case 'c':
			return 'g';
		case 'g':
			return 'c';
		case 't':
			return 'a';
		case 'n':
			return 'n';
		case '_':
			return '_';
		case '*':
			return '*';
		default:
			std::string errorString = "Utilities::reverse_char_nucleotide: nucleotide not existing!";
			errorString.push_back(c);
			throw std::runtime_error(errorString);
    }
}



bool Utilities::oneBernoulliTrial(double p)
{
	//TODO activate at later point
	// srand ( time(NULL) );

	assert(RAND_MAX != 0);
	double f = (double)rand() / RAND_MAX;
	// double resolution = 1 / RAND_MAX;
	assert(f >= 0);
	assert(f <= 1);

	if(f <= p)
	{
		return true;
	}
	else
	{
		return false;
	}
}


bool Utilities::oneBernoulliTrial(double p, boost::mt19937& rng)
{
	boost::random::uniform_01<boost::mt19937&> f_gen(rng);
	double f = f_gen();
	assert(f >= 0);
	assert(f <= 1);

	if(f <= p)
	{
		return true;
	}
	else
	{
		return false;
	}
}


char Utilities::randomNucleotide(boost::mt19937& rng)
{
	boost::random::uniform_int_distribution<> nucleotide_gen (0,3);
	char nucleotides[4] = {'A', 'C', 'G', 'T'};
	int n = nucleotide_gen(rng);
	assert((n >= 0) && (n <= 3));
	return nucleotides[n];
}



std::string Utilities::removeFROM(std::string readID)
{
	if(readID.find(":FROM:") != std::string::npos)
	{
		size_t cutFrom = readID.find(":FROM:");
		return readID.substr(0, cutFrom);
	}
	else
	{
		return readID;
	}
}

void Utilities::printToFile(std::string file, std::string what)
{
	std::ofstream stream;
	stream.open(file.c_str());
	assert(stream.is_open());
	stream << what;
	stream.close();
}


std::vector<std::string> Utilities::getAllLines(std::string file)
{
	std::vector<std::string> forReturn;
    assert(Utilities::fileExists(file));
	std::ifstream stream;
	stream.open(file.c_str());
	assert(stream.is_open());
	std::string line;
	while(stream.good())
	{
		std::getline(stream, line);
		Utilities::eraseNL(line);
		if(line.length())
		{
			forReturn.push_back(line);
		}
	}

	stream.close();
	return forReturn;
}


std::string Utilities::getFirstLine(std::string file)
{
	if(!Utilities::fileExists(file))
	{
		std::cerr << "File " << file << " not existing.\n" << std::flush;
	}
    assert(Utilities::fileExists(file));
	std::ifstream stream;
	stream.open(file.c_str());
	assert(stream.is_open());
	std::string line;
	assert(stream.good());
	std::getline(stream, line);
	Utilities::eraseNL(line);
	stream.close();
	return line;
}

std::string Utilities::removeGaps(std::string in)
{
	std::string out;
	out.reserve(in.size());
	for(size_t i = 0; i < in.size(); i++)
	{
		if(in.at(i) != '_')
		{
			out.push_back(in.at(i));
		}
	}
	return out;
}

std::vector<std::string> Utilities::filesInDirectory(std::string path)
{
	std::vector<std::string> forReturn;

    struct dirent *pDirent;
    DIR *pDir;

    pDir = opendir (path.c_str());
    if (pDir == NULL) {
        throw std::runtime_error("Cannot open directory " + path);
    }

    while ((pDirent = readdir(pDir)) != NULL) {
    	std::string fileName(pDirent->d_name);

    	if((fileName == ".") || (fileName == ".."))
    	{
    		continue;
    	}

    	fileName = path + "/" + fileName;

    	// std::cout << "\t" << fileName << "\n";

    	forReturn.push_back(fileName);
    }
    closedir (pDir);

    return forReturn;
}

double Utilities::logAvg(double a, double b)
{
	if(a > b)
	{
		return(log(0.5) + (log(1 + exp(b - a)) + a));

	} else
	{
		return(log(0.5) + (log(1 + exp(a - b)) + b));
	}
}


std::vector<std::string> Utilities::partitionStringIntokMers(std::string str, int k)
{
	std::vector<std::string> forReturn;
	if((int)str.length() >= k)
	{
		for(int i = 0; i <= (str.length() - k); i++)
		{
			std::string kMer = str.substr(i, k);
			assert((int)kMer.length() == k);
			forReturn.push_back(kMer);
		}
	}
	return forReturn;
}

std::vector<std::string> Utilities::get_map_keys_sorted_by_value(const std::map<std::string, int>& m)
{
	std::vector<std::string> keys;
	for(auto e : m)
	{
		keys.push_back(e.first);
	}
	std::sort( keys.begin( ), keys.end( ), [&]( const std::string& lhs, const std::string& rhs )
	{
	   return (m.at(lhs) < m.at(rhs));
	});

	std::reverse(keys.begin(), keys.end());

	if(keys.size() > 1)
	{
		assert(m.at(keys.at(0)) >= m.at(keys.at(1)));
	}
	return keys;
}

